/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// SUBWORKFLOWS
include { DIFF_ABUNDANCE as DIFFABUNDANCE }       from '../subworkflows/local/differential_abundance/main'
include { DIFF_ABUNDANCE_REAL as DIFFABUNDANCE_REAL } from '../subworkflows/local/differentiaL_abundance_real/main'
include { SIMULATION }                              from '../subworkflows/local/simulation/main'
include { DA_SCORING_WF as DA_SCORING }             from '../subworkflows/local/da_scoring/main'
include { APPLY_WEIGHTED_CONSENSUS_WF }                 from '../subworkflows/local/apply_weighted_consensus/main'

// MODULES
include { EXTRACT_CONTROL } from '../modules/local/extract_control/main'
include { K_INTERSECTION }  from '../modules/local/k_intersection/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DAAMICROBIOME {

    log.info """    ─────────────────────────────────────────────────────────
     nf-core/daamicrobiome - Microbiome DA Consensus Pipeline
            runName   : ${workflow.runName}
            profile   : ${workflow.profile}
            version   : ${workflow.manifest.version ?: 'dev'}
    ─────────────────────────────────────────────────────────
    """.stripIndent()

    def INDENT = '        '
    log.info " Pipeline parameters:"
    def keys  = params.keySet().sort()
    def width = keys*.size().max() ?: 0
    keys.each { k ->
        def v = params[k]
        if (v instanceof List || v instanceof Map) {
            v = groovy.json.JsonOutput.toJson(v)
        }
        log.info "${INDENT}${k.padRight(width)}   : ${v} "
    }
    log.info "─────────────────────────────────────────────────────────\n"

    // Validate common inputs
    if (!params.input) error "params.input must be set (path to the full phyloseq RDS)."

    // =========================================================================
    // PATH A: No simulation -- real-only k-intersection consensus
    // =========================================================================
    if (!params.simulate) {

        log.info " MODE: Real-only (no simulation). Running DA tools and k-intersection."

        def real_ch = channel
            .fromPath(params.input, checkIfExists: true)
            .map { rds -> tuple(rds.baseName, rds) }

        // Run DA tools on real data
        def da_real = DIFFABUNDANCE_REAL(real_ch)

        // Collect DA result paths grouped by rep_id for k-intersection
        def per_rep = da_real.da_results
            .map { rid, tool, p -> tuple(rid, p) }
            .groupTuple()

        // Build XLSX with k=1..n_tools intersection sheets
        K_INTERSECTION(per_rep)
    }

    // =========================================================================
    // PATH B: Simulation-trained weighted consensus
    // =========================================================================
    else {

        log.info " MODE: Simulation-trained weighted consensus."

        // -----------------------------------------------------------------
        // Step 1: Determine control-only phyloseq for simulation
        // -----------------------------------------------------------------
        def control_rds_ch

        if (params.input_control) {
            log.info " Using provided control-only phyloseq: ${params.input_control}"
            control_rds_ch = channel.fromPath(params.input_control, checkIfExists: true).first()
        } else {
            log.info " Extracting control samples from input (condition='${params.condition}', base='${params.base_level}')"
            def full_rds = channel.fromPath(params.input, checkIfExists: true).first()
            control_rds_ch = EXTRACT_CONTROL(full_rds).control_rds
        }

        // -----------------------------------------------------------------
        // Step 2: Simulate datasets from control-only phyloseq
        // -----------------------------------------------------------------
        def sim = SIMULATION(control_rds_ch)

        // -----------------------------------------------------------------
        // Step 3: Run DA tools on simulated replicates
        // -----------------------------------------------------------------
        def da_sim = DIFFABUNDANCE(sim.rds)

        // -----------------------------------------------------------------
        // Step 4: Score tools using weighted consensus scoring
        // -----------------------------------------------------------------
        def sim_root_dir_ch = channel.value(file("${params.outdir}/da_tools"))

        def scoring = DA_SCORING(da_sim.da_results, sim_root_dir_ch, sim.truth_dir)
        def scoring_dir_ch = scoring.scoring_dir

        // -----------------------------------------------------------------
        // Step 5: Run DA tools on full real dataset
        // -----------------------------------------------------------------
        def real_ch = channel
            .fromPath(params.input, checkIfExists: true)
            .map { rds -> tuple("REAL_${rds.baseName}", rds) }

        def da_real = DIFFABUNDANCE_REAL(real_ch)

        // Barrier: wait for all real DA tasks to complete
        def done_real_ch = da_real.da_results.collect().map { "REAL_DA_DONE" }

        // Path to the published real DA results
        def real_da_root_ch = real_ch.map { rid, rds ->
            file("${params.outdir}/da_tools_real/${rid}")
        }.first()

        // -----------------------------------------------------------------
        // Step 6: Apply weighted consensus threshold to real results
        // -----------------------------------------------------------------
        APPLY_WEIGHTED_CONSENSUS_WF(done_real_ch, scoring_dir_ch, real_da_root_ch)
    }
}
