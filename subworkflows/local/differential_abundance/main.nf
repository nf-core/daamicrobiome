// Runs all enabled DA tools on simulated replicates.
// Each tool is conditionally executed based on params.tool_* toggles.
// Emits a unified stream of (rep_id, tool_id, result_path) tuples.

include { ADAPT         } from '../../../modules/local/da_tools/adapt/main'
include { CORNCOB       } from '../../../modules/local/da_tools/corncob/main'
include { LINDA         } from '../../../modules/local/da_tools/linda/main'
include { LOCOM         } from '../../../modules/local/da_tools/locom/main'
include { MAASLIN2      } from '../../../modules/local/da_tools/maaslin2/main'
include { METAGENOMESEQ } from '../../../modules/local/da_tools/metagenomeseq/main'

workflow DIFF_ABUNDANCE {

    take:
    rep_rds_ch   // tuples: (rep_id, rds_path)

    main:
    results = Channel.empty()

    if (params.tool_adapt) {
        ch = ADAPT(rep_rds_ch)
        ch = ch.map { rid, p -> tuple(rid, 'ADAPT', p) }
        results = results.mix(ch)
    }

    if (params.tool_corncob) {
        ch = CORNCOB(rep_rds_ch)
        ch = ch.map { rid, p -> tuple(rid, 'CORNCOB', p) }
        results = results.mix(ch)
    }

    if (params.tool_linda) {
        ch = LINDA(rep_rds_ch)
        ch = ch.map { rid, p -> tuple(rid, 'LINDA', p) }
        results = results.mix(ch)
    }

    if (params.tool_locom) {
        def loc = LOCOM(rep_rds_ch)
        ch = loc.da.map { rid, p -> tuple(rid, 'LOCOM', p) }
        results = results.mix(ch)
    }

    if (params.tool_maaslin2) {
        ch = MAASLIN2(rep_rds_ch)
        ch = ch.map { rid, p -> tuple(rid, 'MAASLIN2', p) }
        results = results.mix(ch)
    }

    if (params.tool_metagenomeseq) {
        METAGENOMESEQ(rep_rds_ch)
        def ch_da = METAGENOMESEQ.out.da
            .map { rid, p -> tuple(rid, 'METAGENOMESEQ', p) }
        results = results.mix(ch_da)
    }

    if (![params.tool_adapt, params.tool_corncob, params.tool_linda,
          params.tool_locom, params.tool_maaslin2, params.tool_metagenomeseq].any { it }) {
        error "No DA tools enabled. Set one or more params.tool_* = true."
    }

    emit:
    da_results = results
}
