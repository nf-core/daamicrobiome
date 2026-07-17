// Runs all enabled DA tools on real (non-simulated) data.
// Uses params.condition / params.base_level for group comparison.
// Emits a unified stream of (rep_id, tool_id, result_path) tuples.

include { ADAPT_REAL         } from '../../../modules/local/da_tools_real/adapt_real/main'
include { CORNCOB_REAL       } from '../../../modules/local/da_tools_real/corncob_real/main'
include { LINDA_REAL         } from '../../../modules/local/da_tools_real/linda_real/main'
include { LOCOM_REAL         } from '../../../modules/local/da_tools_real/locom_real/main'
include { MAASLIN2_REAL      } from '../../../modules/local/da_tools_real/maaslin2_real/main'
include { METAGENOMESEQ_REAL } from '../../../modules/local/da_tools_real/metagenomeseq_real/main'

workflow DIFF_ABUNDANCE_REAL {

    take:
    rep_rds_ch   // tuples: (rep_id, rds_path)

    main:
    results = channel.empty()

    if (params.tool_adapt) {
        ch = ADAPT_REAL(rep_rds_ch)
        ch = ch.map { rid, p -> tuple(rid, 'ADAPT_real', p) }
        results = results.mix(ch)
    }

    if (params.tool_corncob) {
        ch = CORNCOB_REAL(rep_rds_ch)
        ch = ch.map { rid, p -> tuple(rid, 'CORNCOB_real', p) }
        results = results.mix(ch)
    }

    if (params.tool_linda) {
        ch = LINDA_REAL(rep_rds_ch)
        ch = ch.map { rid, p -> tuple(rid, 'LINDA_real', p) }
        results = results.mix(ch)
    }

    if (params.tool_locom) {
        def loc = LOCOM_REAL(rep_rds_ch)
        ch = loc.da.map { rid, p -> tuple(rid, 'LOCOM_real', p) }
        results = results.mix(ch)
    }

    if (params.tool_maaslin2) {
        ch = MAASLIN2_REAL(rep_rds_ch)
        ch = ch.map { rid, p -> tuple(rid, 'MAASLIN2_real', p) }
        results = results.mix(ch)
    }

    if (params.tool_metagenomeseq) {
        METAGENOMESEQ_REAL(rep_rds_ch)
        def ch_da = METAGENOMESEQ_REAL.out.da
            .map { rid, p -> tuple(rid, 'METAGENOMESEQ_real', p) }
        results = results.mix(ch_da)
    }

    if (![params.tool_adapt, params.tool_corncob, params.tool_linda,
          params.tool_locom, params.tool_maaslin2, params.tool_metagenomeseq].any { it }) {
        error "No DA tools enabled. Set one or more params.tool_* = true."
    }

    emit:
    da_results = results
}
