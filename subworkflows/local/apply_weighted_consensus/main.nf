include { APPLY_WEIGHTED_CONSENSUS } from '../../../modules/local/apply_weighted_consensus/main'

workflow APPLY_WEIGHTED_CONSENSUS_WF {

  take:
    done_real_ch
    scoring_dir_ch
    real_da_root_ch

  main:
    result = APPLY_WEIGHTED_CONSENSUS(done_real_ch, scoring_dir_ch, real_da_root_ch)

  emit:
    consensus = result.consensus_dir
}
