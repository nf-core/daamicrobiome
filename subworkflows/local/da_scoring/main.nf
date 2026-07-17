include { DA_SCORING } from '../../../modules/local/da_scoring/main'

workflow DA_SCORING_WF {

  take:
  da_results_ch
  sim_root_dir_ch
  truth_root_dir_ch

  main:
  done_ch = da_results_ch.count().map { n -> "DA_DONE_${n}" }

  def tools = []
  if (params.tool_adapt)         tools << 'adapt'
  if (params.tool_corncob)       tools << 'corncob'
  if (params.tool_linda)         tools << 'linda'
  if (params.tool_locom)         tools << 'locom'
  if (params.tool_maaslin2)      tools << 'maaslin2'
  if (params.tool_metagenomeseq) tools << 'metagenomeseq'
  def tools_csv = tools.join(',')

  scored = DA_SCORING(
    done_ch,
    sim_root_dir_ch,
    truth_root_dir_ch,
    channel.value(tools_csv)
  )

  emit:
  scoring_dir = scored.scoring_dir
}
