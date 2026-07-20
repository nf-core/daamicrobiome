process APPLY_WEIGHTED_CONSENSUS {
  tag { "APPLY_WEIGHTED_CONSENSUS" }

  input:
  val(done_real)
  path(scoring_dir)
  path(da_dir)

  output:
  path("weighted_consensus_out"), emit: consensus_dir
  path("versions.yml"), emit: versions

  script:
  """
  set -euo pipefail

  mkdir -p weighted_consensus_out

  Rscript ${projectDir}/bin/apply_weighted_consensus.R \
    --scoring_dir "${scoring_dir}" \
    --da_dir      "${da_dir}" \
    --outdir      "weighted_consensus_out" \
    --alpha       "${params.eval_alpha ?: 0.05}" \
    --lfc_min     "${params.consensus_lfc_min ?: 0}"

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      R: \$(Rscript -e "cat(R.version.string)")
  END_VERSIONS
  """
}
