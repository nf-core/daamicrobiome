process CORNCOB {
  tag {"CORNCOB:${rep_id}"}

  input:
  tuple val(rep_id), path(rds_file)

  output:
   tuple val(rep_id), path("corncob_da.tsv")
   path("versions.yml"), emit: versions

  script:
    // optional flags from params (fallbacks if not set)
  def condFlag = params.metadata_condition ? "--condition '${params.metadata_condition}'" : ""
  def baseFlag = params.metadata_base      ? "--base '${params.metadata_base}'"         : ""
  def confounderFlag = params.metadata_confounder ? "--confounder '${params.metadata_confounder}'" : ""
  """
  set -euo pipefail
  echo "[CORNCOB] Input: '${rds_file}'" >&2

  run_corncob.R --input "${rds_file}" --output_dir "." ${condFlag} ${baseFlag} ${confounderFlag}

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      corncob: \$(Rscript -e "cat(as.character(packageVersion('corncob')))")
      R: \$(Rscript -e "cat(R.version.string)")
  END_VERSIONS
  """
}
