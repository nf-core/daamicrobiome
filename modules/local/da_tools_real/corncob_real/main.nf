process CORNCOB_REAL {
  tag {"CORNCOB_REAL:${rep_id}"}

  input:
  tuple val(rep_id), path(rds_file)

  output:
   tuple val(rep_id), path("corncob_da.tsv")
   path("versions.yml"), emit: versions

  script:
    // optional flags from params (fallbacks if not set)
  def condFlag = params.condition  ? "--condition '${params.condition}'" : ""
  def baseFlag = params.base_level ? "--base '${params.base_level}'"    : ""
  def conf = params.confounder
  def conf_str = (conf instanceof List) ? conf.join(',') : (conf ? conf.toString() : "")
  conf_str = conf_str.trim()
  def confounderFlag = conf_str ? "--confounder '${conf_str}'" : ""
  """
  set -euo pipefail
  echo "[CORNCOB_REAL] Input: '${rds_file}'" >&2

  run_corncob.R --input "${rds_file}" --output_dir "." ${condFlag} ${baseFlag} ${confounderFlag}

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      corncob: \$(Rscript -e "cat(as.character(packageVersion('corncob')))")
      R: \$(Rscript -e "cat(R.version.string)")
  END_VERSIONS
  """
}