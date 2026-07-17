process MAASLIN2_REAL {
  tag {"MAASLIN2_REAL:${rep_id}"}

  input:
  tuple val(rep_id), path(rds_file)

  output:
   tuple val(rep_id), path("maaslin2_da.tsv")
   path("versions.yml"), emit: versions

  script:
  // optional flags from params (fallbacks if not set)
  def condFlag = params.condition  ? "--condition '${params.condition}'" : ""
  def baseFlag = params.base_level ? "--base '${params.base_level}'"    : ""
  def conf = params.confounder
  def conf_str = (conf instanceof List) ? conf.join(',') : (conf ? conf.toString() : "")
  conf_str = conf_str.trim()
  def confounderFlag = conf_str ? "--confounder '${conf_str}'" : ""
  def bc = params.bc_confounder
  def bc_str = (bc instanceof List) ? bc[0].toString() : (bc ? bc.toString() : "")
  bc_str = bc_str.trim()
  def bc_confounderFlag = bc_str ? "--bc_confounder '${bc_str}'" : ""

  """
  set -euo pipefail
  echo "[MAASLIN2_REAL] Input: '${rds_file}'" >&2

  run_maaslin2.R --input "${rds_file}" --output_dir "." ${condFlag} ${baseFlag} ${confounderFlag} ${bc_confounderFlag}

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      Maaslin2: \$(Rscript -e "cat(as.character(packageVersion('Maaslin2')))")
      R: \$(Rscript -e "cat(R.version.string)")
  END_VERSIONS
  """
}
