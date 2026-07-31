process MAASLIN2 {
  tag {"MAASLIN2:${rep_id}"}

  conda "${moduleDir}/environment.yml"

  input:
  tuple val(rep_id), path(rds_file)

  output:
   tuple val(rep_id), path("maaslin2_da.tsv")
   path("versions.yml"), emit: versions

  script:
  // optional flags from params (fallbacks if not set)
  def condFlag = params.metadata_condition ? "--condition '${params.metadata_condition}'" : ""
  def baseFlag = params.metadata_base      ? "--base '${params.metadata_base}'"         : ""
  def confounderFlag = params.metadata_confounder ? "--confounder '${params.metadata_confounder}'" : ""
  def bc_confounderFlag = params.metadata_bc_confounder ? "--bc_confounder '${params.metadata_bc_confounder}'" : ""
  """
  set -euo pipefail
  echo "[MAASLIN2] Input: '${rds_file}'" >&2

  run_maaslin2.R --input "${rds_file}" --output_dir "." ${condFlag} ${baseFlag} ${confounderFlag} ${bc_confounderFlag}

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      Maaslin2: \$(Rscript -e "cat(as.character(packageVersion('Maaslin2')))")
      R: \$(Rscript -e "cat(R.version.string)")
  END_VERSIONS
  """
}
