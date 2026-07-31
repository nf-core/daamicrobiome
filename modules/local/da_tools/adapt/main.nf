process ADAPT {
  tag { "ADAPT:${rep_id}" }

  conda "${moduleDir}/environment.yml"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'oras://community.wave.seqera.io/library/bioconductor-adapt_bioconductor-phyloseq_r-optparse:bd485df825838000' :
      'community.wave.seqera.io/library/bioconductor-adapt_bioconductor-phyloseq_r-optparse:7fbd22e8854667d5' }"

  input:
  tuple val(rep_id), path(rds_file)

  output:
  tuple val(rep_id), path("adapt_da.tsv")
  path("versions.yml"), emit: versions

  script:
  // optional flags from params (fallbacks if not set)
  def condFlag = params.metadata_condition ? "--condition '${params.metadata_condition}'" : ""
  def baseFlag = params.metadata_base      ? "--base '${params.metadata_base}'"         : ""
  def confounderFlag = params.metadata_confounder ? "--confounder '${params.metadata_confounder}'" : ""
  """
  set -euo pipefail
  echo "[ADAPT] Input: '${rds_file}'" >&2

  run_adapt.R --input "${rds_file}" --output_dir "." ${condFlag} ${baseFlag} ${confounderFlag}

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      ADAPT: \$(Rscript -e "cat(as.character(packageVersion('ADAPT')))")
      R: \$(Rscript -e "cat(R.version.string)")
  END_VERSIONS
  """
}