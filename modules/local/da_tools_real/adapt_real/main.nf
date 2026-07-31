process ADAPT_REAL {
  tag { "ADAPT_REAL:${rep_id}" }

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
  def condFlag = params.condition  ? "--condition '${params.condition}'" : ""
  def baseFlag = params.base_level ? "--base '${params.base_level}'"    : ""
  def conf = params.confounder
  def conf_str = (conf instanceof List) ? conf.join(',') : (conf ? conf.toString() : "")
  conf_str = conf_str.trim()
  def confounderFlag = conf_str ? "--confounder '${conf_str}'" : ""

  """
  set -euo pipefail
  echo "[ADAPT_REAL] Input: '${rds_file}'" >&2

  run_adapt.R --input "${rds_file}" --output_dir "." ${condFlag} ${baseFlag} ${confounderFlag}

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      ADAPT: \$(Rscript -e "cat(as.character(packageVersion('ADAPT')))")
      R: \$(Rscript -e "cat(R.version.string)")
  END_VERSIONS
  """
}