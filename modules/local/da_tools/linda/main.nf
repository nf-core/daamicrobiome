process LINDA {
  tag "LINDA:${rep_id}"

  conda "${moduleDir}/environment.yml"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'oras://community.wave.seqera.io/library/bioconductor-phyloseq_r-microbiomestat_r-optparse:495b95c5f31d5ae2' :
      'community.wave.seqera.io/library/bioconductor-phyloseq_r-microbiomestat_r-optparse:96e50fbd93fbbf1a' }"

  input:
  tuple val(rep_id), path(rds_file)

  output:
  tuple val(rep_id), path("linda_da.tsv")
  path("versions.yml"), emit: versions

  script:
  // optional flags from params (fallbacks if not set)
  def condFlag = params.metadata_condition ? "--condition '${params.metadata_condition}'" : ""
  def baseFlag = params.metadata_base      ? "--base '${params.metadata_base}'"         : ""
  def confounderFlag = params.metadata_confounder ? "--confounder '${params.metadata_confounder}'" : ""

  """
  set -euo pipefail
  echo "[LINDA] Input: '${rds_file}'" >&2

  run_linda.R --input "${rds_file}" --output_dir "." ${condFlag} ${baseFlag} ${confounderFlag}

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      MicrobiomeStat: \$(Rscript -e "cat(as.character(packageVersion('MicrobiomeStat')))")
      R: \$(Rscript -e "cat(R.version.string)")
  END_VERSIONS
  """
}
