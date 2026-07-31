process EXTRACT_CONTROL {
  tag { "EXTRACT_CONTROL" }

  conda "${moduleDir}/environment.yml"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'oras://community.wave.seqera.io/library/bioconductor-phyloseq_r-optparse:0dd0584c7eb781d3' :
      'community.wave.seqera.io/library/bioconductor-phyloseq_r-optparse:3b069fdafb1cf2ae' }"

  input:
  path(full_phyloseq)

  output:
  path("control_phyloseq.RDS"), emit: control_rds
  path("versions.yml"), emit: versions

  script:
  """
  set -euo pipefail

  Rscript ${projectDir}/bin/extract_control.R \
    --input    "${full_phyloseq}" \
    --output   "control_phyloseq.RDS" \
    --condition "${params.condition}" \
    --reference "${params.base_level}"

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      phyloseq: \$(Rscript -e "cat(as.character(packageVersion('phyloseq')))")
      R: \$(Rscript -e "cat(R.version.string)")
  END_VERSIONS
  """
}
