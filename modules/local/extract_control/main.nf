process EXTRACT_CONTROL {
  tag { "EXTRACT_CONTROL" }

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
