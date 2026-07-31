process METAGENOMESEQ {
  tag "METAGENOMESEQ:${rep_id}"

  conda "${moduleDir}/environment.yml"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'oras://community.wave.seqera.io/library/bioconductor-metagenomeseq_bioconductor-phyloseq_r-optparse:52e10cf8b0086f95' :
      'community.wave.seqera.io/library/bioconductor-metagenomeseq_bioconductor-phyloseq_r-optparse:f1a7e0d8232664ef' }"

  input:
  tuple val(rep_id), path(rds_file)

  output:
  tuple val(rep_id), path("metagenomeseq_da.tsv"), emit: da
  tuple val(rep_id), path("metagenomeseq_status.tsv"), emit: status, optional: true
  path("versions.yml"), emit: versions

  script:
  def condFlag = params.metadata_condition ? "--condition '${params.metadata_condition}'" : ""
  def baseFlag = params.metadata_base      ? "--base '${params.metadata_base}'"         : ""
  def confounderFlag = params.metadata_confounder ? "--confounder '${params.metadata_confounder}'" : ""
  def contrastModeFlag = params.metagenomeseq_contrast_mode ? "--contrast_mode '${params.metagenomeseq_contrast_mode}'" : ""
  def contrastFlag     = params.metagenomeseq_contrast      ? "--contrast '${params.metagenomeseq_contrast}'"           : ""

  """
  set -euo pipefail
  echo "[METAGENOMESEQ] Input: '${rds_file}'" >&2

  Rscript "${projectDir}/bin/run_metagenomeseq.R" --input "${rds_file}" --output_dir "." \
    ${condFlag} ${baseFlag} ${confounderFlag} ${contrastModeFlag} ${contrastFlag}

  if [[ -f metagenomeseq_status.tsv ]]; then
    awk -F'\\t' 'NR==2 && \$1==0 {print "[METAGENOMESEQ] WARNING: "\$3" (stage="\$2")" > "/dev/stderr"}' metagenomeseq_status.tsv
  fi

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      metagenomeSeq: \$(Rscript -e "cat(as.character(packageVersion('metagenomeSeq')))")
      R: \$(Rscript -e "cat(R.version.string)")
  END_VERSIONS
  """
}
