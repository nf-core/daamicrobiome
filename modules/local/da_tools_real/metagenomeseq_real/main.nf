process METAGENOMESEQ_REAL {
  tag "METAGENOMESEQ_REAL:${rep_id}"

  input:
  tuple val(rep_id), path(rds_file)

  output:
  // tuple val(rep_id), path("metagenomeseq_da.tsv")
  tuple val(rep_id), path("metagenomeseq_da.tsv"), emit: da
  tuple val(rep_id), path("metagenomeseq_status.tsv"), emit: status, optional: true
  path("versions.yml"), emit: versions

  script:
  // optional flags from params (fallbacks if not set)
  def condFlag = params.condition  ? "--condition '${params.condition}'" : ""
  def baseFlag = params.base_level ? "--base '${params.base_level}'"    : ""

  def conf = params.confounder
  def conf_str = (conf instanceof List) ? conf.join(',') : (conf ? conf.toString() : "")
  conf_str = conf_str.trim()
  def confounderFlag = conf_str ? "--confounder '${conf_str}'" : ""

  // metagenomeseq-specific
  def contrastVal  = params.get('metagenomeseq_contrast', "")
  def contrastMode = params.get('metagenomeseq_contrast_mode', "")

  def contrastFlag     = (contrastVal && contrastVal.toString().trim()) ? "--contrast '${contrastVal}'" : ""
  def contrastModeFlag = (contrastMode && contrastMode.toString().trim()) ? "--contrast_mode '${contrastMode}'" : ""

  """
  set -euo pipefail
  echo "[METAGENOMESEQ_REAL] Input: '${rds_file}'" >&2

  Rscript "${projectDir}/bin/run_metagenomeseq.R" --input "${rds_file}" --output_dir "." \
    ${condFlag} ${baseFlag} ${confounderFlag} ${contrastFlag} ${contrastModeFlag}

  # Print status warning if failed internally (process still succeeds)
  if [[ -f metagenomeseq_status.tsv ]]; then
    awk -F'\\t' 'NR==2 && \$1==0 {print "[METAGENOMESEQ_REAL] WARNING: "\$3" (stage="\$2")" > "/dev/stderr"}' metagenomeseq_status.tsv
  fi

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      metagenomeSeq: \$(Rscript -e "cat(as.character(packageVersion('metagenomeSeq')))")
      R: \$(Rscript -e "cat(R.version.string)")
  END_VERSIONS
  """
}
