process LOCOM_REAL {
  tag "LOCOM_REAL:${rep_id}"

  conda "${moduleDir}/environment.yml"

  input:
  tuple val(rep_id), path(rds_file)

  output:
    // main output used by DIFF_ABUNDANCE
    tuple val(rep_id), path("locom_da.tsv"), emit: da

    // sidecar status file (optional)
    tuple val(rep_id), path("locom_status.tsv"), emit: status, optional: true
    path("versions.yml"), emit: versions

  script:
    def condFlag       = params.condition  ? "--condition '${params.condition}'" : ""
    def baseFlag       = params.base_level ? "--base '${params.base_level}'"    : ""
    def conf = params.confounder
    def conf_str = (conf instanceof List) ? conf.join(',') : (conf ? conf.toString() : "")
    conf_str = conf_str.trim()
    def confounderFlag = conf_str ? "--confounder '${conf_str}'" : ""

    """
    set -euo pipefail
    echo "[LOCOM_REAL] Input: '${rds_file}'" >&2

    run_locom.R --input "${rds_file}" --output_dir "." ${condFlag} ${baseFlag} ${confounderFlag}

    # If status exists and indicates failure, print warning
    if [[ -f locom_status.tsv ]]; then
      ok=\$(tail -n +2 locom_status.tsv | head -n 1 | cut -f1)
      msg=\$(tail -n +2 locom_status.tsv | head -n 1 | cut -f2-)
      if [[ "\$ok" == "0" ]]; then
        echo "[LOCOM_REAL] WARNING: \$msg" >&2
      fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        LOCOM2: \$(Rscript -e "cat(as.character(packageVersion('LOCOM2')))")
        R: \$(Rscript -e "cat(R.version.string)")
    END_VERSIONS
    """
}