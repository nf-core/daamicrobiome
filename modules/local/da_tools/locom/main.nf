process LOCOM {
  tag "LOCOM:${rep_id}"

  input:
  tuple val(rep_id), path(rds_file)

  output:
    // main output used by DIFF_ABUNDANCE
    tuple val(rep_id), path("locom_da.tsv"), emit: da

    // sidecar status file (optional)
    tuple val(rep_id), path("locom_status.tsv"), emit: status, optional: true
    path("versions.yml"), emit: versions

  script:
    def condFlag       = params.metadata_condition  ? "--condition '${params.metadata_condition}'" : ""
    def baseFlag       = params.metadata_base       ? "--base '${params.metadata_base}'" : ""
    def confounderFlag = params.metadata_confounder ? "--confounder '${params.metadata_confounder}'" : ""

    """
    set -euo pipefail
    echo "[LOCOM] Input: '${rds_file}'" >&2

    run_locom.R --input "${rds_file}" --output_dir "." ${condFlag} ${baseFlag} ${confounderFlag}

    # If status exists and indicates failure, print warning
    if [[ -f locom_status.tsv ]]; then
      ok=\$(tail -n +2 locom_status.tsv | head -n 1 | cut -f1)
      msg=\$(tail -n +2 locom_status.tsv | head -n 1 | cut -f2-)
      if [[ "\$ok" == "0" ]]; then
        echo "[LOCOM] WARNING: \$msg" >&2
      fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        LOCOM: \$(Rscript -e "cat(as.character(packageVersion('LOCOM')))")
        R: \$(Rscript -e "cat(R.version.string)")
    END_VERSIONS
    """
}
