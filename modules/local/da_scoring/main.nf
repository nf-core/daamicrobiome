process DA_SCORING {
  tag { "DA_SCORING" }

  conda "${moduleDir}/environment.yml"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'oras://community.wave.seqera.io/library/r-fs_r-ggrepel_r-openxlsx_r-optparse_r-tidyverse:0b1743bb97ecab27' :
      'community.wave.seqera.io/library/r-fs_r-ggrepel_r-openxlsx_r-optparse_r-tidyverse:7c83a0ba9460a597' }"

  input:
  val(done)
  path(sim_root_dir)
  path(truth_root_dir)
  val(tools_csv)

  output:
  path("scoring_output"), emit: scoring_dir
  path("versions.yml"), emit: versions

  script:
  """
  set -euo pipefail

  echo "[DA_SCORING] sim_root_dir  : ${sim_root_dir}" >&2
  echo "[DA_SCORING] truth_root_dir: ${truth_root_dir}" >&2
  echo "[DA_SCORING] tools         : ${tools_csv}" >&2

  Rscript ${projectDir}/bin/run_scoring.R \
    --base_results   "${sim_root_dir}" \
    --base_truth     "${truth_root_dir}" \
    --outdir         "scoring_output" \
    --tools          "${tools_csv}" \
    --alpha          "${params.eval_alpha ?: 0.05}" \
    --fdr_target     "${params.eval_fdr_target ?: 0.05}" \
    --fdr_reference  "${params.eval_fdr_reference ?: 0.10}" \
    --retention_at_ref "${params.eval_retention_at_ref ?: 0.25}" \
    --prefix         "${params.sim_prefix ?: ''}"

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      R: \$(Rscript -e "cat(R.version.string)")
  END_VERSIONS
  """
}
