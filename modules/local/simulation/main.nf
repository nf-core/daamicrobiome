process SIMULATE_MIDASIM {
  tag { "SIMULATE_MIDASIM:${template_rds.getBaseName()}" }

  input:
  path(template_rds)

  output:
  path("rds/*.RDS"),   emit: rds_files
  path("truth"),       emit: truth_dir
  path("meta"),        emit: meta_dir
  path("index.csv"),   emit: index_csv
  path("versions.yml"), emit: versions

  script:
  def n_ctrl   = params.sim_n_ctrl  != null ? "--n_ctrl ${params.sim_n_ctrl}" : ""
  def n_case   = params.sim_n_case  != null ? "--n_case ${params.sim_n_case}" : ""
  def p_da     = (params.sim_p_da instanceof List) ? params.sim_p_da.join(',') : "${params.sim_p_da}"
  def abs_lfc  = (params.sim_abs_logfc instanceof List) ? params.sim_abs_logfc.join(',') : "${params.sim_abs_logfc}"
  def design   = (params.sim_design instanceof List) ? params.sim_design.join(',') : "${params.sim_design}"
  def subset_v = params.sim_subset_var  ? "--subset_var '${params.sim_subset_var}'" : ""
  def subset_k = params.sim_subset_keep ? "--subset_keep '${(params.sim_subset_keep instanceof List) ? params.sim_subset_keep.join(',') : params.sim_subset_keep}'" : ""
  def batch_c  = params.sim_batch_col   ? "--batch_col '${params.sim_batch_col}'" : ""
  def batch_sd = params.sim_batch_sd    != null ? "--batch_sd ${params.sim_batch_sd}" : "--batch_sd 0.3"
  def seed0    = params.sim_seed0       != null ? "--seed0 ${params.sim_seed0}" : "--seed0 123"
  def prefix   = params.sim_prefix      ? "--prefix '${params.sim_prefix}'" : "--prefix 'ps_midasim'"

  """
  set -euo pipefail
  mkdir -p rds truth meta

  Rscript ${projectDir}/bin/run_simulation.R \
    --input "${template_rds}" \
    --outdir "." \
    ${n_ctrl} ${n_case} \
    --p_da "${p_da}" \
    --abs_logfc "${abs_lfc}" \
    --n_rep "${params.sim_n_rep}" \
    --design "${design}" \
    ${subset_v} ${subset_k} \
    ${batch_c} ${batch_sd} \
    ${prefix} \
    ${seed0}

  # Ensure outputs are exactly where Nextflow expects them
  test -f index.csv
  test -d rds
  test -d truth
  test -d meta

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      MIDASim: \$(Rscript -e "cat(as.character(packageVersion('MIDASim')))")
      R: \$(Rscript -e "cat(R.version.string)")
  END_VERSIONS
  """
}
