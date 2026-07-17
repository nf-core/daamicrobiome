include { SIMULATE_MIDASIM } from '../../../modules/local/simulation/main'

workflow SIMULATION {

  take:
  template_rds_ch    // a channel emitting one template .RDS

  main:
  sim = SIMULATE_MIDASIM(template_rds_ch)

  // build tuples: (rep_id, rds_path)
  rep_rds_ch = sim.rds_files
    .flatten()                         // turn List<Path> into stream of Path
    .map { f -> tuple(f.baseName, f) } // baseName = filename without extension

  emit:
  rds       = rep_rds_ch
  truth_dir = sim.truth_dir
  meta_dir  = sim.meta_dir
  index_csv = sim.index_csv
}
