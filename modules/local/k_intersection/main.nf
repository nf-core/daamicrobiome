process K_INTERSECTION {
  tag { "K_INTERSECTION:${rep_id}" }

  conda "${moduleDir}/environment.yml"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'oras://community.wave.seqera.io/library/r-fs_r-ggrepel_r-openxlsx_r-optparse_r-tidyverse:0b1743bb97ecab27' :
      'community.wave.seqera.io/library/r-fs_r-ggrepel_r-openxlsx_r-optparse_r-tidyverse:7c83a0ba9460a597' }"

  input:
  tuple val(rep_id), path(da_tsv_files)

  output:
  path("k_intersection_out"), emit: intersection_dir
  path("versions.yml"), emit: versions

  script:
  """
  set -euo pipefail

  mkdir -p k_intersection_out

  # Build filelist from staged DA result files
  ls *.tsv > filelist.txt

  Rscript ${projectDir}/bin/run_k_intersection.R \
    --input  "filelist.txt" \
    --outdir "k_intersection_out" \
    --alpha  "${params.consensus_alpha ?: 0.05}" \
    --lfc_min "${params.consensus_lfc_min ?: 0}"

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      openxlsx: \$(Rscript -e "cat(as.character(packageVersion('openxlsx')))")
      R: \$(Rscript -e "cat(R.version.string)")
  END_VERSIONS
  """
}
