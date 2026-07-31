#!/usr/bin/env Rscript
source('normalize_confounder.R', chdir = TRUE)

library(MicrobiomeStat); packageVersion("MicrobiomeStat")
library(phyloseq)

run_linda_analysis <- function(ps_object, condition_var, base_condition, confounder = NULL, output_dir) {
  
  confounder <- normalize_confounders(confounder)
  if (length(confounder) == 0) confounder <- NULL

  # Metadata
  metadata_linda <- data.frame(sample_data(ps_object))
  # Ensure condition is a factor with the correct reference level
  metadata_linda[[condition_var]] <- relevel(factor(metadata_linda[[condition_var]]), ref = base_condition)
  
  # build formula
  rhs <- condition_var
  if (!is.null(confounder) && length(confounder) > 0)
    rhs <- paste(rhs, paste(confounder, collapse = "+"), sep = "+")
  fml <- paste("~", rhs)

  # Run LinDA
  results_linda <- linda(
    feature.dat = otu_table(ps_object),
    meta.dat = metadata_linda,
    formula = fml
  )
  
  # Get contrasts (all non-reference levels)
  linda_contrast_cols <- setdiff(levels(metadata_linda[[condition_var]]), base_condition)
  
  # Ensure output directory exists
  dir.create(output_dir, showWarnings = FALSE)
  
  # Loop through contrasts and save separate files for each
  for(contrast in linda_contrast_cols) {
    
    linda_contrast_name <- paste0(condition_var, contrast)
    
    if(!(linda_contrast_name %in% names(results_linda$output))) {
      warning("Contrast ", linda_contrast_name, " not found. Skipping.")
      next
    }
    
    # Save standardized results for heatmap including p-value
    if(nrow(results_linda$output[[linda_contrast_name]]) == 0) {
      heatmap_ready <- data.frame(
        taxon_id = character(0),
        lfc = numeric(0),
        pvalue = numeric(0),
        adj_pval = numeric(0),
        tool = character(0)
      )
    } else {
      heatmap_ready <- data.frame(
        taxon_id = rownames(results_linda$output[[linda_contrast_name]]),
        lfc = results_linda$output[[linda_contrast_name]]$log2FoldChange,
        pvalue = results_linda$output[[linda_contrast_name]]$pvalue,
        adj_pval = results_linda$output[[linda_contrast_name]]$padj,
        tool = "linda",
        contrast = contrast # explicitly indicate contrast
    )
    }

    heatmap_ready <- heatmap_ready[, -which(colnames(heatmap_ready) == "contrast")]

    # Save results
    write.table(
      heatmap_ready, 
      file = file.path(output_dir, "linda_da.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
  }

  return(results_linda)
}

main <- function() {
  suppressPackageStartupMessages(library(optparse))

  option_list <- list(
    make_option(c("-i","--input"),      type = "character", help = "Input RDS path"),
    make_option(c("-c","--condition"),  type = "character", help = "Condition variable name"),
    make_option(c("-b","--base"),       type = "character", help = "Reference level"),
    make_option(c("-n","--confounder"), type = "character", default = NULL, help = "Confounder variable name [default %default]"),
    make_option(c("-o","--output_dir"), type = "character", default = ".", help = "Output directory [default %default]")
  )

  args <- parse_args(OptionParser(option_list = option_list))

  ps_object <- readRDS(args$input)
  conf_cols <- if (is.null(args$confounder) || args$confounder == "") NULL
              else trimws(strsplit(args$confounder, ",")[[1]])
  run_linda_analysis(ps_object, args$condition, args$base, confounder = conf_cols, output_dir = args$output_dir)
}

# Run only when executed as a script (not when sourced)
if (sys.nframe() == 0) main()
