#!/usr/bin/env Rscript
source('normalize_confounder.R', chdir = TRUE)

library(ADAPT)
library(phyloseq)

run_adapt_analysis <- function(ps_object, condition_var, base_condition, confounder, output_dir) {

  if (!inherits(ps_object, "phyloseq")) {
    stop("The input data must be a phyloseq object.")
  }

  confounder <- normalize_confounders(confounder)
  if (length(confounder) == 0) confounder <- NULL
  # If confounders provided, ensure they exist
  if (!is.null(confounder)) {
    samp <- as(sample_data(ps_object), "data.frame")
    if (!all(confounder %in% colnames(samp))) {
      stop("Confounder columns not found: ",
           paste(setdiff(confounder, colnames(samp)), collapse = ", "))
    }
  }
  sample_data(ps_object)[[condition_var]] <- relevel(factor(sample_data(ps_object)[[condition_var]]), ref = base_condition)

  cat("Running ADAPT...\n")
  adapt_result <- adapt(ps_object, condition_var, base.cond = base_condition, 
                        adj.var = confounder,
                        prev.filter = 0, depth.filter = 200, alpha = 0.4)

  DA_taxa <- summary(adapt_result, select = "all")
  if (is.null(DA_taxa) && is.null(DA_taxa)) {
    heatmap_ready <- data.frame(
      taxon_id = character(0),
      lfc = numeric(0),
      pvalue = numeric(0),
      adj_pval = numeric(0),
      tool = character(0)
      )
    } else {
    heatmap_ready <- data.frame(
      taxon_id = rownames(DA_taxa),
      lfc = DA_taxa$log10foldchange,
      pvalue = DA_taxa$pval,
      adj_pval = DA_taxa$adjusted_pval,
      tool = "ADAPT"
      )
    }

  dir.create("ADAPT_output", showWarnings = FALSE)
  write.table(heatmap_ready,
              file = file.path(output_dir, "adapt_da.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  return(adapt_result)
}


main <- function() {
  suppressPackageStartupMessages(library(optparse))

  option_list <- list(
    make_option(c("-i","--input"),      type = "character", help = "Input RDS path"),
    make_option(c("-c","--condition"),  type = "character", help = "Condition variable name"),
    make_option(c("-b","--base"),       type = "character", help = "Reference level"),
    make_option(c("-n","--confounder"), type = "character", default = NULL,
                help = "Confounder variable name (e.g. 'BioProject') or comma-separated list (e.g. 'BioProject,Sex,Age')"),
    make_option(c("-o","--output_dir"), type = "character", default = ".", help = "Output directory [default %default]")
  )

  args <- parse_args(OptionParser(option_list = option_list))

  ps_object <- readRDS(args$input)
  run_adapt_analysis(ps_object, args$condition, args$base, args$confounder, args$output_dir)
}

# Run only when executed as a script (not when sourced)
if (sys.nframe() == 0) main()