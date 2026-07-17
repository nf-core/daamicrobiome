#!/usr/bin/env Rscript
source('normalize_confounder.R', chdir = TRUE)

library(corncob); packageVersion("corncob") # 0.4.2
library(phyloseq)

run_corncob_analysis <- function(ps_object, condition_var, base_condition, confounder, output_dir){
  
  conf_vec <- normalize_confounders(confounder)
  has_confounder <- length(conf_vec) > 0

  # Check condition variable exists
  if (!(condition_var %in% colnames(ps_object@sam_data))) {
    stop(sprintf("Condition variable '%s' not found in sample data", condition_var))
  }
  # Check confounder variable(s) exist
  if (has_confounder) {
    missing <- setdiff(conf_vec, colnames(ps_object@sam_data))
  if (length(missing) > 0) {
    stop(sprintf(
      "Confounder variable(s) not found in sample data: %s",
      paste(missing, collapse = ", ")
    ))
  }
}

  # Ensure condition is a factor with the correct reference level
  ps_object@sam_data[[condition_var]] <- relevel(factor(ps_object@sam_data[[condition_var]]), ref = base_condition)

  if (has_confounder) {
    rhs <- c(condition_var, conf_vec)
    full_formula <- as.formula(paste0("~", paste0(rhs, collapse = " + ")))
    null_formula <- as.formula(paste0("~", paste0(conf_vec, collapse = " + ")))
  } else {
    full_formula <- as.formula(paste0("~", condition_var))
    null_formula <- ~ 1
  }

  corncob_results <- differentialTest(
    formula = full_formula,
    formula_null = null_formula,
    phi.formula = ~ 1,
    phi.formula_null = ~ 1,
    test = "Wald",
    boot = FALSE,
    data = ps_object,
    fdr_cutoff = 1)

  # Save standardized results for heatmap including p-value
  if(length(corncob_results$significant_taxa) == 0) {
    heatmap_ready <- data.frame(
      taxon_id = character(0),
      lfc = numeric(0),
      pvalue = numeric(0),
      adj_pval = numeric(0),
      tool = character(0)
      )
    } else {

    data_for_plot <- plot(corncob_results, data_only = T)
    value <- data_for_plot$x
    tool <- rep("corncob", length(corncob_results$significant_taxa))  # repeat "corncob" for each row

    heatmap_ready <- data.frame(
      taxon_id = corncob_results$significant_taxa,
      lfc = data_for_plot$x,
      pvalue = corncob_results$p[corncob_results$significant_taxa],
      adj_pval = corncob_results$p_fdr[corncob_results$significant_taxa],
      tool = "corncob"
    )
    }

  write.table(heatmap_ready,
              file = file.path(output_dir, "corncob_da.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  return(corncob_results)
}

main <- function() {
  suppressPackageStartupMessages(library(optparse))

  option_list <- list(
    make_option(c("-i","--input"),      type = "character", help = "Input RDS path"),
    make_option(c("-c","--condition"),  type = "character", help = "Condition variable name"),
    make_option(c("-b","--base"),       type = "character", help = "Reference level"),
    make_option(c("-n","--confounder"), type = "character", default = NULL, help = "Confounder variable name"),
    make_option(c("-o","--output_dir"), type = "character", default = ".", help = "Output directory")
  )

  args <- parse_args(OptionParser(option_list = option_list))

  ps_object <- readRDS(args$input)
  run_corncob_analysis(ps_object, args$condition, args$base, args$confounder, args$output_dir)
}

# Run only when executed as a script (not when sourced)
if (sys.nframe() == 0) main()
