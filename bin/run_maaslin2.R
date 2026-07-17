#!/usr/bin/env Rscript
source('normalize_confounder.R', chdir = TRUE)
library(Maaslin2)
library(phyloseq)

run_maaslin2_analysis <- function(ps_object, condition_var, base_condition, confounder = NULL,
                                  bc_confounder = NULL, output_dir,
                                  min_prevalence = 0, max_significance = 0.5) {

  metadata_df <- data.frame(sample_data(ps_object), check.names = FALSE, stringsAsFactors = FALSE)
  
  confounders <- normalize_confounders(confounder)
  if (length(confounders) == 0) confounders <- character(0)

  if (!(condition_var %in% colnames(metadata_df))) {
    stop(sprintf("Condition variable '%s' not found in sample data", condition_var))
  }

  if (length(confounders) > 0) {
    missing <- setdiff(confounders, colnames(metadata_df))
    if (length(missing) > 0) {
      stop("Confounder variable(s) not found: ", paste(missing, collapse = ", "))
    }
  }

  # Relevel the condition variable
  metadata_df[[condition_var]] <- relevel(factor(metadata_df[[condition_var]]), ref = base_condition)

  # Extract OTU and metadata tables
  data_file <- data.frame(t(otu_table(ps_object)), check.names = FALSE)
  # Map sanitized names back to original taxa names to avoid Maaslin2 make.names changes
  orig_taxa <- colnames(data_file)
  safe_taxa <- make.names(orig_taxa, unique = TRUE)
  colnames(data_file) <- safe_taxa
  name_map <- setNames(orig_taxa, safe_taxa)
  metadata_file <- metadata_df

  # build formula (fixed effect, reference)
  fixed_effects <- c(condition_var, confounders)
  reference <- paste(condition_var, base_condition, sep = ",")
  if (!is.null(bc_confounder) && length(confounders) > 0) {
    reference <- c(reference, paste(confounders[1], bc_confounder, sep = ","))
  }

  # Run Maaslin2 (disable plots to avoid ggplot2 incompatibilities)
  maas_args <- list(
    input_data = data_file,
    input_metadata = metadata_file,
    min_prevalence = min_prevalence,
    max_significance = max_significance,
    fixed_effects = fixed_effects,
    reference = reference,
    output = output_dir
  )
  maas_formals <- names(formals(Maaslin2::Maaslin2))
  plot_formals <- maas_formals[grepl("^plot", maas_formals)]
  if (length(plot_formals) > 0) {
    for (p in plot_formals) maas_args[[p]] <- FALSE
  }

  fit <- do.call(Maaslin2, maas_args)
  
  # Load significant results
  all_results_path <- file.path(output_dir, "all_results.tsv")
  if (!file.exists(all_results_path)) {
    stop("No all_results.tsv found. Check if Maaslin2 ran successfully.")
  }
  
  maaslin2_all <- read.table(all_results_path, header = TRUE, sep = "\t")
  # Keep only rows for the main condition (avoid duplicate taxa from confounders)
  maaslin2_all <- subset(maaslin2_all, metadata == condition_var)
  # Restore original taxa names
  if (nrow(maaslin2_all) > 0 && "feature" %in% colnames(maaslin2_all)) {
    maaslin2_all$feature <- unname(name_map[maaslin2_all$feature])
  }

  # Save standardized results
  if (nrow(maaslin2_all) == 0) {
    heatmap_ready <- data.frame(
      taxon_id = character(0),
      lfc = numeric(0),
      pvalue = numeric(0),
      adj_pval = numeric(0),
      tool = character(0)
      )
    } else {
    heatmap_ready <- data.frame(
      taxon_id = maaslin2_all$feature,
      lfc = maaslin2_all$coef,
      pvalue = maaslin2_all$pval,
      adj_pval = maaslin2_all$qval,
      tool = "maaslin2"
      )
    }
  write.table(heatmap_ready,
              file = file.path(output_dir, "maaslin2_da.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  
  return(fit)
}

main <- function() {
  suppressPackageStartupMessages(library(optparse))

  option_list <- list(
    make_option(c("-i","--input"),      type = "character", help = "Input RDS path"),
    make_option(c("-c","--condition"),  type = "character", help = "Condition variable name"),
    make_option(c("-b","--base"),       type = "character", help = "Reference level"),
    make_option(c("-n","--confounder"), type = "character", default = NULL,
                help = "Confounder variable name (e.g. 'BioProject') or comma-separated list (e.g. 'BioProject1,BioProject2')"),
    make_option(c("-f","--bc_confounder"), type = "character", default = NULL,
                help = "Base categoryname of the confounder (e.g. 'BioProject number with more samples')"),
    make_option(c("-o","--output_dir"), type = "character", default = ".", help = "Output directory [default %default]")
  )

  args <- parse_args(OptionParser(option_list = option_list))

  ps_object <- readRDS(args$input)
  conf_cols <- NULL
  if (!is.null(args$confounder) && nzchar(args$confounder)) {
    conf_cols <- trimws(unlist(strsplit(args$confounder, ",")))
    conf_cols <- conf_cols[nzchar(conf_cols)]
    if (length(conf_cols) == 0) conf_cols <- NULL
  }
  run_maaslin2_analysis(ps_object, args$condition, args$base, conf_cols, args$bc_confounder, args$output_dir)
}

# Run only when executed as a script (not when sourced)
if (sys.nframe() == 0) main()
