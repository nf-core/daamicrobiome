#!/usr/bin/env Rscript

# Inputs:
#   - ps_object:       A `phyloseq` object containing the microbiome dataset.
#   - condition_var:   A string indicating the name of the metadata variable to test (e.g. "Group").
#   - base_condition:  A string indicating the reference level within the condition variable (e.g. "Healthy").
#   - output_dir:      Output directory
#   - otu_filter:      value btw 0 and 1 to filter out OTUs present in less than value * tot samples
# Output:
#   - LOCOM results
#   - Also saves a TSV file to "locom_results.tsv"
source('normalize_confounder.R', chdir = TRUE)
library(LOCOM2); packageVersion("LOCOM2")
library(phyloseq)

run_locom_analysis <- function(ps_object, condition_var, base_condition,
                               confounder = NULL, output_dir,
                               otu_filter = 0, fdr.nominal = 0.5) {

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Always write a status file so users can see if LOCOM failed
  status_path <- file.path(output_dir, "locom_status.tsv")
  write_status <- function(ok, err_msg = "", n_samples = NA_integer_, n_taxa = NA_integer_) {
    status_df <- data.frame(
      ok = as.integer(ok),
      error_message = err_msg,
      n_samples = n_samples,
      n_taxa = n_taxa,
      condition_var = condition_var,
      base_condition = base_condition,
      confounders = if (is.null(confounder) || length(confounder) == 0) "" else paste(confounder, collapse = ","),
      otu_filter = otu_filter,
      fdr_nominal = fdr.nominal,
      stringsAsFactors = FALSE
    )
    write.table(status_df, status_path, sep = "\t", quote = FALSE, row.names = FALSE)
  }

  # Always write a valid locom_da.tsv (even if empty)
  out_tsv <- file.path(output_dir, "locom_da.tsv")
  write_empty_locom_tsv <- function() {
    empty <- data.frame(
      taxon_id = character(0),
      lfc      = numeric(0),
      pvalue   = numeric(0),
      adj_pval = numeric(0),
      tool     = character(0),
      stringsAsFactors = FALSE
    )
    write.table(empty, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  }

  # -------------------------
  # Extract count table + metadata
  # -------------------------
  otu_mat <- t(as(otu_table(ps_object), "matrix"))

  # drop empty samples (this is safe/necessary; many tools effectively do this)
  otu_mat <- otu_mat[rowSums(otu_mat) > 0, , drop = FALSE]

  sample_df <- as.data.frame(sample_data(ps_object))

  if (!(condition_var %in% colnames(sample_df))) {
    stop(sprintf("Condition variable '%s' not found in sample data", condition_var))
  }

  # IMPORTANT: reorder sample_df to match otu_mat row order
  common_ids <- intersect(rownames(sample_df), rownames(otu_mat))
  otu_mat    <- otu_mat[common_ids, , drop = FALSE]
  sample_df  <- sample_df[common_ids, , drop = FALSE]

  n_samples <- nrow(otu_mat)
  n_taxa    <- ncol(otu_mat)

  if (n_samples < 4 || n_taxa < 2) {
    msg <- sprintf("Too few samples/taxa after dropping empty samples (n_samples=%d, n_taxa=%d).", n_samples, n_taxa)
    message("[LOCOM2] WARNING: ", msg)
    write_empty_locom_tsv()
    write_status(ok = FALSE, err_msg = msg, n_samples = n_samples, n_taxa = n_taxa)
    return(NULL)
  }

  # Ensure condition is a factor with the correct reference level
  sample_df[[condition_var]] <- relevel(factor(sample_df[[condition_var]]), ref = base_condition)

  # Binary outcome: 0 = base condition, 1 = all other
  Y <- ifelse(sample_df[[condition_var]] == base_condition, 0, 1)

  # -------------------------
  # Build confounder matrix C (same as your code)
  # -------------------------
  build_C <- function(df, conf_cols) {
    if (is.null(conf_cols) || length(conf_cols) == 0) return(NULL)
    df <- as.data.frame(df, stringsAsFactors = FALSE)

    if (!all(conf_cols %in% colnames(df))) {
      stop("Confounder columns not found: ",
           paste(setdiff(conf_cols, colnames(df)), collapse = ", "))
    }

    conf_df <- df[, conf_cols, drop = FALSE]
    for (nm in names(conf_df)) {
      if (inherits(conf_df[[nm]], "Rle")) conf_df[[nm]] <- as.vector(conf_df[[nm]])
      if (is.list(conf_df[[nm]]) && !is.factor(conf_df[[nm]]))
        conf_df[[nm]] <- unlist(conf_df[[nm]], use.names = FALSE)
      if (is.character(conf_df[[nm]]) || is.logical(conf_df[[nm]]))
        conf_df[[nm]] <- factor(conf_df[[nm]])
    }

    conf_df <- data.frame(conf_df, check.names = TRUE)
    C <- model.matrix(~ 0 + ., data = conf_df)

    if (ncol(C) > 0) {
      keep <- apply(C, 2, function(x) var(as.numeric(x), na.rm = TRUE) > 0)
      C <- C[, keep, drop = FALSE]
      if (ncol(C) == 0) C <- NULL
    } else C <- NULL

    C
  }

  confounder <- normalize_confounders(confounder)
  if (length(confounder) == 0) confounder <- NULL

  C <- build_C(sample_df, confounder)

  # -------------------------
  # Run LOCOM safely
  # -------------------------
  output_locom <- NULL
  ok <- TRUE
  err_msg <- ""

  output_locom <- tryCatch({
    locom2(
      otu.table    = otu_mat,
      Y            = Y,
      C            = C,
      seed         = 123,
      filter       = (otu_filter > 0),
      fdr.nominal  = fdr.nominal,
      n.perm.max   = 1000,
      n.rej.stop   = 100,
      n.cores      = 1
    )
  }, error = function(e) {
    ok <<- FALSE
    err_msg <<- conditionMessage(e)
    message("[LOCOM2] WARNING: LOCOM failed with error: ", err_msg)
    NULL
  })

  # -------------------------
  # Write output TSV (always)
  # -------------------------
  if (is.null(output_locom) || is.null(output_locom$beta) || length(output_locom$beta) == 0) {
    write_empty_locom_tsv()
  } else {
    heatmap_ready <- data.frame(
      taxon_id = names(output_locom$beta),
      lfc      = as.numeric(output_locom$beta),
      pvalue   = as.numeric(output_locom$p.otu.Wald),
      adj_pval = as.numeric(output_locom$q.otu.Wald),
      tool     = "locom",
      stringsAsFactors = FALSE
    )
    write.table(heatmap_ready, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  }

  # Write status (always)
  write_status(ok = ok, err_msg = err_msg, n_samples = n_samples, n_taxa = n_taxa)

  return(output_locom)
}


main <- function() {
  suppressPackageStartupMessages(library(optparse))

  option_list <- list(
    make_option(c("-i","--input"),      type = "character", help = "Input RDS path"),
    make_option(c("-c","--condition"),  type = "character", help = "Condition variable name"),
    make_option(c("-b","--base"),       type = "character", help = "Reference level"),
    make_option(c("-n","--confounder"), type = "character", default = NULL,
                help = "Confounder variable name (e.g. 'BioProject') or comma-separated list (e.g. 'BioProject,Sex,Age')"),
    make_option(c("-o","--output_dir"), type = "character", default = ".", help = "Output directory [default %default]"),
    make_option(c("-f","--otu_filter"),  type="double",   default=0,
                help="LOCOM filter.thresh in [0,1]. OTUs present in fewer than this *100%% of samples are filtered [default %default]")
  )

  args <- parse_args(OptionParser(option_list = option_list))

  ps_object <- readRDS(args$input)
  conf_cols <- NULL
  if (!is.null(args$confounder)) {
    conf_cols <- trimws(unlist(strsplit(args$confounder, ",")))
    conf_cols <- conf_cols[nzchar(conf_cols)]
    if (length(conf_cols) == 0) conf_cols <- NULL
  }
  
  run_locom_analysis(ps_object, args$condition, args$base, conf_cols, args$output_dir)
}

# Run only when executed as a script (not when sourced)
if (sys.nframe() == 0) main()
