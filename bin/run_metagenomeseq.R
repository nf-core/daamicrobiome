#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(metagenomeSeq)
  library(phyloseq)
  library(optparse)
  library(Biobase)
})

write_empty_metagenomeseq_tsv <- function(out_tsv) {
  empty <- data.frame(
    taxon_id  = character(0),
    lfc       = numeric(0),
    pvalue    = numeric(0),
    adj_pval  = numeric(0),
    tool      = character(0),
    stringsAsFactors = FALSE
  )
  write.table(empty, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
}

write_status <- function(status_path, ok, stage, err_msg,
                         n_samples = NA_integer_, n_taxa = NA_integer_,
                         condition_var = NA_character_,
                         base_condition = NA_character_,
                         contrast_used = NA_character_,
                         confounders = NA_character_,
                         n_levels = NA_integer_,
                         levels = NA_character_,
                         n_unique_removed = NA_integer_) {
  df <- data.frame(
    ok = as.integer(ok),
    stage = as.character(stage),
    error_message = as.character(err_msg),
    n_samples = n_samples,
    n_taxa = n_taxa,
    condition_var = as.character(condition_var),
    base_condition = as.character(base_condition),
    contrast_used = as.character(contrast_used),
    confounders = as.character(confounders),
    n_levels = n_levels,
    levels = as.character(levels),
    n_unique_removed = n_unique_removed,
    stringsAsFactors = FALSE
  )
  write.table(df, file = status_path, sep = "\t", quote = FALSE, row.names = FALSE)
}

# robust extraction of indices to remove from uniqueFeatures output
extract_unique_feature_indices <- function(uf) {
  if (is.null(uf) || length(uf) == 0) return(integer(0))

  rm_idx <- integer(0)
  if (is.matrix(uf) || is.data.frame(uf)) {
    if ("featureIndices" %in% colnames(uf)) {
      rm_idx <- uf[, "featureIndices", drop = TRUE]
    } else {
      # fall back to first column if structure differs
      rm_idx <- uf[, 1, drop = TRUE]
    }
  } else {
    # uf sometimes is a vector in edge cases
    rm_idx <- uf
  }

  rm_idx <- unique(as.integer(rm_idx))
  rm_idx <- rm_idx[!is.na(rm_idx)]
  rm_idx
}

# find the effect size column corresponding to condition_var + contrast level
detect_lfc_column <- function(results_df, condition_var, contrast_level) {
  cn <- colnames(results_df)

  # model.matrix uses make.names(level) inside column names
  lvl_s <- make.names(contrast_level)
  exact1 <- paste0(condition_var, lvl_s)
  if (exact1 %in% cn) return(exact1)

  # some outputs may have separator-like patterns; try fuzzy matching
  hits <- grep(paste0("^", condition_var), cn, value = TRUE)
  if (length(hits) == 0) return(NA_character_)

  # prefer those that contain the sanitized level
  hits2 <- hits[grepl(lvl_s, hits, fixed = TRUE)]
  if (length(hits2) == 1) return(hits2)

  # if still ambiguous, return NA and handle upstream
  NA_character_
}

# pull pvalue/adjp columns robustly
detect_p_columns <- function(results_df) {
  cn <- colnames(results_df)
  pick <- function(cands) {
    hit <- cands[cands %in% cn]
    if (length(hit) > 0) return(hit[[1]])
    # case-insensitive fallback
    hit2 <- cn[tolower(cn) %in% tolower(cands)]
    if (length(hit2) > 0) return(hit2[[1]])
    NA_character_
  }
  p_col    <- pick(c("pvalues", "pvalue", "p.value", "p_val", "p"))
  adjp_col <- pick(c("adjPvalues", "adj_pvalues", "padj", "adj.p.value", "qvalues", "qvalue"))
  list(p_col = p_col, adjp_col = adjp_col)
}

run_metagenomeseq_analysis <- function(ps_object, condition_var, base_condition,
                                       confounders = NULL,
                                       output_dir = ".",
                                       contrast_condition = NULL,
                                       contrast_mode = c("error", "most_common", "first")) {

  contrast_mode <- match.arg(contrast_mode)

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  out_tsv     <- file.path(output_dir, "metagenomeseq_da.tsv")
  status_path <- file.path(output_dir, "metagenomeseq_status.tsv")

  # defaults for status
  stage <- "init"
  ok <- TRUE
  err_msg <- ""
  contrast_used <- NA_character_
  n_unique_removed <- NA_integer_

  # always ensure a valid output exists
  write_empty_metagenomeseq_tsv(out_tsv)

  # normalize confounders vector
  if (is.null(confounders)) confounders <- character(0)
  confounders <- confounders[nzchar(confounders)]

  # wrap EVERYTHING so we never crash the process
  res <- tryCatch({

    stage <- "convert_phyloseq"
    MRexp <- phyloseq_to_metagenomeSeq(ps_object)

    stage <- "check_metadata"
    pdata_df <- Biobase::pData(MRexp)
    if (!(condition_var %in% colnames(pdata_df))) {
      stop(sprintf("Condition variable '%s' not found in sample data", condition_var))
    }
    if (length(confounders) > 0) {
      missing <- setdiff(confounders, colnames(pdata_df))
      if (length(missing) > 0) {
        stop("Confounder variable(s) not found: ", paste(missing, collapse = ", "))
      }
    }

    # factorize condition early
    pdata_df[[condition_var]] <- factor(pdata_df[[condition_var]])

    # handle levels / base existence
    if (!(base_condition %in% levels(pdata_df[[condition_var]]))) {
      stop("Base condition not present in data: ", base_condition)
    }
    if (nlevels(pdata_df[[condition_var]]) < 2) {
      stop("Condition has <2 levels; cannot test.")
    }

    # Choose contrast robustly
    levs <- levels(pdata_df[[condition_var]])
    non_base <- setdiff(levs, base_condition)

    # If contrast_condition is invalid but looks like a mode keyword, ignore it
    if (!is.null(contrast_condition) && nzchar(contrast_condition) &&
        !(contrast_condition %in% levs)) {
      if (contrast_condition %in% c("most_common", "first", "error")) {
        contrast_condition <- NULL
      } else {
        stop("Specified contrast level not present: ", contrast_condition)
      }
    }

    if (is.null(contrast_condition) || !nzchar(contrast_condition)) {
      if (length(non_base) == 1) {
        contrast_used <- non_base[[1]]
      } else {
        if (contrast_mode == "error") {
          stop("Multiple non-base levels detected (", paste(non_base, collapse = ", "),
               "). Provide --contrast or use --contrast_mode most_common/first.")
        } else if (contrast_mode == "first") {
          contrast_used <- non_base[[1]]
        } else if (contrast_mode == "most_common") {
          tab <- table(pdata_df[[condition_var]])
          tab <- tab[names(tab) != base_condition]
          contrast_used <- names(which.max(tab))
        }
      }
    } else {
      contrast_used <- contrast_condition
    }

    # uniqueFeatures filter (robust to uf type)
    stage <- "uniqueFeatures"
    uf <- uniqueFeatures(MRexp, cl = pdata_df[[condition_var]])
    rm_idx <- extract_unique_feature_indices(uf)
    n_unique_removed <- length(rm_idx)

    MR_filtered <- if (length(rm_idx) > 0) MRexp[-rm_idx, ] else MRexp

    # normalization can fail with tiny feature sets; guard it
    stage <- "cumNorm"
    p <- cumNormStat(MR_filtered)
    MR_filtered <- cumNorm(MR_filtered, p = p)

    # update pdata after filtering
    pdata_df2 <- Biobase::pData(MR_filtered)
    pdata_df2[[condition_var]] <- factor(pdata_df2[[condition_var]])
    pdata_df2[[condition_var]] <- relevel(pdata_df2[[condition_var]], ref = base_condition)
    Biobase::pData(MR_filtered) <- pdata_df2

    # quick sanity: sample/taxa counts
    stage <- "dimensions"
    n_samples <- ncol(MR_filtered)
    n_taxa <- nrow(MR_filtered)

    if (is.na(n_samples) || is.na(n_taxa) || n_samples < 4 || n_taxa < 2) {
      stop(sprintf("Too few samples/taxa after filtering (n_samples=%s, n_taxa=%s).", n_samples, n_taxa))
    }

    # Build formula safely
    stage <- "design_matrix"
    if (length(confounders) == 0) {
      form <- reformulate(condition_var)
    } else {
      form <- reformulate(c(condition_var, confounders))
    }

    mod <- model.matrix(form, data = pdata_df2)
    if (ncol(mod) < 2) {
      stop("Design matrix has <2 columns (nothing to test).")
    }

    # fitZig
    stage <- "fitZig"
    fit <- fitZig(MR_filtered, mod = mod)

    # MRtable: prefer full table (do not apply eff threshold; scoring expects full)
    stage <- "MRtable"
    results <- MRtable(fit, number = Inf)

    if (is.null(results) || nrow(results) == 0) {
      # empty but not a crash
      write_empty_metagenomeseq_tsv(out_tsv)
      write_status(
        status_path, ok = TRUE, stage = "MRtable_empty", err_msg = "",
        n_samples = n_samples, n_taxa = n_taxa,
        condition_var = condition_var, base_condition = base_condition,
        contrast_used = contrast_used,
        confounders = paste(confounders, collapse = ","),
        n_levels = nlevels(pdata_df2[[condition_var]]),
        levels = paste(levels(pdata_df2[[condition_var]]), collapse = ","),
        n_unique_removed = n_unique_removed
      )
      return(invisible(results))
    }

    # detect columns
    stage <- "detect_columns"
    lfc_col <- detect_lfc_column(results, condition_var, contrast_used)
    pcs <- detect_p_columns(results)

    if (is.na(lfc_col)) {
      stop("Could not identify effect-size column for contrast '", contrast_used,
           "'. Available columns: ", paste(colnames(results), collapse = ", "))
    }
    if (is.na(pcs$p_col) || is.na(pcs$adjp_col)) {
      stop("Could not identify pvalue/adjp columns. Available columns: ",
           paste(colnames(results), collapse = ", "))
    }

    # write standardized output
    stage <- "write_output"
    heatmap_ready <- data.frame(
      taxon_id = rownames(results),
      lfc      = suppressWarnings(as.numeric(results[[lfc_col]])),
      pvalue   = suppressWarnings(as.numeric(results[[pcs$p_col]])),
      adj_pval = suppressWarnings(as.numeric(results[[pcs$adjp_col]])),
      tool     = "metagenomeseq",
      stringsAsFactors = FALSE
    )

    # Ensure no NA columns due to coercion issues
    # (keep as NA numeric; scoring should handle)
    write.table(heatmap_ready, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

    # status OK
    write_status(
      status_path, ok = TRUE, stage = "ok", err_msg = "",
      n_samples = n_samples, n_taxa = n_taxa,
      condition_var = condition_var, base_condition = base_condition,
      contrast_used = contrast_used,
      confounders = paste(confounders, collapse = ","),
      n_levels = nlevels(pdata_df2[[condition_var]]),
      levels = paste(levels(pdata_df2[[condition_var]]), collapse = ","),
      n_unique_removed = n_unique_removed
    )

    invisible(results)

  }, error = function(e) {
    ok <<- FALSE
    err_msg <<- conditionMessage(e)

    # keep empty output already written
    write_status(
      status_path, ok = FALSE, stage = stage, err_msg = err_msg,
      condition_var = condition_var, base_condition = base_condition,
      contrast_used = contrast_used,
      confounders = paste(confounders, collapse = ","),
      n_unique_removed = n_unique_removed
    )
    invisible(NULL)
  })

  return(res)
}

main <- function() {
  option_list <- list(
    make_option(c("-i","--input"),      type = "character", help = "Input RDS path"),
    make_option(c("-c","--condition"),  type = "character", help = "Condition variable name"),
    make_option(c("-b","--base"),       type = "character", help = "Reference level"),
    make_option(c("-n","--confounder"), type = "character", default = NULL,
                help = "Confounders: single name or comma-separated list"),
    make_option(c("-o","--output_dir"), type = "character", default = ".", help = "Output directory"),
    make_option(c("-t","--contrast"),   type = "character", default = NULL,
                help = "Optional: specific condition level to contrast against base"),
    make_option(c("--contrast_mode"),   type = "character", default = "most_common",
                help = "If contrast not provided and >2 non-base levels exist: error|most_common|first")
  )

  args <- parse_args(OptionParser(option_list = option_list))

  ps_object <- readRDS(args$input)

  contrast_raw <- args$contrast
  if (!is.null(contrast_raw)) {
    contrast_raw <- trimws(contrast_raw)
    contrast_raw <- gsub("^\\s*\\[|\\]\\s*$", "", contrast_raw)   # strip surrounding [ ]
    contrast_raw <- gsub("[\"']", "", contrast_raw)              # strip quotes
    if (!nzchar(contrast_raw)) contrast_raw <- NULL
  }

  contrast_mode <- args$contrast_mode
  if (is.null(contrast_mode) || !nzchar(contrast_mode)) {
    contrast_mode <- "most_common"
  }

  # If contrast looks like a mode keyword (or equals the mode), ignore it
  if (!is.null(contrast_raw) &&
      (contrast_raw %in% c("most_common", "first", "error") || identical(contrast_raw, contrast_mode))) {
    contrast_raw <- NULL
  }

  args$contrast <- contrast_raw
  args$contrast_mode <- contrast_mode

  conf_cols <- character(0)
  if (!is.null(args$confounder) && nzchar(args$confounder)) {
    raw <- args$confounder
    raw <- gsub("^\\s*\\[|\\]\\s*$", "", raw)   # strip surrounding [ ]
    raw <- gsub("[\"']", "", raw)              # strip quotes
    conf_cols <- trimws(unlist(strsplit(raw, ",")))
    conf_cols <- conf_cols[nzchar(conf_cols)]
  }

  run_metagenomeseq_analysis(
    ps_object = ps_object,
    condition_var = args$condition,
    base_condition = args$base,
    confounders = conf_cols,
    output_dir = args$output_dir,
    contrast_condition = args$contrast,
    contrast_mode = args$contrast_mode
  )

  # Exit 0 even on internal errors — status is reported via metagenomeseq_status.tsv
  quit(status = 0)
}

if (sys.nframe() == 0) main()
