#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(phyloseq)
  library(MIDASim)
  library(tibble)
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(jsonlite)
})

# --- helper: largest-remainder allocation (Hamilton) ---
.apportion <- function(total, weights) {
  w <- as.numeric(weights)
  q <- total * w / sum(w)
  base <- floor(q)
  rem <- total - sum(base)
  if (rem > 0) {
    ord <- order(q - base, decreasing = TRUE)
    base[ord[seq_len(rem)]] <- base[ord[seq_len(rem)]] + 1L
  }
  if (!is.null(names(weights))) base <- setNames(base, names(weights))
  base
}

# --- helper: taxa rank–abundance dataframe from otu matrix (taxa x samples) ---
.rank_abundance_df <- function(otu_mat, eps = 1e-12) {
  if (is.null(dim(otu_mat)) || nrow(otu_mat) == 0 || ncol(otu_mat) == 0) {
    return(tibble(taxon = character(), mean_rel = numeric(), rank = integer()))
  }

  lib <- colSums(otu_mat)
  keep <- lib > 0
  otu_mat <- otu_mat[, keep, drop = FALSE]
  lib <- lib[keep]

  rel <- sweep(otu_mat, 2, lib, "/")
  mean_rel <- rowMeans(rel, na.rm = TRUE)

  tibble(
    taxon = rownames(otu_mat),
    mean_rel = pmax(mean_rel, eps)
  ) %>%
    arrange(desc(mean_rel)) %>%
    mutate(rank = row_number())
}

.save_rank_abundance_plot <- function(plot, png_path, width = 8, height = 5, dpi = 150) {
  # ensure directory exists
  dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)

  # Prefer ragg if available (best for headless)
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggplot2::ggsave(
      filename = png_path, plot = plot, device = ragg::agg_png,
      width = width, height = height, units = "in", dpi = dpi
    )
    return(png_path)
  }

  # Try cairo PNG if available
  if (capabilities("cairo")) {
    ggplot2::ggsave(
      filename = png_path, plot = plot,
      device = function(...) grDevices::png(..., type = "cairo"),
      width = width, height = height, units = "in", dpi = dpi
    )
    return(png_path)
  }

  # Fallback: save PDF if PNG isn't possible
  pdf_path <- sub("\\.png$", ".pdf", png_path)
  ggplot2::ggsave(
    filename = pdf_path, plot = plot, device = "pdf",
    width = width, height = height, units = "in"
  )
  pdf_path
}

# --- helper: rank–abundance plot (mean relative abundance, log scale) ---
.plot_rank_abundance <- function(df, title = NULL) {
  ggplot(df, aes(x = rank, y = mean_rel)) +
    geom_line() +
    scale_y_log10() +
    labs(
      title = title %||% "Taxa rank–abundance (mean relative abundance)",
      x = "Taxon rank (sorted by mean relative abundance)",
      y = "Mean relative abundance (log10)"
    ) +
    theme_bw()
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# --- core simulator for ONE scenario tag ---
.simulate_one <- function(fit, base_mu_len, tax_table_df, batch_levels, batch_fac,
                          lib_by_batch, n_ctrl, n_case, p_da, abs_lfc,
                          design_tag, rep_id, batch_sd, outdir,
                          prefix = "ps_midasim", seed0 = 123) {

  # deterministic seed per scenario
  scen_hash <- abs(as.integer(1e6 * ((design_tag == "confounded") + (p_da*100) + (abs_lfc*10) + rep_id )))
  set.seed(seed0 + scen_hash)

  taxa_names <- names(fit$mean.rel.abund)
  p <- length(taxa_names)
  stopifnot(p == base_mu_len)

  # 1) sample counts per batch + split groups
  healthy_n_by_batch <- table(factor(batch_fac, levels = batch_levels))
  prop_by_batch <- as.numeric(healthy_n_by_batch) / sum(healthy_n_by_batch)

  n_total <- n_ctrl + n_case
  n_by_batch <- .apportion(n_total, prop_by_batch); names(n_by_batch) <- batch_levels

  if (design_tag == "balanced") {
    ctrl_by_batch <- .apportion(n_ctrl, n_by_batch)
    case_by_batch <- n_by_batch - ctrl_by_batch
  } else {
    main_b <- names(which.max(n_by_batch))
    case_by_batch <- pmin(n_by_batch, as.integer(round(n_by_batch * 0.3)))
    case_by_batch[main_b] <- min(n_by_batch[main_b],
                                 case_by_batch[main_b] + as.integer(0.5 * n_by_batch[main_b]))

    adjust_to <- function(x, target, caps) {
      diff <- target - sum(x)
      if (diff > 0) {
        spare <- caps - x; ord <- order(spare, decreasing = TRUE); i <- 1L
        while (diff > 0) {
          b <- ord[i]
          if (spare[b] > 0) { x[b] <- x[b] + 1L; spare[b] <- spare[b] - 1L; diff <- diff - 1L }
          i <- if (i < length(ord)) i + 1L else 1L
        }
      } else if (diff < 0) {
        ord <- order(x, decreasing = TRUE); i <- 1L
        while (diff < 0) {
          b <- ord[i]
          if (x[b] > 0) { x[b] <- x[b] - 1L; diff <- diff + 1L }
          i <- if (i < length(ord)) i + 1L else 1L
        }
      }
      x
    }
    case_by_batch <- adjust_to(case_by_batch, n_case, caps = n_by_batch)
    ctrl_by_batch <- n_by_batch - case_by_batch
  }
  stopifnot(sum(ctrl_by_batch) == n_ctrl, sum(case_by_batch) == n_case)

  # 2) materialize per-sample batch + group vectors
  batch_vec <- group_vec <- character(0)
  for (b in batch_levels) {
    nc <- as.integer(ctrl_by_batch[[b]])
    ns <- as.integer(case_by_batch[[b]])
    if (is.na(nc) || is.na(ns) || nc < 0L || ns < 0L) {
      stop("Bad counts for batch ", b, ": nc=", nc, " ns=", ns)
    }
    batch_vec <- c(batch_vec, rep(b, nc), rep(b, ns))
    group_vec <- c(group_vec, rep("control", nc), rep("case", ns))
  }
  n <- length(group_vec)

  # 3) choose DA taxa and signed logFC
  n_da <- max(1L, floor(p_da * p)); da_idx <- sample.int(p, n_da)
  logFC <- numeric(p); logFC[da_idx] <- abs_lfc * sample(c(-1, 1), n_da, TRUE)
  names(logFC) <- taxa_names

  # 4) per-sample mean compositions (batch offsets + DA)
  base_mu <- fit$mean.rel.abund
  M <- matrix(rep(base_mu, each = n), nrow = n, ncol = p)

  batch_offsets <- sapply(batch_levels, function(.) rnorm(p, 0, batch_sd))
  colnames(batch_offsets) <- batch_levels; rownames(batch_offsets) <- taxa_names

  for (i in seq_len(n)) {
    b <- batch_vec[i]
    M[i, ] <- M[i, ] * exp(batch_offsets[, b])
  }
  case_rows <- which(group_vec == "case")
  if (length(case_rows)) {
    M[case_rows, da_idx] <- sweep(M[case_rows, da_idx], 2, exp(logFC[da_idx]), `*`)
  }

  # renormalize to simplex
  rs <- rowSums(M); stopifnot(all(is.finite(M)), all(rs > 0)); M <- M / rs

  # 5) library sizes per sample (resample within batch)
  lib_sizes <- vapply(seq_len(n), function(i) {
    b <- batch_vec[i]
    sample(lib_by_batch[[b]], 1)
  }, numeric(1))

  # 6) MIDASim simulate counts
  fit2 <- MIDASim.modify(fitted = fit, lib.size = lib_sizes, individual.rel.abund = M)
  sim  <- MIDASim(fitted.modified = fit2, only.rel = FALSE)
  Y <- sim$sim_count  # rows = samples, cols = taxa

  # sample IDs (more robust than factor->tabulate ambiguity)
  b_int <- as.integer(factor(batch_vec, levels = batch_levels))
  rn_idx <- ave(b_int, b_int, FUN = seq_along)
  rownames(Y) <- sprintf("%s_%03d", batch_vec, rn_idx)

  # 7) build phyloseq (taxa AS ROWS) + persist
  samp <- data.frame(
    group      = factor(group_vec, levels = c("control", "case")),
    Batch      = factor(batch_vec,  levels = batch_levels),
    p_DA       = p_da,
    abs_logFC  = abs_lfc,
    design     = design_tag,
    replicate  = rep_id,
    row.names  = rownames(Y)
  )

  tax_df_use <- NULL
  if (!is.null(tax_table_df)) {
    tax_df_use <- as.data.frame(tax_table_df)
    tax_df_use <- tax_df_use[intersect(rownames(tax_df_use), colnames(Y)), , drop = FALSE]
    missing <- setdiff(colnames(Y), rownames(tax_df_use))
    if (length(missing)) {
      add <- matrix(NA_character_, nrow = length(missing), ncol = ncol(tax_df_use),
                    dimnames = list(missing, colnames(tax_df_use)))
      tax_df_use <- rbind(tax_df_use, add)
    }
    tax_df_use <- tax_df_use[colnames(Y), , drop = FALSE]
  }

  OTU <- otu_table(t(Y), taxa_are_rows = TRUE)
  SAM <- sample_data(samp)
  PS  <- if (is.null(tax_df_use)) phyloseq(OTU, SAM) else phyloseq(OTU, SAM, tax_table(as.matrix(tax_df_use)))

  # ---- NAMING (important for downstream scoring) ----
  tag_base <- sprintf("pda%02d_logfc%.1f_%s_rep%02d",
                      as.integer(100*p_da), abs_lfc, design_tag, rep_id)
  tag <- paste(prefix, tag_base, sep = "_")

  rds_dir   <- file.path(outdir, "rds")
  truth_dir <- file.path(outdir, "truth")
  meta_dir  <- file.path(outdir, "meta")
  dir.create(rds_dir,   showWarnings = FALSE, recursive = TRUE)
  dir.create(truth_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(meta_dir,  showWarnings = FALSE, recursive = TRUE)

  ps_path <- file.path(rds_dir, paste0(tag, ".RDS"))
  saveRDS(PS, ps_path)

  truth <- tibble(
    taxon = colnames(Y),
    is_DA = colnames(Y) %in% taxa_names[da_idx],
    logFC = logFC[colnames(Y)]
  )
  truth_path <- file.path(truth_dir, paste0("truth_", tag_base, ".csv"))
  write.csv(truth, truth_path, row.names = FALSE)

  meta <- list(
    tag = tag,
    tag_base = tag_base,
    prefix = prefix,
    design = design_tag,
    p_DA = p_da,
    abs_logFC = abs_lfc,
    replicate = rep_id,
    n_ctrl = n_ctrl,
    n_case = n_case,
    batch_sd = batch_sd,
    seed_int = seed0 + scen_hash
  )
  meta_path <- file.path(meta_dir, paste0("meta_", tag_base, ".json"))
  write_json(meta, meta_path, auto_unbox = TRUE)

  otu_mat <- as(otu_table(PS), "matrix")
  # drop empty samples/taxa for safety
  otu_mat <- otu_mat[rowSums(otu_mat) > 0, colSums(otu_mat) > 0, drop = FALSE]
# 8 taxa distribution plot (rank–abundance) per simulated dataset
df_rank <- .rank_abundance_df(otu_mat)
plot_title <- sprintf("%s | %s | p_DA=%.2f | abs_logFC=%.1f | rep=%02d",
                      tag_base, design_tag, p_da, abs_lfc, rep_id)

p_rank <- .plot_rank_abundance(df_rank, title = plot_title)

plot_path <- file.path(meta_dir, paste0("taxa_rank_abundance_", tag_base, ".png"))

plot_path <- tryCatch(
  .save_rank_abundance_plot(p_rank, plot_path, width = 8, height = 5, dpi = 150),
  error = function(e) {
    # write a visible error log in meta_dir so it gets published too
    errfile <- file.path(meta_dir, "taxa_dist_plot_errors.log")
    writeLines(sprintf("[%s] %s", tag_base, conditionMessage(e)), errfile, sep = "\n", useBytes = TRUE)
    NA_character_
  }
)

if (!is.na(plot_path)) message("Saved taxa distribution plot: ", plot_path)
  tibble(
    tag = tag,
    tag_base = tag_base,
    design = design_tag,
    p_DA = p_da,  
    abs_logFC = abs_lfc,
    replicate = rep_id,
    rds = ps_path,
    truth_csv = truth_path,
    meta_json = meta_path,
    plot_taxa_dist = plot_path,
    ntaxa = nrow(otu_mat),
    nsamples = ncol(otu_mat),
    sparsity = 1 - mean(otu_mat > 0),
    seed_int = seed0 + scen_hash
  )
}

# --- high-level wrapper ---
simulate_midasim_from_phyloseq <- function(
  ps,
  outdir,
  n_ctrl = NULL, n_case = NULL,
  p_da_vec       = c(0.05, 0.10),
  abs_logfc_vec  = c(0.5, 1.0, 2.0),
  n_rep          = 3,
  design_vec     = c("balanced", "confounded"),
  subset_var     = NULL,
  subset_keep    = NULL,
  batch_col      = NULL,
  batch_sd       = 0.3,
  prefix         = "ps_midasim",
  seed0          = 123
) {
  stopifnot(inherits(ps, "phyloseq"))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  if (!is.null(subset_var) && !is.null(subset_keep)) {
    meta0 <- data.frame(sample_data(ps), check.names = FALSE)
    if (!subset_var %in% names(meta0)) {
      stop("subset_var '", subset_var, "' not found in sample_data.")
    }
    keep_ids <- rownames(meta0)[meta0[[subset_var]] %in% subset_keep]
    ps <- prune_samples(keep_ids, ps)
  }

  if (is.null(n_ctrl) || is.null(n_case)) {
    n_total <- phyloseq::nsamples(ps)
    n_ctrl  <- floor(n_total / 2)
    n_case  <- n_total - n_ctrl
  }

  M0 <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) M0 <- t(M0)
  meta <- data.frame(sample_data(ps), check.names = FALSE)
  tax_df <- if (!is.null(tax_table(ps))) as.matrix(tax_table(ps)) else NULL

  if (!is.null(batch_col) && batch_col %in% names(meta)) {
    batch_fac <- factor(meta[[batch_col]])
  } else {
    batch_fac <- factor(rep("Batch1", nrow(meta)))
  }
  batch_levels <- levels(batch_fac)

  fit <- MIDASim.setup(otu.tab = M0, mode = "parametric")

  lib_sizes <- rowSums(M0)
  lib_by_batch <- split(lib_sizes, batch_fac)

  idx <- list()
  k <- 0L
  for (design_tag in design_vec) {
    for (p_da in p_da_vec) {
      for (abs_lfc in abs_logfc_vec) {
        for (rep_id in seq_len(n_rep)) {
          k <- k + 1L
          idx[[k]] <- .simulate_one(
            fit = fit, base_mu_len = ncol(M0), tax_table_df = tax_df,
            batch_levels = batch_levels, batch_fac = batch_fac,
            lib_by_batch = lib_by_batch,
            n_ctrl = n_ctrl, n_case = n_case,
            p_da = p_da, abs_lfc = abs_lfc,
            design_tag = design_tag, rep_id = rep_id,
            batch_sd = batch_sd, outdir = outdir,
            prefix = prefix,
            seed0 = seed0
          )
        }
      }
    }
  }

  index_tbl <- bind_rows(idx)
  write_csv(index_tbl, file.path(outdir, "index.csv"))
  message("Wrote ", nrow(index_tbl), " scenarios to: ", file.path(outdir, "index.csv"))
  invisible(index_tbl)
}

# --------------------------
# CLI runner
# --------------------------
args <- commandArgs(trailingOnly = TRUE)

.get <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[[i + 1]])
  default
}

.as_num_vec <- function(x) {
  if (is.null(x) || is.na(x) || x == "") return(NULL)
  as.numeric(strsplit(x, ",")[[1]])
}

.as_chr_vec <- function(x) {
  if (is.null(x) || is.na(x) || x == "") return(NULL)
  trimws(strsplit(x, ",")[[1]])
}

in_rds  <- .get("--input")
outdir  <- .get("--outdir")
prefix  <- .get("--prefix", "ps_midasim")

if (is.null(in_rds) || is.null(outdir)) {
  stop("Usage:\n",
       "  Rscript run_simulation.R --input TEMPLATE.RDS --outdir OUTDIR [--prefix ps_midasim] ...\n")
}

ps <- readRDS(in_rds)

n_ctrl <- .get("--n_ctrl", NA_character_)
n_case <- .get("--n_case", NA_character_)
n_ctrl <- if (!is.na(n_ctrl)) as.integer(n_ctrl) else NULL
n_case <- if (!is.na(n_case)) as.integer(n_case) else NULL

p_da_vec      <- .as_num_vec(.get("--p_da", "0.05,0.10"))
abs_logfc_vec <- .as_num_vec(.get("--abs_logfc", "0.5,1.0,2.0"))
n_rep         <- as.integer(.get("--n_rep", "3"))
design_vec    <- .as_chr_vec(.get("--design", "balanced,confounded"))

subset_var    <- .get("--subset_var", NA_character_)
subset_keep   <- .get("--subset_keep", NA_character_)
subset_var  <- if (!is.na(subset_var) && subset_var != "") subset_var else NULL
subset_keep <- if (!is.na(subset_keep) && subset_keep != "") .as_chr_vec(subset_keep) else NULL

batch_col <- .get("--batch_col", NA_character_)
batch_col <- if (!is.na(batch_col) && batch_col != "") batch_col else NULL
batch_sd  <- as.numeric(.get("--batch_sd", "0.3"))

seed0    <- as.integer(.get("--seed0", "123"))

simulate_midasim_from_phyloseq(
  ps            = ps,
  outdir        = outdir,
  n_ctrl        = n_ctrl,
  n_case        = n_case,
  p_da_vec      = p_da_vec,
  abs_logfc_vec = abs_logfc_vec,
  n_rep         = n_rep,
  design_vec    = design_vec,
  subset_var    = subset_var,
  subset_keep   = subset_keep,
  batch_col     = batch_col,
  batch_sd      = batch_sd,
  prefix        = prefix,
  seed0         = seed0
)
