#!/usr/bin/env Rscript
# ======================================================================
# Apply the soft-consensus threshold learned from simulations to real
# DA results. Reads tool weights and optimal threshold from run_scoring.R
# outputs, then applies them to *_da.tsv files from real data.
# ======================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(fs)
})

args <- commandArgs(trailingOnly = TRUE)
.get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[[i + 1]])
  default
}

scoring_dir <- .get_arg("--scoring_dir")
da_dir      <- .get_arg("--da_dir")
outdir      <- .get_arg("--outdir", ".")
alpha       <- as.numeric(.get_arg("--alpha", "0.05"))
lfc_min     <- as.numeric(.get_arg("--lfc_min", "0"))

if (is.null(scoring_dir) || is.null(da_dir)) {
  stop("Usage: apply_soft_consensus.R --scoring_dir <dir> --da_dir <dir> ",
       "[--outdir .] [--alpha 0.05] [--lfc_min 0]")
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -------------------------
# Read scoring outputs
# -------------------------
scores_file <- file.path(scoring_dir, "soft_tool_scores.csv")
threshold_file <- file.path(scoring_dir, "optimal_threshold.txt")

if (!file.exists(scores_file)) stop("Missing: ", scores_file)
if (!file.exists(threshold_file)) stop("Missing: ", threshold_file)

scores_tbl <- read_csv(scores_file, show_col_types = FALSE)
T_opt <- as.numeric(readLines(threshold_file, n = 1))

cat("[apply_soft_consensus] Scoring dir  :", scoring_dir, "\n")
cat("[apply_soft_consensus] DA dir       :", da_dir, "\n")
cat("[apply_soft_consensus] Threshold    :", T_opt, "\n")
cat("[apply_soft_consensus] Alpha        :", alpha, "\n")
cat("[apply_soft_consensus] LFC min      :", lfc_min, "\n")

weight_map <- setNames(scores_tbl$score, tolower(scores_tbl$tool))
cat("[apply_soft_consensus] Tool weights :\n")
for (t in names(weight_map)) cat("  ", t, "=", round(weight_map[[t]], 4), "\n")

# -------------------------
# Read real DA results
# -------------------------
da_files <- fs::dir_ls(da_dir, regexp = "_da\\.tsv$", type = "file")
if (length(da_files) == 0) stop("No *_da.tsv files found in: ", da_dir)

read_da <- function(path) {
  df <- read_tsv(path, show_col_types = FALSE, col_types = cols(.default = col_character()))
  tool_name <- tolower(
    if ("tool" %in% names(df)) df$tool[1]
    else sub("_da\\.tsv$", "", basename(path))
  )
  df %>%
    transmute(
      taxon_id = taxon_id,
      lfc      = as.numeric(lfc),
      pvalue   = as.numeric(pvalue),
      adj_pval = as.numeric(adj_pval),
      tool     = tool_name
    )
}

da_long <- bind_rows(lapply(da_files, read_da))
cat("[apply_soft_consensus] Read", nrow(da_long), "rows from", length(da_files), "DA files\n")

# -------------------------
# Compute per-taxon weighted scores
# -------------------------
da_calls <- da_long %>%
  mutate(
    called = !is.na(adj_pval) & adj_pval < alpha & !is.na(lfc) & abs(lfc) >= lfc_min,
    sign = case_when(
      !called        ~ NA_integer_,
      lfc > 0        ~ 1L,
      lfc < 0        ~ -1L,
      TRUE           ~ NA_integer_
    ),
    weight = weight_map[tool]
  ) %>%
  filter(called, !is.na(sign), !is.na(weight), weight > 0)

all_tools <- sort(unique(da_calls$tool))

taxon_scores <- da_calls %>%
  group_by(taxon_id) %>%
  summarise(
    score_up = sum(weight[sign > 0]),
    score_dn = sum(weight[sign < 0]),
    n_up     = sum(sign > 0),
    n_dn     = sum(sign < 0),
    tools_up   = paste(sort(unique(tool[sign > 0])), collapse = ";"),
    tools_down = paste(sort(unique(tool[sign < 0])), collapse = ";"),
    tool_list  = list(tool),
    .groups  = "drop"
  ) %>%
  # Majority-direction: assign to the direction with higher weighted score.
  # Ties (equal weighted scores) are dropped as ambiguous.
  mutate(
    direction = case_when(
      score_up > score_dn ~ "up",
      score_dn > score_up ~ "down",
      TRUE                ~ "ambiguous"
    ),
    ensemble_score = pmax(score_up, score_dn),
    n_tools_agreeing = ifelse(direction == "up", n_up, n_dn)
  ) %>%
  filter(direction != "ambiguous", ensemble_score >= T_opt) %>%
  arrange(desc(ensemble_score), taxon_id)

# Add per-tool detection columns
for (t in all_tools) {
  taxon_scores[[paste0("det_", t)]] <- as.integer(
    sapply(taxon_scores$tool_list, function(tl) t %in% tl)
  )
}

taxon_scores <- taxon_scores %>%
  select(taxon_id, direction, ensemble_score, n_tools_agreeing,
         tools_up, tools_down, starts_with("det_"))

out_file <- file.path(outdir, "soft_consensus_results.csv")
write_csv(taxon_scores, out_file)

n_ambig <- sum(da_calls %>% group_by(taxon_id) %>%
  summarise(su = sum(weight[sign > 0]), sd = sum(weight[sign < 0]), .groups = "drop") %>%
  { .$su == .$sd & .$su > 0 })
if (n_ambig > 0)
  cat("[apply_soft_consensus]", n_ambig, "taxa dropped (equal weighted score in both directions)\n")
cat("[apply_soft_consensus] Taxa passing threshold:", nrow(taxon_scores), "\n")
cat("[apply_soft_consensus] Wrote:", out_file, "\n")

# Also write a summary
summary_tbl <- tibble(
  optimal_threshold = T_opt,
  alpha = alpha,
  lfc_min = lfc_min,
  n_taxa_passing = nrow(taxon_scores),
  n_up = sum(taxon_scores$direction == "up"),
  n_down = sum(taxon_scores$direction == "down")
)
write_csv(summary_tbl, file.path(outdir, "soft_consensus_summary.csv"))
