#!/usr/bin/env Rscript
# ======================================================================
# Benchmark simulated DA results against truth (weighted consensus scoring)
#
# What it does:
# 1) Reads all simulated tool results and truth tables
# 2) Computes TP / FP / FN / TN, FDR and F1 for each tool and dataset
# 3) Summarizes methods across scenarios
# 4) Selects the best single tool under FDR < fdr_target
# 5) Computes weighted tool scores based on mean F1 and mean FDR
# 6) Finds the optimal weighted consensus threshold on simulations
# 7) Saves CSVs, plot, optimal_threshold.txt, and RDS
#
# Assumptions:
# - tool result files have columns: taxon_id, lfc, pvalue, adj_pval, tool
# - truth files have columns: taxon, is_DA, logFC
# ======================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
})

# ----------------------------------------------------------------------
# CLI ARGUMENT PARSING
# ----------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
.get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[[i + 1]])
  default
}

base_results     <- .get_arg("--base_results")
base_truth       <- .get_arg("--base_truth")
outdir           <- .get_arg("--outdir", "scoring_output")
tools_csv        <- .get_arg("--tools", "adapt,corncob,linda,locom,maaslin2")
alpha_results    <- as.numeric(.get_arg("--alpha", "0.05"))
fdr_target       <- as.numeric(.get_arg("--fdr_target", "0.05"))
fdr_reference    <- as.numeric(.get_arg("--fdr_reference", "0.10"))
retention_at_ref <- as.numeric(.get_arg("--retention_at_ref", "0.25"))
prefix_filter    <- .get_arg("--prefix", NULL)

if (is.null(base_results) || is.null(base_truth)) {
  stop("Usage: run_scoring.R --base_results <dir> --base_truth <dir> ",
       "[--outdir scoring_output] [--tools adapt,corncob,...] ",
       "[--alpha 0.05] [--fdr_target 0.05] [--fdr_reference 0.10] ",
       "[--retention_at_ref 0.25] [--prefix BV]")
}

tools <- trimws(unlist(strsplit(tools_csv, ",")))
tool_files <- setNames(paste0(tools, "_da.tsv"), tools)

round_digits       <- 8

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("[run_scoring] base_results :", base_results, "\n")
cat("[run_scoring] base_truth   :", base_truth, "\n")
cat("[run_scoring] outdir       :", outdir, "\n")
cat("[run_scoring] tools        :", paste(tools, collapse = ", "), "\n")
cat("[run_scoring] alpha        :", alpha_results, "\n")
cat("[run_scoring] fdr_target   :", fdr_target, "\n")
cat("[run_scoring] fdr_reference:", fdr_reference, "\n")
cat("[run_scoring] retention    :", retention_at_ref, "\n")

# ----------------------------------------------------------------------
# DYNAMIC SCENARIO DISCOVERY
# ----------------------------------------------------------------------
discover_scenarios <- function(sim_root, prefix_filter = NULL) {
  all_dirs <- list.dirs(sim_root, full.names = FALSE, recursive = FALSE)
  pat <- "^(.*?)_?pda(\\d{2})_logfc([0-9.]+)_([^_]+)_rep0?(\\d+)$"
  m <- stringr::str_match(all_dirs, pat)
  ok <- !is.na(m[, 1])

  if (!any(ok)) {
    stop("No scenario directories matched expected pattern in: ", sim_root, "\n",
         "Expected pattern like: PREFIX_pda05_logfc0.5_balanced_rep01")
  }

  tbl <- tibble(
    folder  = all_dirs[ok],
    prefix  = m[ok, 2],
    pda     = m[ok, 3],
    logfc   = m[ok, 4],
    design  = m[ok, 5],
    rep     = as.integer(m[ok, 6])
  )

  if (!is.null(prefix_filter) && nzchar(prefix_filter)) {
    tbl <- tbl %>% filter(prefix == prefix_filter)
    if (nrow(tbl) == 0) stop("No scenarios found with prefix '", prefix_filter, "'")
  }

  tbl$scenario <- paste0("pda", tbl$pda, "_logfc", tbl$logfc, "_", tbl$design)
  tbl <- tbl %>% arrange(scenario, rep)
  tbl
}

scenarios <- discover_scenarios(base_results, prefix_filter)
dataset_ids <- scenarios$folder
cat("[run_scoring] Discovered", nrow(scenarios), "scenario folders\n")
cat("[run_scoring] Unique scenarios:", length(unique(scenarios$scenario)), "\n")

# ----------------------------------------------------------------------
# HELPERS
# ----------------------------------------------------------------------
direction_from_lfc <- function(x) {
  ifelse(is.na(x), NA_character_,
         ifelse(x > 0, "higher_in_disease",
                ifelse(x < 0, "higher_in_control", "none")))
}

read_sim_result <- function(tool, folder, base_results, tool_files) {
  fpath <- file.path(base_results, folder, tool_files[[tool]])
  if (!file.exists(fpath)) {
    warning("Missing result file: ", fpath)
    return(data.frame(taxon_id = character(0), lfc = numeric(0),
                      pvalue = numeric(0), adj_pval = numeric(0),
                      tool = character(0)))
  }
  read.table(fpath, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
}

read_sim_truth <- function(folder, base_truth) {
  base_suffix <- stringr::str_extract(folder, "pda\\d{2}_logfc[0-9.]+_[^_]+_rep0?\\d+")
  if (is.na(base_suffix)) stop("Cannot extract scenario suffix from folder: ", folder)
  fpath <- file.path(base_truth, paste0("truth_", base_suffix, ".csv"))
  if (!file.exists(fpath)) stop("Missing truth file: ", fpath)
  read.csv(fpath, stringsAsFactors = FALSE)
}

score_directional_prediction <- function(pred_df, truth_df) {
  merged <- truth_df %>%
    transmute(
      taxon = taxon,
      true_DA = is_DA,
      true_direction = ifelse(is_DA, direction_from_lfc(logFC), "none")
    ) %>%
    left_join(pred_df, by = "taxon") %>%
    mutate(
      pred_DA = ifelse(is.na(pred_DA), FALSE, pred_DA),
      pred_direction = ifelse(is.na(pred_direction), "none", pred_direction)
    )

  correct_direction <- merged$pred_direction == merged$true_direction

  TP <- sum(merged$pred_DA & merged$true_DA & correct_direction)
  FP <- sum(merged$pred_DA & (!merged$true_DA | !correct_direction))
  FN <- sum(merged$true_DA & (!merged$pred_DA | !correct_direction))
  TN <- sum(!merged$pred_DA & !merged$true_DA)

  FDR <- if ((TP + FP) == 0) 0 else FP / (TP + FP)
  F1  <- if ((2 * TP + FP + FN) == 0) 0 else 2 * TP / (2 * TP + FP + FN)

  data.frame(TP = TP, FP = FP, FN = FN, TN = TN,
             discoveries = TP + FP, FDR = FDR, F1 = F1)
}

score_single_tool_dataset <- function(res_df, truth_df, alpha = 0.05,
                                       significance_col = "adj_pval") {
  pred_df <- res_df %>%
    transmute(
      taxon = taxon_id,
      pred_DA = .data[[significance_col]] <= alpha,
      pred_direction = direction_from_lfc(lfc)
    )
  score_directional_prediction(pred_df, truth_df)
}

compute_weighted_tool_scores <- function(tool_summary_table,
                                      fdr_target = 0.05,
                                      fdr_reference = 0.10,
                                      retention_at_reference = 0.25) {
  if (fdr_reference <= fdr_target) stop("fdr_reference must be > fdr_target")
  if (retention_at_reference <= 0 || retention_at_reference >= 1) {
    stop("retention_at_reference must be in (0,1)")
  }

  df <- tool_summary_table %>% select(tool, mean_F1, mean_FDR)
  lambda <- -log(retention_at_reference) / (fdr_reference - fdr_target)

  df <- df %>%
    mutate(
      penalty_factor = ifelse(mean_FDR <= fdr_target, 1,
                              exp(-lambda * (mean_FDR - fdr_target))),
      unnormalized_weight = mean_F1 * penalty_factor
    )

  if (sum(df$unnormalized_weight, na.rm = TRUE) == 0) {
    df$score <- 1 / nrow(df)
  } else {
    df$score <- df$unnormalized_weight / sum(df$unnormalized_weight, na.rm = TRUE)
  }

  list(lambda = lambda, summary_table = df %>% arrange(desc(score)))
}

build_weighted_ensemble_dataset <- function(result_list_all, truth_df, tool_scores,
                                         alpha = 0.05,
                                         significance_col = "adj_pval",
                                         round_digits = 8) {
  use_tools <- intersect(names(tool_scores), names(result_list_all))
  tool_scores <- tool_scores[use_tools]

  ensemble_df <- truth_df %>%
    transmute(
      taxon = taxon,
      true_DA = is_DA,
      true_logFC = logFC,
      true_direction = ifelse(is_DA, direction_from_lfc(logFC), "none")
    )

  for (tool in use_tools) {
    res_small <- result_list_all[[tool]] %>%
      transmute(
        taxon = taxon_id,
        !!paste0("called_pos_", tool) := (.data[[significance_col]] <= alpha) & (lfc > 0),
        !!paste0("called_neg_", tool) := (.data[[significance_col]] <= alpha) & (lfc < 0)
      )
    ensemble_df <- ensemble_df %>% left_join(res_small, by = "taxon")
  }

  pos_cols <- paste0("called_pos_", use_tools)
  neg_cols <- paste0("called_neg_", use_tools)
  ensemble_df[pos_cols] <- lapply(ensemble_df[pos_cols], function(x) ifelse(is.na(x), FALSE, x))
  ensemble_df[neg_cols] <- lapply(ensemble_df[neg_cols], function(x) ifelse(is.na(x), FALSE, x))

  pos_mat <- as.matrix(ensemble_df[, pos_cols, drop = FALSE]); storage.mode(pos_mat) <- "numeric"
  neg_mat <- as.matrix(ensemble_df[, neg_cols, drop = FALSE]); storage.mode(neg_mat) <- "numeric"

  ensemble_df$positive_support <- round(rowSums(sweep(pos_mat, 2, tool_scores, `*`)), round_digits)
  ensemble_df$negative_support <- round(rowSums(sweep(neg_mat, 2, tool_scores, `*`)), round_digits)
  ensemble_df$ensemble_score   <- pmax(ensemble_df$positive_support, ensemble_df$negative_support)
  ensemble_df$pred_direction   <- ifelse(
    ensemble_df$positive_support == 0 & ensemble_df$negative_support == 0,
    "none",
    ifelse(
      ensemble_df$positive_support > ensemble_df$negative_support,
      "higher_in_disease",
      ifelse(
        ensemble_df$negative_support > ensemble_df$positive_support,
        "higher_in_control",
        "ambiguous"
      )
    )
  )

  ensemble_df
}

find_optimal_weighted_threshold <- function(raw_results, truth_tables, tool_scores,
                                         dataset_ids, alpha = 0.05,
                                         significance_col = "adj_pval",
                                         fdr_target = 0.05,
                                         round_digits = 8) {
  ensemble_tables <- lapply(dataset_ids, function(ds) {
    build_weighted_ensemble_dataset(
      result_list_all = lapply(names(raw_results), function(tool) raw_results[[tool]][[ds]]) |> setNames(names(raw_results)),
      truth_df = truth_tables[[ds]],
      tool_scores = tool_scores,
      alpha = alpha, significance_col = significance_col, round_digits = round_digits
    )
  })
  names(ensemble_tables) <- dataset_ids

  candidate_thresholds <- sort(unique(unlist(lapply(ensemble_tables, function(x) x$ensemble_score))))
  candidate_thresholds <- candidate_thresholds[candidate_thresholds > 0]

  threshold_results <- bind_rows(lapply(dataset_ids, function(ds) {
    bind_rows(lapply(candidate_thresholds, function(thr) {
      pred_df <- ensemble_tables[[ds]] %>%
        transmute(
          taxon = taxon,
          pred_DA = (ensemble_score >= thr) & pred_direction %in% c("higher_in_disease", "higher_in_control"),
          pred_direction = pred_direction
        )
      out <- score_directional_prediction(pred_df, truth_tables[[ds]])
      out$dataset <- ds
      out$threshold <- thr
      out
    }))
  }))

  threshold_summary <- threshold_results %>%
    group_by(threshold) %>%
    summarise(
      mean_FDR = mean(FDR, na.rm = TRUE),
      mean_F1 = mean(F1, na.rm = TRUE),
      sd_FDR = sd(FDR, na.rm = TRUE),
      sd_F1 = sd(F1, na.rm = TRUE),
      mean_discoveries = mean(discoveries, na.rm = TRUE),
      .groups = "drop"
    )

  valid <- threshold_summary %>% filter(mean_FDR < fdr_target)

  if (nrow(valid) == 0) {
    warning("No weighted threshold achieved mean FDR < ", fdr_target,
            ". Using fallback: minimum mean FDR, then maximum mean F1, then higher threshold.")
    best_row <- threshold_summary %>%
      arrange(mean_FDR, desc(mean_F1), desc(threshold)) %>%
      dplyr::slice_head(n = 1)
  } else {
    best_row <- valid %>%
      arrange(desc(mean_F1), mean_FDR, desc(threshold)) %>%
      dplyr::slice_head(n = 1)
  }

  list(
    ensemble_tables = ensemble_tables,
    threshold_results = threshold_results,
    threshold_summary = threshold_summary,
    best_threshold = best_row$threshold,
    best_row = best_row
  )
}

# ----------------------------------------------------------------------
# PLOTTING
# ----------------------------------------------------------------------
plot_global_performance <- function(tool_summary_global, best_threshold_row,
                                    fdr_target = 0.05, outdir = ".") {
  tools_df <- tool_summary_global %>%
    transmute(label = tool, method_type = "Single tool", mean_FDR = mean_FDR, mean_F1 = mean_F1)
  thr_df <- best_threshold_row %>%
    transmute(label = paste0("Weighted T = ", round(threshold, 4)), method_type = "Weighted threshold",
              mean_FDR = mean_FDR, mean_F1 = mean_F1)

  plot_df <- bind_rows(tools_df, thr_df)

  p <- ggplot(plot_df, aes(x = mean_FDR, y = mean_F1, color = method_type, label = label)) +
    geom_vline(xintercept = fdr_target, linetype = "dashed", color = "grey40") +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(size = 3.5, show.legend = FALSE,
                             box.padding = 0.4, point.padding = 0.3, max.overlaps = Inf) +
    scale_x_reverse(name = "Mean FDR", limits = c(1, 0)) +
    scale_y_continuous(name = "Mean F1", limits = c(0, 1)) +
    labs(title = "Global performance: individual tools vs weighted consensus threshold",
         subtitle = paste0("Dashed line = FDR target (", fdr_target, ")"), color = NULL) +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

  ggsave(file.path(outdir, "global_performance_tools_threshold.png"), p, width = 8, height = 5, dpi = 300)
  p
}

# ======================================================================
# MAIN ANALYSIS
# ======================================================================

# 1) Read all simulated results and truth tables
cat("[run_scoring] Reading truth tables...\n")
truth_tables <- setNames(
  lapply(dataset_ids, read_sim_truth, base_truth = base_truth),
  dataset_ids
)

cat("[run_scoring] Reading tool results...\n")
raw_results <- lapply(tools, function(tool) {
  setNames(
    lapply(dataset_ids, function(ds) {
      read_sim_result(tool, ds, base_results = base_results, tool_files = tool_files)
    }),
    dataset_ids
  )
})
names(raw_results) <- tools

# 2) Score each tool on each dataset
cat("[run_scoring] Scoring individual tools...\n")
tool_metrics_by_dataset <- bind_rows(lapply(tools, function(tool) {
  bind_rows(lapply(dataset_ids, function(ds) {
    out <- score_single_tool_dataset(
      res_df = raw_results[[tool]][[ds]],
      truth_df = truth_tables[[ds]],
      alpha = alpha_results, significance_col = "adj_pval"
    )
    out$tool <- tool
    out$dataset <- ds
    scen_info <- scenarios %>% filter(folder == ds)
    out$scenario <- scen_info$scenario[1]
    out$replicate <- scen_info$rep[1]
    out
  }))
}))

# 3) Summarize by scenario and globally
tool_summary_by_scenario <- tool_metrics_by_dataset %>%
  group_by(tool, scenario) %>%
  summarise(mean_FDR = mean(FDR, na.rm = TRUE),
            mean_F1 = mean(F1, na.rm = TRUE), .groups = "drop")

tool_summary_global <- tool_metrics_by_dataset %>%
  group_by(tool) %>%
  summarise(mean_FDR = mean(FDR, na.rm = TRUE),
            mean_F1 = mean(F1, na.rm = TRUE), .groups = "drop") %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  mutate(tool = as.character(tool)) %>%
  tibble::as_tibble()

F1_table <- tool_summary_by_scenario %>%
  select(tool, scenario, mean_F1) %>%
  pivot_wider(names_from = tool, values_from = mean_F1) %>%
  as.data.frame()
rownames(F1_table) <- F1_table$scenario
F1_table <- F1_table[, setdiff(colnames(F1_table), "scenario"), drop = FALSE]
F1_table <- rbind(F1_table, Average = tool_summary_global$mean_F1)
colnames(F1_table) <- tool_summary_global$tool

FDR_table <- tool_summary_by_scenario %>%
  select(tool, scenario, mean_FDR) %>%
  pivot_wider(names_from = tool, values_from = mean_FDR) %>%
  as.data.frame()
rownames(FDR_table) <- FDR_table$scenario
FDR_table <- FDR_table[, setdiff(colnames(FDR_table), "scenario"), drop = FALSE]
FDR_table <- rbind(FDR_table, Average = tool_summary_global$mean_FDR)
colnames(FDR_table) <- tool_summary_global$tool

# 4) Best single tool under FDR < fdr_target
best_single_tool <- tool_summary_global %>%
  filter(mean_FDR < fdr_target) %>%
  arrange(desc(mean_F1), mean_FDR, tool)

if (nrow(best_single_tool) == 0) {
  warning("No single tool achieved mean FDR < ", fdr_target,
          ". Using fallback: minimum mean FDR, then maximum mean F1.")
  best_single_tool <- tool_summary_global %>%
    arrange(mean_FDR, desc(mean_F1), tool) %>%
    dplyr::slice_head(n = 1)
} else {
  best_single_tool <- best_single_tool %>% dplyr::slice_head(n = 1)
}

# 5) Weighted tool scores and optimal threshold
cat("[run_scoring] Computing weighted tool scores...\n")
weighted_scores <- compute_weighted_tool_scores(
  tool_summary_table = tool_summary_global,
  fdr_target = fdr_target,
  fdr_reference = fdr_reference,
  retention_at_reference = retention_at_ref
)

weighted_tool_scores <- weighted_scores$summary_table$score
names(weighted_tool_scores) <- weighted_scores$summary_table$tool

cat("[run_scoring] Finding optimal weighted threshold...\n")
weighted_threshold <- find_optimal_weighted_threshold(
  raw_results = raw_results,
  truth_tables = truth_tables,
  tool_scores = weighted_tool_scores,
  dataset_ids = dataset_ids,
  alpha = alpha_results, significance_col = "adj_pval",
  fdr_target = fdr_target, round_digits = round_digits
)

# 6) Plot
cat("[run_scoring] Generating plot...\n")
tryCatch({
  plot_global_performance(tool_summary_global, weighted_threshold$best_row,
                          fdr_target = fdr_target, outdir = outdir)
}, error = function(e) {
  warning("Plot generation failed (non-fatal): ", conditionMessage(e))
})

# 7) Save outputs
cat("[run_scoring] Saving outputs...\n")
simulation_analysis <- list(
  settings = list(
    alpha_results = alpha_results, fdr_target = fdr_target,
    fdr_reference = fdr_reference, retention_at_ref = retention_at_ref,
    round_digits = round_digits, tools = tools
  ),
  scenarios = scenarios,
  dataset_ids = dataset_ids,
  tool_metrics_by_dataset = tool_metrics_by_dataset,
  tool_summary_by_scenario = tool_summary_by_scenario,
  tool_summary_global = tool_summary_global,
  F1_table = F1_table, FDR_table = FDR_table,
  best_single_tool = best_single_tool,
  weighted_scores = weighted_scores,
  weighted_threshold = weighted_threshold
)

saveRDS(simulation_analysis, file.path(outdir, "simulation_analysis.rds"))
write.csv(tool_metrics_by_dataset, file.path(outdir, "tool_metrics_by_dataset.csv"), row.names = FALSE)
write.csv(tool_summary_by_scenario, file.path(outdir, "tool_summary_by_scenario.csv"), row.names = FALSE)
write.csv(tool_summary_global, file.path(outdir, "tool_summary_global.csv"), row.names = FALSE)
write.csv(best_single_tool, file.path(outdir, "best_single_tool.csv"), row.names = FALSE)
write.csv(weighted_scores$summary_table, file.path(outdir, "tool_scores.csv"), row.names = FALSE)
write.csv(weighted_threshold$threshold_summary, file.path(outdir, "threshold_summary.csv"), row.names = FALSE)

writeLines(as.character(weighted_threshold$best_threshold),
           file.path(outdir, "optimal_threshold.txt"))

cat("[run_scoring] Optimal weighted threshold:", weighted_threshold$best_threshold, "\n")
cat("[run_scoring] Best single tool:", best_single_tool$tool[1], "\n")
cat("[run_scoring] Done. Results in:", outdir, "\n")
