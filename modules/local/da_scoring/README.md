# DA Scoring / Benchmarking Outputs

This directory contains benchmarking results comparing differential abundance (DA) tool outputs against
simulated ground truth. Scenarios vary by:
- `design`: balanced vs confounded
- `p_DA`: proportion of truly DA taxa
- `abs_logFC`: effect size
- `replicate`: simulation replicate

All "byDesign" files are stratified by `design`.

## Core input summaries

### `index_tbl.csv`
Scenario index discovered from scenario folder names.
Columns: `tag`, `p_DA`, `abs_logFC`, `design`, `replicate`.

### `sim_taxa_unified.csv` / `sim_taxa_unified.rds`
Per-taxon unified table across all scenarios.
Includes:
- truth: `is_DA`, `logFC_true`
- scenario metadata: `tag`, `design`, `p_DA`, `abs_logFC`, `replicate`
- per-tool outputs: `p_<tool>`, `padj_<tool>`, `logFC_<tool>`

## Single-tool performance (pooled)

### `perf_tbl_singletools.csv`
Global pooled metrics per tool: TP/FP/FN/TN, precision, recall, F1, FDR.
Also includes F1-derived weights (`w_F1`).

### `coef_tbl_logitweights.csv` (+ optional `logit_fit.rds`)
Logistic-regression weights derived from -log10(p) features.
Outputs per tool: `beta`, `w_logit`.

### `weight_comp.csv`, `F1weight_vs_logitWeight.png`
Comparison scatter of F1-based vs logit-based weights.

## Consensus (k-of-n tools)

### `consensus_k_perScenario.csv`
Per-scenario metrics for each k.

### `consensus_k_summary.csv`, `consensus_k_summary.png`
Pooled mean metrics vs k.

### `consensus_k_summary_byDesign.csv`, `consensus_k_summary_byDesign.png`
Mean metrics vs k split by design (faceted).

### `union_intersection_summary.csv`
Union = k=1, Intersection = k=n_tools (pooled).

## Weighted consensus (threshold T)

Weighted consensus calls a taxon if:
`sum_i w_i * I(padj_i < alpha) >= T`

### F1 weights
- `wcons_F1_perScenario.csv`
- `wcons_F1_summary.csv`
- `wcons_F1_summary_byDesign.csv`
- `wconsensus_F1_summary.png`
- `wconsensus_F1_summary_byDesign.png`

### Logit weights
- `wcons_logit_perScenario.csv`
- `wcons_logit_summary.csv`
- `wcons_logit_summary_byDesign.csv`

### Comparison plots
Pooled:
- `weighttype_F1vsT.png`
- `weighttype_FDRvsT.png`

By design:
- `weighttype_F1vsT_byDesign.png`
- `weighttype_FDRvsT_byDesign.png`

## Per-tool performance vs effect size (abs_logFC)

Tables:
- `tool_perScenario_perf.csv`
- `tool_F1_FDR_by_absLFC.csv`
- `tool_F1_FDR_by_absLFC_byDesign.csv`

Plots:
- `perTool_F1_lfc.png`, `perTool_FDR_lfc.png`
- `perTool_F1_lfc_byDesign.png`, `perTool_FDR_lfc_byDesign.png`

## Best-threshold selection + "with consensus" curves

### `chosen_thresholds.csv`
Selected thresholds:
- `k_star`
- `T_star_F1`
- `T_star_logit`

Selection controlled by:
- `alpha`
- `target_FDR` (if strategy uses it)
- `threshold_strategy`

### With-consensus performance vs abs_logFC (pooled)
- `lfc_F1_FDR_with_consensus_table.csv`
- `lfc_F1_with_consensus.png`
- `lfc_FDR_with_consensus.png`

### With-consensus performance vs abs_logFC (by design)
- `lfc_F1_FDR_with_consensus_table_byDesign.csv`
- `lfc_F1_with_consensus_byDesign.png`
- `lfc_FDR_with_consensus_byDesign.png`

## Global method comparison (pooled)

### `global_method_compare.csv`
Pooled metrics for:
- each single tool
- best k-of-n consensus (k_star)
- best weighted consensus (F1/logit)

### `globalPerformance_FDRvsF1.png`
FDR vs F1 scatter for all pooled methods.
