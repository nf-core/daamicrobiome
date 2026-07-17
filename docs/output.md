# nf-core/daamicrobiome: Output

## Path A output (`--simulate false`)

### k-intersection consensus

Output directory: `results/k_intersection/`

| File | Description |
|------|-------------|
| `k_intersection.xlsx` | XLSX workbook with one sheet per k (k = 1 to n_tools). Each sheet lists taxa called DA by at least k tools in the same direction. |
| `k_<N>_intersection.csv` | Individual CSV files for each k value |

Columns in each sheet/CSV:

- `taxon` — taxon identifier
- `direction` — majority direction (`up` or `down`)
- `n_tools` — number of tools calling this taxon DA
- `tools_up` — tools calling this taxon DA in the up direction
- `tools_down` — tools calling this taxon DA in the down direction

### DA tool results (real)

Output directory: `results/da_tools_real/<dataset_id>/`

Each enabled DA tool produces a TSV file with per-taxon p-values, adjusted p-values, and log-fold change estimates.

## Path B output (`--simulate true`)

### Soft consensus results

Output directory: `results/weighted_consensus/`

| File | Description |
|------|-------------|
| `weighted_consensus_da_taxa.csv` | Final list of consensus DA taxa with direction and per-tool details |

Columns:

- `taxon` — taxon identifier
- `direction` — consensus direction (`up` or `down`)
- `weighted_score` — weighted score across tools
- `tools_up` — tools calling this taxon up
- `tools_down` — tools calling this taxon down
- `det_<tool>` — per-tool detection indicator (1/0)

### Scoring results

Output directory: `results/da_scoring/`

| File | Description |
|------|-------------|
| `soft_tool_scores.csv` | Per-tool weights learned from simulations |
| `optimal_threshold.txt` | The optimal score threshold |
| `soft_threshold_summary.csv` | Summary of threshold optimisation |
| `global_performance_tools_threshold.png` | Performance plot across tools and thresholds |
| `tool_metrics_by_dataset.csv` | Per-tool, per-simulation-scenario metrics |
| `simulation_analysis.rds` | Full R object with all scoring results |

### Simulated data

Output directory: `results/simulated_data/`

| File/Directory | Description |
|----------------|-------------|
| `rds/*.RDS` | Simulated phyloseq objects |
| `truth/` | Ground-truth DA status per simulation scenario |
| `meta/` | Simulation metadata |
| `index.csv` | Index mapping scenario IDs to simulation parameters |

### DA tool results (simulated)

Output directory: `results/da_tools/<scenario_id>/`

Per-tool TSV results for each simulated scenario, used by the scoring step.

### DA tool results (real)

Output directory: `results/da_tools_real/<dataset_id>/`

Same as Path A — per-tool TSV files for the real dataset.

### Control extraction

Output directory: `results/control_extract/` (only if `--input_control` was not provided)

| File | Description |
|------|-------------|
| `control_phyloseq.RDS` | Phyloseq object containing only control/healthy samples |
