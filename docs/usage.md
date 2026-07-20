# nf-core/daamicrobiome: Usage

## Input

The pipeline requires a **phyloseq RDS object** (`.RDS`) containing:

- An OTU/ASV table (taxa as rows, samples as columns)
- A sample metadata data frame with at least a condition/group column
- A taxonomy table

```bash
--input /path/to/phyloseq.RDS
```

### Optional: pre-extracted control phyloseq (Path B only)

If you have already subset the control/healthy samples into a separate phyloseq object, provide it with `--input_control`. Otherwise, the pipeline extracts them automatically using `--condition` and `--base_level`.

```bash
--input_control /path/to/control_phyloseq.RDS
```

## Execution paths

### Path A: k-intersection consensus (`--simulate false`)

Runs DA tools on the real dataset and computes a k-intersection consensus.

```bash
nextflow run nf-core/daamicrobiome \
  -profile docker \
  --input phyloseq.RDS \
  --simulate false \
  --condition "study_condition" \
  --base_level "healthy" \
  --outdir results/
```

### Path B: simulation-trained weighted consensus (`--simulate true`)

Default mode. Extracts control samples, simulates ground-truth datasets via MIDASim, scores DA tools, and applies the optimal threshold to real data.

```bash
nextflow run nf-core/daamicrobiome \
  -profile docker \
  --input phyloseq.RDS \
  --simulate true \
  --condition "study_condition" \
  --base_level "healthy" \
  --outdir results/
```

## Parameters

### Core parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | (required) | Path to the full phyloseq RDS object |
| `--input_control` | `null` | Path to a control-only phyloseq RDS (optional, Path B) |
| `--outdir` | `results/` | Output directory |
| `--simulate` | `true` | `false` = Path A (k-intersection), `true` = Path B (weighted consensus) |

### Metadata parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--condition` | `study_condition` | Column in sample_data identifying the group variable |
| `--base_level` | `healthy` | Reference/control level in the condition column |
| `--confounder` | `null` | Optional confounder column(s) for real data DA tools |
| `--bc_confounder` | `null` | Batch-correction confounder (MaAsLin2 only) |

### DA tool toggles

All tools are enabled by default except metagenomeSeq.

| Parameter | Default |
|-----------|---------|
| `--tool_adapt` | `true` |
| `--tool_corncob` | `true` |
| `--tool_linda` | `true` |
| `--tool_locom` | `true` |
| `--tool_maaslin2` | `true` |
| `--tool_metagenomeseq` | `false` |

### Simulation parameters (Path B only)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--sim_n_ctrl` | `null` (auto) | Number of control samples per simulation |
| `--sim_n_case` | `null` (auto) | Number of case samples per simulation |
| `--sim_p_da` | `[0.05, 0.10]` | Proportion of DA taxa |
| `--sim_abs_logfc` | `[0.5, 1.0, 2.0]` | Absolute log-fold change values |
| `--sim_n_rep` | `3` | Number of replicates per scenario |
| `--sim_design` | `["balanced"]` | Simulation design |
| `--sim_seed0` | `123` | Random seed for reproducibility |

### Scoring parameters (Path B only)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--eval_alpha` | `0.05` | Significance threshold |
| `--eval_fdr_target` | `0.05` | Target FDR for threshold optimisation |
| `--eval_fdr_reference` | `0.10` | Reference FDR level |
| `--eval_retention_at_ref` | `0.25` | Minimum retention at reference FDR |

### Path A parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--consensus_alpha` | `0.05` | Significance threshold for k-intersection |
| `--consensus_lfc_min` | `0` | Minimum absolute log-fold change filter |

### Resource limits

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--max_cpus` | `16` | Maximum CPUs per process |
| `--max_memory` | `128.GB` | Maximum memory per process |
| `--max_time` | `240.h` | Maximum wall time per process |

## Profiles

| Profile | Description |
|---------|-------------|
| `docker` | Run with Docker containers |
| `hpc` | Run with Singularity on a SLURM cluster |
