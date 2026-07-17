#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(phyloseq)
  library(optparse)
})

option_list <- list(
  make_option(c("-i", "--input"),     type = "character", help = "Input phyloseq RDS (full dataset with case + control)"),
  make_option(c("-o", "--output"),    type = "character", help = "Output RDS path for control-only phyloseq"),
  make_option(c("-c", "--condition"), type = "character", help = "Sample data column containing the group variable"),
  make_option(c("-r", "--reference"), type = "character", help = "Value in the condition column representing control/healthy samples")
)

args <- parse_args(OptionParser(option_list = option_list))

if (is.null(args$input) || is.null(args$output) ||
    is.null(args$condition) || is.null(args$reference)) {
  stop("All arguments are required: --input, --output, --condition, --reference")
}

ps <- readRDS(args$input)
if (!inherits(ps, "phyloseq")) stop("Input is not a phyloseq object.")

sdata <- as(sample_data(ps), "data.frame")
if (!(args$condition %in% colnames(sdata))) {
  stop("Column '", args$condition, "' not found in sample_data. ",
       "Available columns: ", paste(colnames(sdata), collapse = ", "))
}

vals <- as.character(sdata[[args$condition]])
if (!(args$reference %in% vals)) {
  stop("Reference level '", args$reference, "' not found in column '", args$condition, "'. ",
       "Available levels: ", paste(sort(unique(vals)), collapse = ", "))
}

keep <- sample_names(ps)[vals == args$reference]
ps_ctrl <- prune_samples(keep, ps)
ps_ctrl <- prune_taxa(taxa_sums(ps_ctrl) > 0, ps_ctrl)

cat(sprintf("Total samples in input  : %d\n", nsamples(ps)))
cat(sprintf("Control samples retained: %d\n", nsamples(ps_ctrl)))
cat(sprintf("Taxa retained (non-zero): %d (of %d)\n", ntaxa(ps_ctrl), ntaxa(ps)))

saveRDS(ps_ctrl, args$output)
cat("Saved control-only phyloseq to:", args$output, "\n")
