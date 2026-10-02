#!/usr/bin/env Rscript
# Dependencies: install.packages(c("glmnet", "ggplot2", "patchwork", "ragg"))
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("mixOmics", ask = FALSE, update = FALSE)
args <- commandArgs(trailingOnly = TRUE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
root <- dirname(normalizePath(script))
source(file.path(root, "benchmark.R"))
config <- list(seed = 20261001L, repeats = 3L, outer_folds = 4L,
               inner_folds = 3L, budgets = c(5L, 10L), components = c(2L, 3L),
               designs = c(0.1, 0.5), alpha = c(0, 0.5, 1),
               lambda = c(1, 0.1, 0.01))
# Custom input: Rscript run.R /path/to/data /path/to/results
# Three CSVs: samples.csv (sample_id, outcome, stratum), block_a.csv, block_b.csv.
# Omics CSVs contain sample_id followed by numeric, already transformed features.
# stratum controls fold balancing only; it is never used as a predictor.
if ("--prepare" %in% args) prepare_demo(file.path(root, "data"))
if ("--check" %in% args) { check_contracts(); quit(status = 0) }
input <- if (length(args) && !startsWith(args[1], "--")) args[1] else file.path(root, "data")
output <- if (length(args) > 1 && !startsWith(args[2], "--")) args[2] else file.path(root, "results")
run_benchmark(input, output, config)
