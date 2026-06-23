#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(UpSetR)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grep(file_arg, args)])
script_dir <- if (length(script_path) == 0) getwd() else dirname(normalizePath(script_path[1]))
setwd(script_dir)

results_dir <- "ipf_analysis_results"
sens_dir <- "ipf_sensitivity_results"
figures_dir <- file.path(results_dir, "figures")
tables_dir <- file.path(results_dir, "tables")

sig_set <- function(path) {
  df <- read.csv(path, stringsAsFactors = FALSE)
  df$gene_id[df$adj.P.Val < 0.05 & abs(df$logFC) >= 1]
}

sets <- list(
  Baseline = sig_set(file.path(results_dir, "tables", "de_disease_average.csv")),
  No_age = sig_set(file.path(sens_dir, "tables", "de_disease_average_no_age.csv")),
  Exclude_GSE184316 = sig_set(file.path(sens_dir, "tables", "de_disease_average_exclude_GSE184316.csv")),
  Leave_out_GSE213001 = sig_set(file.path(sens_dir, "tables", "de_disease_average_leave_out_GSE213001.csv"))
)

membership <- fromList(sets)
write.csv(membership, file.path(tables_dir, "disease_average_overlap_upset_membership.csv"), row.names = FALSE)

pdf(file.path(figures_dir, "upset_disease_average_sensitivity.pdf"), width = 10.6, height = 7.1, useDingbats = FALSE)
upset(
  membership,
  sets = c("Baseline", "No_age", "Exclude_GSE184316", "Leave_out_GSE213001"),
  order.by = c("freq", "degree"),
  keep.order = TRUE,
  sets.bar.color = c("#4D4D4D", "#D55E00", "#009E73", "#0072B2"),
  main.bar.color = "#4D4D4D",
  matrix.color = "#4D4D4D",
  mainbar.y.label = "Intersection size",
  sets.x.label = "Genes in each set",
  text.scale = 1.3
)
dev.off()
