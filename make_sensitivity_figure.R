#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grep(file_arg, args)])
script_dir <- if (length(script_path) == 0) getwd() else dirname(normalizePath(script_path[1]))
setwd(script_dir)

sens <- read.csv("ipf_sensitivity_results/tables/sensitivity_vs_baseline.csv", stringsAsFactors = FALSE)
summ <- read.csv("ipf_sensitivity_results/tables/sensitivity_summary.csv", stringsAsFactors = FALSE)

palette_publication <- c("Primary sensitivity" = "#D55E00", "Leave-one-study-out" = "#0072B2")

if ("contrast" %in% names(sens)) {
  sens <- sens[sens$contrast == "disease_average", , drop = FALSE]
}

plot_df <- merge(sens, summ[, c("analysis", "samples")], by = "analysis", all.x = TRUE)
plot_df$analysis_label <- plot_df$analysis
plot_df$analysis_label[plot_df$analysis_label == "exclude_GSE184316"] <- "Exclude GSE184316"
plot_df$analysis_label[plot_df$analysis_label == "no_age"] <- "No age covariate"
plot_df$analysis_label[plot_df$analysis_label == "leave_out_GSE138239"] <- "Leave out GSE138239"
plot_df$analysis_label[plot_df$analysis_label == "leave_out_GSE184316"] <- "Leave out GSE184316"
plot_df$analysis_label[plot_df$analysis_label == "leave_out_GSE199949"] <- "Leave out GSE199949"
plot_df$analysis_label[plot_df$analysis_label == "leave_out_GSE213001"] <- "Leave out GSE213001"
plot_df$analysis_label[plot_df$analysis_label == "leave_out_GSE73189"] <- "Leave out GSE73189"

plot_df$group <- ifelse(grepl("^leave_out_", plot_df$analysis), "Leave-one-study-out", "Primary sensitivity")
plot_df$analysis_label <- factor(plot_df$analysis_label, levels = rev(plot_df$analysis_label[order(plot_df$pearson_logFC)]))

p1 <- ggplot(plot_df, aes(x = pearson_logFC, y = analysis_label, fill = group)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", pearson_logFC)), hjust = -0.15, size = 3.2) +
  scale_x_continuous(limits = c(0, 1.08), expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = palette_publication) +
  labs(x = "Pearson correlation of disease logFC vs baseline", y = NULL, title = "Effect-size stability") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

p2 <- ggplot(plot_df, aes(x = shared_sig, y = analysis_label, fill = group)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = shared_sig), hjust = -0.1, size = 3.2) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  scale_fill_manual(values = palette_publication) +
  labs(x = "Shared significant genes vs baseline\n(FDR < 0.05 and |logFC| >= 1)", y = NULL, title = "Signature overlap") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

pdf("ipf_sensitivity_results/figures_sensitivity_summary.pdf", width = 12, height = 6.5, useDingbats = FALSE)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(p1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
invisible(dev.off())

write.csv(plot_df, "ipf_sensitivity_results/tables/sensitivity_plot_data.csv", row.names = FALSE)
