#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(igraph)
  library(tidygraph)
  library(ggraph)
  library(ggplot2)
  library(scales)
})

palette_publication <- list(
  blue = "#0072B2",
  orange = "#D55E00",
  gray = "#4D4D4D",
  cream = "#F7F4EC",
  navy = "#16324F"
)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grep(file_arg, args)])
script_dir <- if (length(script_path) == 0) getwd() else dirname(normalizePath(script_path[1]))
setwd(script_dir)

edges <- read.csv("ipf_analysis_results/network/differential_coexpression_edges.csv", stringsAsFactors = FALSE)
nodes <- read.csv("ipf_analysis_results/network/differential_coexpression_nodes.csv", stringsAsFactors = FALSE)
de_tbl <- read.csv("ipf_analysis_results/tables/de_disease_average.csv", stringsAsFactors = FALSE)

edge_keep <- edges[order(-abs(edges$delta_cor), -edges$max_abs_cor), , drop = FALSE]
edge_keep <- head(edge_keep, 24)
edge_keep$edge_group <- ifelse(edge_keep$delta_cor > 0, "Stronger in IPF", "Stronger in Healthy")
edge_keep$weight_abs <- abs(edge_keep$delta_cor)

node_ids <- unique(c(edge_keep$gene1, edge_keep$gene2))
node_df <- nodes[nodes$gene_id %in% node_ids, , drop = FALSE]
node_df <- node_df[match(node_ids, node_df$gene_id), , drop = FALSE]
node_df$label_group <- ifelse(node_df$degree >= 3, "label", "nolabel")
node_df$display_label <- ifelse(node_df$label_group == "label", node_df$gene_label, "")
node_df$logFC <- de_tbl$logFC[match(node_df$gene_id, de_tbl$gene_id)]

graph_tbl <- tbl_graph(
  nodes = node_df[, c("gene_id", "gene_label", "degree", "label_group", "display_label", "logFC")],
  edges = edge_keep[, c("gene1", "gene2", "delta_cor", "weight_abs", "edge_group")],
  directed = FALSE
)

set.seed(42)
p <- ggraph(graph_tbl, layout = "fr") +
  geom_edge_link(aes(width = weight_abs, color = edge_group), alpha = 0.8, show.legend = TRUE) +
  scale_edge_width(range = c(0.8, 2.8), guide = "none") +
  scale_edge_color_manual(values = c("Stronger in IPF" = palette_publication$orange, "Stronger in Healthy" = palette_publication$blue)) +
  geom_node_point(aes(size = degree), color = palette_publication$navy, fill = palette_publication$cream, shape = 21, stroke = 1.1) +
  scale_size_continuous(range = c(4, 11)) +
  geom_node_text(aes(label = display_label), repel = TRUE, size = 4, color = palette_publication$gray) +
  labs(
    title = "Exploratory Differential Coexpression Network",
    subtitle = "Top differential-correlation edges among disease-associated genes",
    edge_color = NULL
  ) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11)
  )

ggsave(
  filename = "ipf_analysis_results/figures/network_differential_coexpression_graph.pdf",
  plot = p,
  width = 10.5,
  height = 8.5,
  device = cairo_pdf
)

p_logfc <- ggraph(graph_tbl, layout = "fr") +
  geom_edge_link(aes(width = weight_abs, color = edge_group), alpha = 0.75, show.legend = TRUE) +
  scale_edge_width(range = c(0.8, 2.8), guide = "none") +
  scale_edge_color_manual(values = c("Stronger in IPF" = palette_publication$orange, "Stronger in Healthy" = palette_publication$blue)) +
  geom_node_point(aes(size = degree, fill = logFC), color = palette_publication$navy, shape = 21, stroke = 1.0) +
  scale_size_continuous(range = c(4, 11)) +
  scale_fill_gradient2(
    low = muted(palette_publication$blue),
    mid = "white",
    high = muted(palette_publication$orange),
    midpoint = 0,
    name = "Disease logFC"
  ) +
  geom_node_text(aes(label = display_label), repel = TRUE, size = 4, color = palette_publication$gray) +
  guides(size = "none") +
  labs(
    title = "Exploratory Differential Coexpression Network",
    subtitle = "Top differential-correlation edges with nodes colored by disease logFC",
    edge_color = NULL
  ) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11)
  )

ggsave(
  filename = "ipf_analysis_results/figures/network_differential_coexpression_graph_logfc.pdf",
  plot = p_logfc,
  width = 10.5,
  height = 8.5,
  device = cairo_pdf
)
