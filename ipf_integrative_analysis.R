#!/usr/bin/env Rscript
# Perform integrative analysis of combined IPF RNA-seq dataset, including:
# - Data loading and preprocessing
# - Differential expression analysis with voom/limma and duplicateCorrelation
# - Visualization of results with volcano plots, heatmaps, and PCA
# - Gene set enrichment analysis with fgsea and msigdbr
# - Differential co-expression analysis for top disease-associated genes


suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(AnnotationDbi)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(fgsea)
  library(msigdbr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

palette_publication <- list(
  blue = "#0072B2",
  orange = "#D55E00",
  teal = "#009E73",
  gold = "#E69F00",
  gray = "#7F7F7F",
  light_gray = "#C7C7C7"
)

to_pdf_path <- function(path) {
  paste0(sub("\\.[^.]+$", "", path), ".pdf")
}

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", args[grep(file_arg, args)])
  if (length(script_path) == 0) {
    return(getwd())
  }
  dirname(normalizePath(script_path[1]))
}

script_dir <- get_script_dir()
setwd(script_dir)

results_dir <- file.path(script_dir, "ipf_analysis_results")
tables_dir <- file.path(results_dir, "tables")
figures_dir <- file.path(results_dir, "figures")
gsea_dir <- file.path(results_dir, "gsea")
network_dir <- file.path(results_dir, "network")
qc_dir <- file.path(results_dir, "qc")

for (path in c(results_dir, tables_dir, figures_dir, gsea_dir, network_dir, qc_dir)) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
}

metadata <- read.csv("combined_metadata.csv", stringsAsFactors = FALSE, check.names = FALSE)
counts <- read.csv("combined_study_matrix.csv", row.names = 1, check.names = FALSE)
annotation <- read.csv("ensembl_to_gene.csv", stringsAsFactors = FALSE)

required_cols <- c("Run", "Condition", "Sex", "Age", "GEO", "LibraryLayout", "Subject")
missing_required <- setdiff(required_cols, names(metadata))
if (length(missing_required) > 0) {
  stop("Missing required metadata columns: ", paste(missing_required, collapse = ", "))
}

metadata <- metadata[metadata$Run %in% colnames(counts), , drop = FALSE]
metadata <- metadata[complete.cases(metadata[, required_cols]), , drop = FALSE]
metadata <- metadata[metadata$Condition %in% c("Healthy", "IPF"), , drop = FALSE]
metadata <- metadata[metadata$Sex %in% c("female", "male"), , drop = FALSE]
metadata <- metadata[order(metadata$Run), , drop = FALSE]
counts <- counts[, metadata$Run, drop = FALSE]
metadata <- metadata[match(colnames(counts), metadata$Run), , drop = FALSE]

stopifnot(identical(metadata$Run, colnames(counts)))

metadata$Condition <- factor(metadata$Condition, levels = c("Healthy", "IPF"))
metadata$Sex <- factor(metadata$Sex, levels = c("female", "male"))
metadata$GEO <- factor(metadata$GEO)
metadata$LibraryLayout <- factor(metadata$LibraryLayout)
metadata$Subject <- factor(metadata$Subject)
metadata$Group <- factor(
  paste(metadata$Condition, metadata$Sex, sep = "_"),
  levels = c("Healthy_female", "Healthy_male", "IPF_female", "IPF_male")
)
metadata$Age_z <- as.numeric(scale(metadata$Age))

annotation$gene_id <- annotation$ensembl_id
annotation$gene_id_noversion <- annotation$ensembl_id_no_dot
missing_symbol <- is.na(annotation$gene_name) | annotation$gene_name == ""
annotation$fallback_symbol <- NA_character_
ens_lookup <- unique(annotation$gene_id_noversion[missing_symbol & !is.na(annotation$gene_id_noversion)])
if (length(ens_lookup) > 0) {
  ens_map <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = ens_lookup,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
  annotation$fallback_symbol[missing_symbol] <- unname(ens_map[annotation$gene_id_noversion[missing_symbol]])
}
annotation$gene_label <- ifelse(
  !is.na(annotation$gene_name) & annotation$gene_name != "",
  annotation$gene_name,
  ifelse(
    !is.na(annotation$fallback_symbol) & annotation$fallback_symbol != "",
    annotation$fallback_symbol,
    annotation$gene_id_noversion
  )
)

gene_annot <- data.frame(
  gene_id = rownames(counts),
  gene_id_noversion = sub("\\..*$", "", rownames(counts)),
  stringsAsFactors = FALSE
)
gene_annot <- merge(gene_annot, annotation[, c("gene_id", "gene_id_noversion", "gene_name", "gene_label", "entrezid")],
                    by = c("gene_id", "gene_id_noversion"), all.x = TRUE, sort = FALSE)
gene_annot$gene_label[is.na(gene_annot$gene_label)] <- gene_annot$gene_id_noversion[is.na(gene_annot$gene_label)]
gene_annot <- gene_annot[match(rownames(counts), gene_annot$gene_id), , drop = FALSE]

message("Building DGEList and filtering low-expression genes...")
dge <- DGEList(counts = as.matrix(counts))
keep <- filterByExpr(dge, group = metadata$Group)
dge <- dge[keep, , keep.lib.sizes = FALSE]
gene_annot <- gene_annot[keep, , drop = FALSE]
dge <- calcNormFactors(dge)

design <- model.matrix(~ 0 + Group + GEO + LibraryLayout + Group:Age_z, data = metadata)
colnames(design) <- make.names(colnames(design))

write.csv(metadata, file.path(qc_dir, "analysis_metadata_used.csv"), row.names = FALSE)
write.csv(
  data.frame(
    metric = c("samples_total", "genes_input", "genes_after_filter", "subjects", "repeated_subjects"),
    value = c(nrow(metadata), nrow(counts), nrow(dge), nlevels(metadata$Subject), sum(table(metadata$Subject) > 1))
  ),
  file.path(qc_dir, "analysis_summary.csv"),
  row.names = FALSE
)

message("Running voom/limma with duplicateCorrelation...")
v0 <- voom(dge, design, plot = FALSE)
var_rank <- order(apply(v0$E, 1, var), decreasing = TRUE)
cor_subset <- seq_len(min(500, length(var_rank)))
dupcor <- duplicateCorrelation(v0[var_rank[cor_subset], ], design = design, block = metadata$Subject)
v <- voom(dge, design, plot = FALSE, block = metadata$Subject, correlation = dupcor$consensus)
fit <- lmFit(v, design = design, block = metadata$Subject, correlation = dupcor$consensus)

contrast_matrix <- makeContrasts(
  disease_female = GroupIPF_female - GroupHealthy_female,
  disease_male = GroupIPF_male - GroupHealthy_male,
  disease_average = ((GroupIPF_female - GroupHealthy_female) + (GroupIPF_male - GroupHealthy_male)) / 2,
  sex_healthy = GroupHealthy_male - GroupHealthy_female,
  sex_ipf = GroupIPF_male - GroupIPF_female,
  disease_sex_interaction = (GroupIPF_male - GroupHealthy_male) - (GroupIPF_female - GroupHealthy_female),
  age_disease_female = GroupIPF_female.Age_z - GroupHealthy_female.Age_z,
  age_disease_male = GroupIPF_male.Age_z - GroupHealthy_male.Age_z,
  age_disease_average = ((GroupIPF_female.Age_z - GroupHealthy_female.Age_z) +
    (GroupIPF_male.Age_z - GroupHealthy_male.Age_z)) / 2,
  age_sex_healthy = GroupHealthy_male.Age_z - GroupHealthy_female.Age_z,
  age_sex_ipf = GroupIPF_male.Age_z - GroupIPF_female.Age_z,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2, robust = TRUE)

write.csv(
  data.frame(consensus_correlation = dupcor$consensus),
  file.path(qc_dir, "duplicate_correlation.csv"),
  row.names = FALSE
)

volcano_plot <- function(tbl, title, outfile, lfc_cutoff = 1, fdr_cutoff = 0.05, top_n_each = 10) {
  tbl$plot_fdr <- pmax(tbl$adj.P.Val, .Machine$double.xmin)
  tbl$significant <- with(tbl, adj.P.Val < fdr_cutoff & abs(logFC) >= lfc_cutoff)
  tbl$label <- ""
  tbl$label_candidate <- !is.na(tbl$gene_label) & tbl$gene_label != "" & !grepl("^ENSG", tbl$gene_label)

  up_hits <- tbl[tbl$significant & tbl$label_candidate & tbl$logFC > 0, , drop = FALSE]
  up_hits <- up_hits[order(up_hits$adj.P.Val, -abs(up_hits$logFC)), , drop = FALSE]
  up_hits <- up_hits[!duplicated(up_hits$gene_label), , drop = FALSE]
  up_hits <- head(up_hits, top_n_each)

  down_hits <- tbl[tbl$significant & tbl$label_candidate & tbl$logFC < 0, , drop = FALSE]
  down_hits <- down_hits[order(down_hits$adj.P.Val, -abs(down_hits$logFC)), , drop = FALSE]
  down_hits <- down_hits[!duplicated(down_hits$gene_label), , drop = FALSE]
  down_hits <- head(down_hits, top_n_each)

  ranked_hits <- rbind(up_hits, down_hits)
  if (nrow(ranked_hits) > 0) {
    tbl$label[match(ranked_hits$gene_id, tbl$gene_id)] <- ranked_hits$gene_label
  }

  p <- ggplot(tbl, aes(x = logFC, y = -log10(plot_fdr), color = significant)) +
    geom_point(alpha = 0.7, size = 1.2) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = 2, color = "gray50") +
    geom_hline(yintercept = -log10(fdr_cutoff), linetype = 2, color = "gray50") +
    scale_color_manual(values = c(`TRUE` = palette_publication$orange, `FALSE` = palette_publication$light_gray)) +
    geom_text_repel(aes(label = label), size = 3, max.overlaps = Inf, box.padding = 0.35, min.segment.length = 0) +
    labs(title = title, x = "log2 fold change", y = "-log10 FDR") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")

  ggsave(to_pdf_path(outfile), p, width = 8, height = 6, device = cairo_pdf)
}

make_heatmap <- function(expr_mat, genes, annotation_df, sample_meta, outfile, title) {
  genes <- unique(genes)
  genes <- genes[genes %in% rownames(expr_mat)]
  if (length(genes) < 2) {
    return(invisible(NULL))
  }
  mat <- expr_mat[genes, , drop = FALSE]
  rownames(mat) <- make.unique(annotation_df$gene_label[match(rownames(mat), annotation_df$gene_id)])
  mat <- t(scale(t(mat)))
  ann_col <- sample_meta[, c("Condition", "Sex", "GEO", "LibraryLayout"), drop = FALSE]
  rownames(ann_col) <- sample_meta$Run
  pheatmap(
    mat,
    annotation_col = ann_col,
    show_colnames = FALSE,
    fontsize_row = 8,
    main = title,
    filename = to_pdf_path(outfile),
    width = 11,
    height = 10
  )
}

plot_pca <- function(expr_mat, sample_meta, outfile, title) {
  vars <- apply(expr_mat, 1, var)
  keep_genes <- names(sort(vars, decreasing = TRUE))[seq_len(min(1000, length(vars)))]
  pca <- prcomp(t(expr_mat[keep_genes, , drop = FALSE]), center = TRUE, scale. = FALSE)
  pca_df <- data.frame(
    sample_meta,
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2]
  )
  percent <- summary(pca)$importance[2, 1:2] * 100
  p <- ggplot(pca_df, aes(PC1, PC2, color = Condition, shape = Sex)) +
    geom_point(size = 3, alpha = 0.9) +
    facet_wrap(~ GEO) +
    scale_color_manual(values = c("Healthy" = palette_publication$blue, "IPF" = palette_publication$orange)) +
    labs(
      title = title,
      x = sprintf("PC1 (%.1f%%)", percent[1]),
      y = sprintf("PC2 (%.1f%%)", percent[2])
    ) +
    theme_bw(base_size = 12)
  ggsave(to_pdf_path(outfile), p, width = 12, height = 8, device = cairo_pdf)
}

clean_result_table <- function(fit_obj, coef_name, annotation_df) {
  tbl <- topTable(fit_obj, coef = coef_name, number = Inf, sort.by = "P")
  tbl$gene_id <- rownames(tbl)
  ann <- annotation_df[match(tbl$gene_id, annotation_df$gene_id), c("gene_id_noversion", "gene_name", "gene_label", "entrezid")]
  tbl <- cbind(tbl, ann)
  tbl <- tbl[order(tbl$P.Value, -abs(tbl$logFC)), , drop = FALSE]
  tbl$direction <- ifelse(tbl$logFC > 0, "up", "down")
  tbl
}

subset_de_hits <- function(tbl, fdr_cutoff = 0.05, lfc_cutoff = 1) {
  tbl[tbl$adj.P.Val < fdr_cutoff & abs(tbl$logFC) >= lfc_cutoff, , drop = FALSE]
}

align_result_table <- function(tbl, gene_ids, table_name) {
  aligned_tbl <- tbl[match(gene_ids, tbl$gene_id), , drop = FALSE]
  if (anyNA(aligned_tbl$gene_id)) {
    stop("Could not align result table by gene_id for ", table_name, ".")
  }
  aligned_tbl
}

logcpm <- cpm(dge, log = TRUE, prior.count = 2)
plot_pca(logcpm, metadata, file.path(qc_dir, "pca_logcpm_raw.png"), "PCA on filtered logCPM")

viz_design <- model.matrix(~ 0 + Group, metadata)
batch_corrected_expr <- removeBatchEffect(
  v$E,
  batch = metadata$GEO,
  batch2 = metadata$LibraryLayout,
  design = viz_design,
  covariates = cbind(Age_z = metadata$Age_z)
)
plot_pca(batch_corrected_expr, metadata, file.path(qc_dir, "pca_batch_corrected.png"), "PCA after batch adjustment")

message("Exporting DE results and plots...")
contrast_names <- colnames(contrast_matrix)
summary_rows <- list()
result_tables <- list()

for (coef_name in contrast_names) {
  tbl <- clean_result_table(fit2, coef_name, gene_annot)
  de_hits <- subset_de_hits(tbl)
  top_hit <- if (nrow(de_hits) > 0) de_hits[1, , drop = FALSE] else tbl[1, , drop = FALSE]
  result_tables[[coef_name]] <- tbl
  write.csv(tbl, file.path(tables_dir, paste0("de_", coef_name, ".csv")), row.names = FALSE)
  summary_rows[[coef_name]] <- data.frame(
    contrast = coef_name,
    genes_tested = nrow(tbl),
    fdr_005 = sum(tbl$adj.P.Val < 0.05),
    fdr_005_lfc1 = sum(tbl$adj.P.Val < 0.05 & abs(tbl$logFC) >= 1),
    top_gene = top_hit$gene_label[1],
    top_logFC = top_hit$logFC[1],
    top_fdr = top_hit$adj.P.Val[1]
  )
}

summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(tables_dir, "de_summary.csv"), row.names = FALSE)

for (coef_name in c("disease_average", "disease_female", "disease_male", "disease_sex_interaction",
                    "age_disease_average", "age_disease_female", "age_disease_male")) {
  volcano_plot(
    result_tables[[coef_name]],
    gsub("_", " ", coef_name),
    file.path(figures_dir, paste0("volcano_", coef_name, ".png"))
  )
}

top_heatmap_genes <- subset_de_hits(result_tables$disease_average)
top_heatmap_genes <- head(top_heatmap_genes$gene_id, 50)
make_heatmap(
  batch_corrected_expr,
  top_heatmap_genes,
  gene_annot,
  metadata,
  file.path(figures_dir, "heatmap_disease_average_top50.png"),
  "Top disease-associated genes"
)

female_tbl <- result_tables$disease_female
male_tbl <- align_result_table(result_tables$disease_male, female_tbl$gene_id, "disease_male")
stopifnot(identical(female_tbl$gene_id, male_tbl$gene_id))

venn_input <- cbind(
  disease_female = female_tbl$adj.P.Val < 0.05 & abs(female_tbl$logFC) >= 1,
  disease_male = male_tbl$adj.P.Val < 0.05 & abs(male_tbl$logFC) >= 1
)
rownames(venn_input) <- female_tbl$gene_id
pdf(file.path(figures_dir, "venn_disease_female_male.pdf"), width = 7.8, height = 6.7, useDingbats = FALSE)
vennDiagram(venn_input, circle.col = c(palette_publication$orange, palette_publication$blue))
invisible(dev.off())

overlap_tbl <- data.frame(
  gene_id = female_tbl$gene_id,
  gene_label = female_tbl$gene_label,
  female_sig = venn_input[, "disease_female"],
  male_sig = venn_input[, "disease_male"],
  female_logFC = female_tbl$logFC,
  male_logFC = male_tbl$logFC
)
write.csv(overlap_tbl, file.path(tables_dir, "disease_female_male_overlap.csv"), row.names = FALSE)

message("Running GSEA...")
msig_all <- tryCatch(msigdbr(species = "Homo sapiens"), error = function(e) msigdbr())
collections <- list(
  Reactome = list(collection = "C2", subcollection = "CP:REACTOME"),
  KEGG = list(collection = "C2", subcollection = c("CP:KEGG", "CP:KEGG_LEGACY")),
  GO_BP = list(collection = "C5", subcollection = "GO:BP"),
  GO_MF = list(collection = "C5", subcollection = "GO:MF"),
  GO_CC = list(collection = "C5", subcollection = "GO:CC")
)

build_stats <- function(tbl, stat_col = "t") {
  stats <- tbl[, c("gene_label", stat_col), drop = FALSE]
  stats <- stats[!is.na(stats$gene_label) & stats$gene_label != "", , drop = FALSE]
  stats <- stats[order(abs(stats[[stat_col]]), decreasing = TRUE), , drop = FALSE]
  stats <- stats[!duplicated(stats$gene_label), , drop = FALSE]
  vals <- stats[[stat_col]]
  names(vals) <- stats$gene_label
  vals <- sort(vals, decreasing = TRUE)
  vals
}

plot_gsea <- function(gsea_tbl, title, outfile) {
  sig <- gsea_tbl[gsea_tbl$padj < 0.05, , drop = FALSE]
  if (nrow(sig) == 0) {
    return(invisible(NULL))
  }
  sig <- sig[order(abs(sig$NES), decreasing = TRUE), , drop = FALSE]
  sig <- head(sig, 12)
  sig$pathway <- factor(sig$pathway, levels = rev(sig$pathway))
  sig <- sig[match(levels(sig$pathway), sig$pathway), , drop = FALSE]
  p <- ggplot(sig, aes(x = pathway, y = NES, fill = NES > 0)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c(`TRUE` = palette_publication$orange, `FALSE` = palette_publication$blue)) +
    labs(title = title, x = NULL, y = "NES") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none")
  ggsave(to_pdf_path(outfile), p, width = 9, height = 6, device = cairo_pdf)
}

gsea_targets <- c("disease_average", "disease_female", "disease_male", "disease_sex_interaction", "age_disease_average")

for (contrast_name in gsea_targets) {
  stats_vec <- build_stats(result_tables[[contrast_name]], stat_col = "t")
  for (collection_name in names(collections)) {
    collection_spec <- collections[[collection_name]]
    selected_subcollection <- collection_spec$subcollection
    if (collection_name == "KEGG") {
      available_kegg <- intersect(collection_spec$subcollection, unique(msig_all$gs_subcollection))
      preferred_order <- c("CP:KEGG", "CP:KEGG_LEGACY")
      selected_subcollection <- preferred_order[preferred_order %in% available_kegg]
      if (length(selected_subcollection) == 0) {
        next
      }
      selected_subcollection <- selected_subcollection[1]
    }
    gs_df <- msig_all[
      msig_all$gs_collection == collection_spec$collection &
        msig_all$gs_subcollection %in% selected_subcollection,
      c("gs_name", "gene_symbol")
    ]
    gs_df <- unique(gs_df[!is.na(gs_df$gene_symbol) & gs_df$gene_symbol != "", , drop = FALSE])
    pathways <- split(gs_df$gene_symbol, gs_df$gs_name)
    pathways <- pathways[lengths(pathways) >= 10 & lengths(pathways) <= 500]
    if (length(pathways) == 0) {
      next
    }
    fgsea_res <- fgsea(pathways = pathways, stats = stats_vec, minSize = 10, maxSize = 500)
    fgsea_res <- fgsea_res[order(fgsea_res$padj, -abs(fgsea_res$NES)), ]
    if ("leadingEdge" %in% names(fgsea_res)) {
      fgsea_res$leadingEdge <- vapply(fgsea_res$leadingEdge, function(x) paste(x, collapse = ";"), character(1))
    }
    out_stub <- paste0("gsea_", contrast_name, "_", collection_name)
    write.csv(as.data.frame(fgsea_res), file.path(gsea_dir, paste0(out_stub, ".csv")), row.names = FALSE)
    plot_gsea(as.data.frame(fgsea_res), paste(contrast_name, collection_name), file.path(figures_dir, paste0(out_stub, ".png")))
  }
}

fibrosis_sets <- msig_all[grepl("fibrosis|pulmonary fibrosis|interstitial lung", msig_all$gs_name, ignore.case = TRUE), ]
if (nrow(fibrosis_sets) > 0) {
  fibrosis_sets <- unique(fibrosis_sets[, c("gs_name", "gene_symbol")])
  write.csv(fibrosis_sets, file.path(gsea_dir, "candidate_fibrosis_gene_sets.csv"), row.names = FALSE)
}

message("Building differential co-expression summary...")
coexp_expr <- removeBatchEffect(
  v$E,
  batch = metadata$GEO,
  batch2 = metadata$LibraryLayout,
  design = model.matrix(~ 0 + Condition, metadata),
  covariates = cbind(Age_z = metadata$Age_z, SexMale = ifelse(metadata$Sex == "male", 1, 0))
)

network_unit <- interaction(metadata$Subject, metadata$Condition, metadata$Sex, drop = TRUE)
coexp_expr_subject_sum <- rowsum(t(coexp_expr), group = network_unit, reorder = FALSE)
unit_levels <- rownames(coexp_expr_subject_sum)
unit_counts <- as.integer(table(factor(as.character(network_unit), levels = unit_levels)))
if (length(unit_counts) != nrow(coexp_expr_subject_sum)) {
  stop("Network-unit counts do not align with collapsed expression matrix.")
}
coexp_expr_subject <- t(coexp_expr_subject_sum / unit_counts)
subject_meta <- metadata[match(unit_levels, as.character(network_unit)), , drop = FALSE]
if (anyNA(subject_meta$Run)) {
  stop("Network-unit metadata do not align with collapsed expression matrix.")
}
subject_meta$network_unit <- unit_levels
subject_meta <- subject_meta[match(colnames(coexp_expr_subject), subject_meta$network_unit), , drop = FALSE]

network_genes <- result_tables$disease_average
network_genes <- network_genes[network_genes$adj.P.Val < 0.05 & abs(network_genes$logFC) >= 1, , drop = FALSE]
if (nrow(network_genes) < 30) {
  network_genes <- head(result_tables$disease_average[order(result_tables$disease_average$adj.P.Val), , drop = FALSE], 100)
} else {
  network_genes <- head(network_genes, 150)
}
network_gene_ids <- network_genes$gene_id[network_genes$gene_id %in% rownames(coexp_expr_subject)]

healthy_idx <- subject_meta$Condition == "Healthy"
ipf_idx <- subject_meta$Condition == "IPF"
if (length(network_gene_ids) >= 10 && sum(healthy_idx) >= 10 && sum(ipf_idx) >= 10) {
  cor_healthy <- cor(t(coexp_expr_subject[network_gene_ids, healthy_idx, drop = FALSE]), method = "spearman")
  cor_ipf <- cor(t(coexp_expr_subject[network_gene_ids, ipf_idx, drop = FALSE]), method = "spearman")
  delta_cor <- cor_ipf - cor_healthy

  upper_idx <- which(upper.tri(delta_cor), arr.ind = TRUE)
  edge_tbl <- data.frame(
    gene1 = network_gene_ids[upper_idx[, 1]],
    gene2 = network_gene_ids[upper_idx[, 2]],
    cor_healthy = cor_healthy[upper_idx],
    cor_ipf = cor_ipf[upper_idx],
    delta_cor = delta_cor[upper_idx]
  )
  edge_tbl$max_abs_cor <- pmax(abs(edge_tbl$cor_healthy), abs(edge_tbl$cor_ipf))
  edge_tbl <- edge_tbl[abs(edge_tbl$delta_cor) >= 0.5 & edge_tbl$max_abs_cor >= 0.6, , drop = FALSE]
  edge_tbl <- edge_tbl[order(-abs(edge_tbl$delta_cor), -edge_tbl$max_abs_cor), , drop = FALSE]
  edge_tbl$gene1_label <- gene_annot$gene_label[match(edge_tbl$gene1, gene_annot$gene_id)]
  edge_tbl$gene2_label <- gene_annot$gene_label[match(edge_tbl$gene2, gene_annot$gene_id)]
  write.csv(edge_tbl, file.path(network_dir, "differential_coexpression_edges.csv"), row.names = FALSE)

  if (nrow(edge_tbl) > 0) {
    degree_tab <- sort(table(c(edge_tbl$gene1, edge_tbl$gene2)), decreasing = TRUE)
    node_tbl <- data.frame(
      gene_id = names(degree_tab),
      degree = as.integer(degree_tab),
      gene_label = gene_annot$gene_label[match(names(degree_tab), gene_annot$gene_id)]
    )
    write.csv(node_tbl, file.path(network_dir, "differential_coexpression_nodes.csv"), row.names = FALSE)

    top_nodes <- head(node_tbl$gene_id, 40)
    heatmap_mat <- delta_cor[top_nodes, top_nodes, drop = FALSE]
    rownames(heatmap_mat) <- make.unique(gene_annot$gene_label[match(rownames(heatmap_mat), gene_annot$gene_id)])
    colnames(heatmap_mat) <- make.unique(gene_annot$gene_label[match(colnames(heatmap_mat), gene_annot$gene_id)])
    pheatmap(
      heatmap_mat,
      color = colorRampPalette(c(palette_publication$blue, "white", palette_publication$orange))(100),
      main = "Differential correlation (IPF - Healthy)",
      filename = file.path(figures_dir, "heatmap_differential_coexpression_top_hubs.pdf"),
      width = 10,
      height = 9
    )
  }
}