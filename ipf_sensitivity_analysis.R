#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

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

results_dir <- file.path(script_dir, "ipf_sensitivity_results")
tables_dir <- file.path(results_dir, "tables")
qc_dir <- file.path(results_dir, "qc")
for (path in c(results_dir, tables_dir, qc_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

metadata_all <- read.csv("combined_metadata.csv", stringsAsFactors = FALSE, check.names = FALSE)
counts_all <- read.csv("combined_study_matrix.csv", row.names = 1, check.names = FALSE)
annotation <- read.csv("ensembl_to_gene.csv", stringsAsFactors = FALSE)
annotation$gene_id <- annotation$ensembl_id
annotation$gene_label <- ifelse(
  is.na(annotation$gene_name) | annotation$gene_name == "",
  annotation$ensembl_id_no_dot,
  annotation$gene_name
)

baseline_tables <- list(
  disease_average = read.csv(file.path("ipf_analysis_results", "tables", "de_disease_average.csv"), stringsAsFactors = FALSE),
  disease_female = read.csv(file.path("ipf_analysis_results", "tables", "de_disease_female.csv"), stringsAsFactors = FALSE),
  disease_male = read.csv(file.path("ipf_analysis_results", "tables", "de_disease_male.csv"), stringsAsFactors = FALSE)
)
baseline_gene_ids <- baseline_tables$disease_average$gene_id

missing_symbol <- is.na(annotation$gene_name) | annotation$gene_name == ""
annotation$fallback_symbol <- NA_character_
ens_lookup <- unique(annotation$ensembl_id_no_dot[missing_symbol & !is.na(annotation$ensembl_id_no_dot)])
if (length(ens_lookup) > 0) {
  ens_map <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = ens_lookup,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
  annotation$fallback_symbol[missing_symbol] <- unname(ens_map[annotation$ensembl_id_no_dot[missing_symbol]])
}
annotation$gene_label <- ifelse(
  !is.na(annotation$gene_name) & annotation$gene_name != "",
  annotation$gene_name,
  ifelse(
    !is.na(annotation$fallback_symbol) & annotation$fallback_symbol != "",
    annotation$fallback_symbol,
    annotation$ensembl_id_no_dot
  )
)

prepare_data <- function(include_age = TRUE, exclude_geo = NULL) {
  required_cols <- c("Run", "Condition", "Sex", "GEO", "LibraryLayout", "Subject")
  if (include_age) {
    required_cols <- c(required_cols, "Age")
  }

  meta <- metadata_all[metadata_all$Run %in% colnames(counts_all), , drop = FALSE]
  meta <- meta[complete.cases(meta[, required_cols]), , drop = FALSE]
  meta <- meta[meta$Condition %in% c("Healthy", "IPF"), , drop = FALSE]
  meta <- meta[meta$Sex %in% c("female", "male"), , drop = FALSE]
  if (!is.null(exclude_geo)) {
    meta <- meta[meta$GEO != exclude_geo, , drop = FALSE]
  }
  meta <- meta[order(meta$Run), , drop = FALSE]
  counts <- counts_all[, meta$Run, drop = FALSE]
  meta <- meta[match(colnames(counts), meta$Run), , drop = FALSE]

  meta$Condition <- factor(meta$Condition, levels = c("Healthy", "IPF"))
  meta$Sex <- factor(meta$Sex, levels = c("female", "male"))
  meta$GEO <- factor(meta$GEO)
  meta$LibraryLayout <- factor(meta$LibraryLayout)
  meta$Subject <- factor(meta$Subject)
  meta$Group <- factor(
    paste(meta$Condition, meta$Sex, sep = "_"),
    levels = c("Healthy_female", "Healthy_male", "IPF_female", "IPF_male")
  )
  if (include_age) {
    meta$Age_z <- as.numeric(scale(meta$Age))
  }

  list(metadata = meta, counts = counts)
}

fit_sensitivity <- function(label, include_age = TRUE, exclude_geo = NULL) {
  dat <- prepare_data(include_age = include_age, exclude_geo = exclude_geo)
  meta <- dat$metadata
  counts <- dat$counts

  counts <- counts[baseline_gene_ids[baseline_gene_ids %in% rownames(counts)], , drop = FALSE]
  dge <- DGEList(counts = as.matrix(counts))
  dge <- calcNormFactors(dge)

  if (include_age) {
    design <- model.matrix(~ 0 + Group + GEO + LibraryLayout + Group:Age_z, data = meta)
  } else {
    design <- model.matrix(~ 0 + Group + GEO + LibraryLayout, data = meta)
  }
  colnames(design) <- make.names(colnames(design))

  v0 <- voom(dge, design, plot = FALSE)
  if (all(table(meta$Subject) == 1)) {
    dupcor <- list(consensus = 0)
  } else {
    var_rank <- order(apply(v0$E, 1, var), decreasing = TRUE)
    cor_subset <- seq_len(min(500, length(var_rank)))
    dupcor <- duplicateCorrelation(v0[var_rank[cor_subset], ], design = design, block = meta$Subject)
  }
  v <- voom(dge, design, plot = FALSE, block = meta$Subject, correlation = dupcor$consensus)
  fit <- lmFit(v, design, block = meta$Subject, correlation = dupcor$consensus)

  contrast_matrix <- makeContrasts(
    disease_female = GroupIPF_female - GroupHealthy_female,
    disease_male = GroupIPF_male - GroupHealthy_male,
    disease_average = ((GroupIPF_female - GroupHealthy_female) + (GroupIPF_male - GroupHealthy_male)) / 2,
    levels = design
  )
  fit2 <- eBayes(contrasts.fit(fit, contrast_matrix), robust = TRUE)

  build_result_table <- function(coef_name) {
    out <- topTable(fit2, coef = coef_name, number = Inf, sort.by = "P")
    out$gene_id <- rownames(out)
    out$gene_label <- annotation$gene_label[match(out$gene_id, annotation$gene_id)]
    out$gene_label[is.na(out$gene_label)] <- sub("\\..*$", "", out$gene_id[is.na(out$gene_label)])
    out[order(out$P.Value, -abs(out$logFC)), , drop = FALSE]
  }

  result_tables <- list(
    disease_average = build_result_table("disease_average"),
    disease_female = build_result_table("disease_female"),
    disease_male = build_result_table("disease_male")
  )

  out <- result_tables$disease_average

  summary_row <- data.frame(
    analysis = label,
    include_age = include_age,
    excluded_geo = ifelse(is.null(exclude_geo), "", exclude_geo),
    samples = nrow(meta),
    healthy = sum(meta$Condition == "Healthy"),
    ipf = sum(meta$Condition == "IPF"),
    genes_tested = nrow(out),
    fdr_005 = sum(out$adj.P.Val < 0.05),
    fdr_005_lfc1 = sum(out$adj.P.Val < 0.05 & abs(out$logFC) >= 1),
    top_gene = out$gene_label[1],
    top_logFC = out$logFC[1],
    top_fdr = out$adj.P.Val[1],
    dupcor = dupcor$consensus,
    stringsAsFactors = FALSE
  )

  for (contrast_name in names(result_tables)) {
    write.csv(
      result_tables[[contrast_name]],
      file.path(tables_dir, paste0("de_", contrast_name, "_", label, ".csv")),
      row.names = FALSE
    )
  }
  write.csv(meta, file.path(qc_dir, paste0("metadata_", label, ".csv")), row.names = FALSE)

  list(summary = summary_row, tables = result_tables)
}

compare_to_baseline <- function(tbl, baseline_tbl, label, contrast_name) {
  merged <- merge(
    baseline_tbl[, c("gene_id", "logFC", "adj.P.Val")],
    tbl[, c("gene_id", "logFC", "adj.P.Val")],
    by = "gene_id",
    suffixes = c("_baseline", paste0("_", label)),
    all = FALSE
  )

  baseline_sig <- with(merged, adj.P.Val_baseline < 0.05 & abs(logFC_baseline) >= 1)
  current_sig <- with(merged, merged[[paste0("adj.P.Val_", label)]] < 0.05 & abs(merged[[paste0("logFC_", label)]]) >= 1)

  data.frame(
    analysis = label,
    contrast = contrast_name,
    genes_compared = nrow(merged),
    pearson_logFC = cor(merged$logFC_baseline, merged[[paste0("logFC_", label)]], method = "pearson"),
    spearman_logFC = cor(merged$logFC_baseline, merged[[paste0("logFC_", label)]], method = "spearman"),
    shared_sig = sum(baseline_sig & current_sig),
    baseline_sig = sum(baseline_sig),
    current_sig = sum(current_sig),
    baseline_only = sum(baseline_sig & !current_sig),
    current_only = sum(!baseline_sig & current_sig),
    stringsAsFactors = FALSE
  )
}

message("Running sensitivity analyses...")

res_exclude_confounded <- fit_sensitivity("exclude_GSE184316", include_age = TRUE, exclude_geo = "GSE184316")
res_no_age <- fit_sensitivity("no_age", include_age = FALSE, exclude_geo = NULL)

loo_geos <- sort(unique(prepare_data(include_age = TRUE)$metadata$GEO))
loo_results <- list()
loo_summaries <- list()
loo_compares <- list()
for (geo in loo_geos) {
  label <- paste0("leave_out_", geo)
  fit_res <- fit_sensitivity(label, include_age = TRUE, exclude_geo = geo)
  loo_results[[label]] <- fit_res
  loo_summaries[[label]] <- fit_res$summary
  loo_compares[[label]] <- do.call(rbind, lapply(names(fit_res$tables), function(contrast_name) {
    compare_to_baseline(fit_res$tables[[contrast_name]], baseline_tables[[contrast_name]], label, contrast_name)
  }))
}

sensitivity_summary <- do.call(rbind, c(
  list(res_exclude_confounded$summary, res_no_age$summary),
  loo_summaries
))
write.csv(sensitivity_summary, file.path(tables_dir, "sensitivity_summary.csv"), row.names = FALSE)

comparison_summary <- do.call(rbind, c(
  list(do.call(rbind, lapply(names(res_exclude_confounded$tables), function(contrast_name) {
    compare_to_baseline(
      res_exclude_confounded$tables[[contrast_name]],
      baseline_tables[[contrast_name]],
      "exclude_GSE184316",
      contrast_name
    )
  }))),
  list(do.call(rbind, lapply(names(res_no_age$tables), function(contrast_name) {
    compare_to_baseline(
      res_no_age$tables[[contrast_name]],
      baseline_tables[[contrast_name]],
      "no_age",
      contrast_name
    )
  }))),
  loo_compares
))
write.csv(comparison_summary, file.path(tables_dir, "sensitivity_vs_baseline.csv"), row.names = FALSE)

top_gene_overlap <- function(tbl, n = 100) {
  head(tbl$gene_id[order(tbl$adj.P.Val, -abs(tbl$logFC))], n)
}

overlap_rows <- list(
  do.call(rbind, lapply(names(res_exclude_confounded$tables), function(contrast_name) {
    data.frame(
      analysis = "exclude_GSE184316",
      contrast = contrast_name,
      top100_overlap = sum(
        top_gene_overlap(baseline_tables[[contrast_name]]) %in% top_gene_overlap(res_exclude_confounded$tables[[contrast_name]])
      ),
      stringsAsFactors = FALSE
    )
  })),
  do.call(rbind, lapply(names(res_no_age$tables), function(contrast_name) {
    data.frame(
      analysis = "no_age",
      contrast = contrast_name,
      top100_overlap = sum(
        top_gene_overlap(baseline_tables[[contrast_name]]) %in% top_gene_overlap(res_no_age$tables[[contrast_name]])
      ),
      stringsAsFactors = FALSE
    )
  }))
)
for (geo in loo_geos) {
  label <- paste0("leave_out_", geo)
  overlap_rows[[label]] <- do.call(rbind, lapply(names(loo_results[[label]]$tables), function(contrast_name) {
    data.frame(
      analysis = label,
      contrast = contrast_name,
      top100_overlap = sum(
        top_gene_overlap(baseline_tables[[contrast_name]]) %in% top_gene_overlap(loo_results[[label]]$tables[[contrast_name]])
      ),
      stringsAsFactors = FALSE
    )
  }))
}
write.csv(do.call(rbind, overlap_rows), file.path(tables_dir, "top100_overlap_vs_baseline.csv"), row.names = FALSE)

message("Sensitivity analyses completed.")
