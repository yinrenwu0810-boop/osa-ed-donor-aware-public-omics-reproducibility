args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: Rscript run_GSE206528_doublet_sensitivity_DE.R <pseudobulk_dir> <output_dir> <workspace_root>")
pseudobulk_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = FALSE)
workspace_root <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(file.path(workspace_root, "analysis", ".Rlib"), .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(msigdbr)
  library(fgsea)
})

set.seed(20260814)
count_table <- fread(file.path(pseudobulk_dir, "pseudobulk_counts.csv"))
sample_table <- fread(file.path(pseudobulk_dir, "pseudobulk_samples.csv"))
setnames(sample_table, c("pseudobulk_id", "sample_accession", "donor", "condition", "cell_type"))
cell_counts <- fread(file.path(pseudobulk_dir, "celltype_counts_by_donor.tsv"))
gene_symbol <- count_table[[1L]]
counts <- as.matrix(count_table[, -1L])
storage.mode(counts) <- "double"
rownames(counts) <- gene_symbol
sample_table <- sample_table[match(colnames(counts), pseudobulk_id)]
stopifnot(identical(sample_table$pseudobulk_id, colnames(counts)))

hallmark_db <- msigdbr(species = "Homo sapiens", collection = "H")
pathways <- split(hallmark_db$gene_symbol, hallmark_db$gs_name)
cell_types <- sort(unique(sample_table$cell_type))
summary_rows <- list()

run_comparison <- function(target_cell_type, case_condition, comparison_name, min_cells = 20L) {
  metadata <- sample_table[cell_type == target_cell_type & condition %in% c("normal", case_condition)]
  eligible <- cell_counts[cell_type == target_cell_type & condition %in% c("normal", case_condition) & cell_count >= min_cells]
  metadata <- metadata[sample_accession %in% eligible$sample_accession]
  expected <- if (case_condition == "organic_ED_nonDM") c(normal = 3L, organic_ED_nonDM = 3L) else c(normal = 3L, organic_ED_DM = 2L)
  observed <- table(metadata$condition)
  if (!all(names(expected) %in% names(observed)) || any(observed[names(expected)] < expected)) return(NULL)
  idx <- match(metadata$pseudobulk_id, colnames(counts))
  group <- factor(metadata$condition, levels = c("normal", case_condition))
  model <- model.matrix(~ group)
  y <- DGEList(counts = counts[, idx, drop = FALSE], group = group)
  keep <- filterByExpr(y, design = model, min.count = 10)
  y <- normLibSizes(y[keep, , keep.lib.sizes = FALSE])
  fit <- glmQLFit(y, model, robust = TRUE)
  test <- glmQLFTest(fit, coef = 2)
  tab <- topTags(test, n = Inf, sort.by = "none")$table
  tab$gene_symbol <- rownames(tab)
  tab <- tab[, c("gene_symbol", "logFC", "logCPM", "F", "PValue", "FDR")]
  tab <- tab[order(tab$FDR, tab$PValue), ]
  safe_type <- gsub("[^A-Za-z0-9]+", "_", target_cell_type)
  fwrite(tab, file.path(output_dir, paste0("differential_", safe_type, "_", comparison_name, ".tsv")), sep = "\t")
  ranks <- test$table$F
  ranks <- sqrt(pmax(ranks, 0)) * sign(test$table$logFC)
  names(ranks) <- rownames(test$table)
  ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)
  fg <- fgseaMultilevel(pathways, ranks, minSize = 15, maxSize = 500, eps = 0)
  fg[, abs_NES := abs(NES)]
  setorder(fg, padj, -abs_NES)
  fg[, leadingEdge := vapply(leadingEdge, paste, collapse = ";", character(1))]
  fg[, abs_NES := NULL]
  fwrite(fg, file.path(output_dir, paste0("GSEA_Hallmark_", safe_type, "_", comparison_name, ".tsv")), sep = "\t")
  data.table(cell_type = target_cell_type, comparison = comparison_name, donors = ncol(y), min_cells_per_donor = min(eligible$cell_count), tested_genes = nrow(tab), FDR_0.05 = sum(tab$FDR < 0.05), hallmark_FDR_0.05 = sum(fg$padj < 0.05, na.rm = TRUE))
}

for (cell_type_name in cell_types) {
  primary <- run_comparison(cell_type_name, "organic_ED_nonDM", "nonDM_ED_vs_normal")
  if (!is.null(primary)) summary_rows[[length(summary_rows) + 1L]] <- primary
  exploratory <- run_comparison(cell_type_name, "organic_ED_DM", "DMED_vs_normal")
  if (!is.null(exploratory)) summary_rows[[length(summary_rows) + 1L]] <- exploratory
}
summary <- rbindlist(summary_rows, fill = TRUE)
fwrite(summary, file.path(output_dir, "analysis_summary.tsv"), sep = "\t")

# Prespecified result-known candidate and Hallmark comparisons against the historical 03A outputs.
historical_dir <- file.path(workspace_root, "analysis", "results", "GSE206528", "pseudobulk_DE")
candidates <- c("TYMS", "EFNB2", "LRRC17")
candidate_rows <- list(); hallmark_rows <- list()
for (new_file in list.files(output_dir, pattern = "^(differential|GSEA_Hallmark)_.*\\.tsv$", full.names = FALSE)) {
  old_file <- file.path(historical_dir, new_file)
  if (!file.exists(old_file)) next
  if (startsWith(new_file, "differential_")) {
    new <- fread(file.path(output_dir, new_file)); old <- fread(old_file)
    merged <- merge(old[gene_symbol %in% candidates, .(gene_symbol, old_logFC = logFC, old_PValue = PValue, old_FDR = FDR)], new[gene_symbol %in% candidates, .(gene_symbol, new_logFC = logFC, new_PValue = PValue, new_FDR = FDR)], by = "gene_symbol", all = TRUE)
    if (nrow(merged)) { merged[, file := new_file]; candidate_rows[[length(candidate_rows) + 1L]] <- merged }
  } else {
    new <- fread(file.path(output_dir, new_file)); old <- fread(old_file)
    merged <- merge(old[, .(pathway, old_NES = NES, old_padj = padj)], new[, .(pathway, new_NES = NES, new_padj = padj)], by = "pathway", all = TRUE)
    merged[, file := new_file]; hallmark_rows[[length(hallmark_rows) + 1L]] <- merged
  }
}
candidate_out <- rbindlist(candidate_rows, fill = TRUE)
hallmark_out <- rbindlist(hallmark_rows, fill = TRUE)
fwrite(candidate_out, file.path(output_dir, "TYMS_EFNB2_LRRC17_sensitivity.tsv"), sep = "\t")
fwrite(hallmark_out, file.path(output_dir, "Hallmark_sensitivity_comparison.tsv"), sep = "\t")

composition <- cell_counts[condition %in% c("normal", "organic_ED_nonDM") & cell_type != "T_fibroblast_mixed", .(
  normal_mean_fraction = mean(cell_fraction[condition == "normal"]),
  ED_mean_fraction = mean(cell_fraction[condition == "organic_ED_nonDM"]),
  log2_fraction_ratio = log2((mean(cell_fraction[condition == "organic_ED_nonDM"]) + 1e-5) / (mean(cell_fraction[condition == "normal"]) + 1e-5)),
  p_value = tryCatch(t.test(cell_fraction ~ condition)$p.value, error = function(e) NA_real_)
), by = cell_type]
composition[, FDR := p.adjust(p_value, method = "BH")]
setorder(composition, FDR, p_value)
fwrite(composition, file.path(output_dir, "celltype_composition_nonDM_ED_vs_normal.tsv"), sep = "\t")

manifest <- data.table(parameter = c("analysis", "seed", "input_pseudobulk", "historical_reference", "library_edgeR", "library_fgsea", "library_msigdbr"), value = c("doublet-exclusion sensitivity; fixed 03A labels", 20260814, pseudobulk_dir, historical_dir, as.character(packageVersion("edgeR")), as.character(packageVersion("fgsea")), as.character(packageVersion("msigdbr"))))
fwrite(manifest, file.path(output_dir, "analysis_manifest.tsv"), sep = "\t")
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
message("Gate 03B edgeR/Hallmark sensitivity complete")
