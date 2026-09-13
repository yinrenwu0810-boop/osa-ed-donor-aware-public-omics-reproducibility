args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: Rscript integrate_GSE206528_doublet_sensitivity_v2.R <ed_dir> <output_dir> <workspace_root>")
ed_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
out_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = FALSE)
root <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(file.path(root, "analysis", ".Rlib"), .libPaths()))
suppressPackageStartupMessages(library(data.table))

integration_dir <- file.path(root, "analysis", "results", "integration")
endo <- fread(file.path(integration_dir, "endothelial_IH_meta_analysis.tsv"))
files <- list.files(ed_dir, pattern = "^differential_.*_nonDM_ED_vs_normal\\.tsv$", full.names = TRUE)
if (!length(files)) stop("No doublet-exclusion nonDM ED differential tables found")
integrated <- rbindlist(lapply(files, function(path) {
  cell_type <- sub("^differential_", "", basename(path)); cell_type <- sub("_nonDM_ED_vs_normal\\.tsv$", "", cell_type)
  ed <- fread(path, select = c("gene_symbol", "logFC", "PValue", "FDR"))
  setnames(ed, c("logFC", "PValue", "FDR"), c("ED_logFC", "ED_p", "ED_FDR"))
  x <- merge(endo, ed, by = "gene_symbol", all = FALSE)
  x[, cell_type := cell_type]
  x[, IH_ED_direction_concordant := sign(endothelial_meta_z) == sign(ED_logFC)]
  x[, evidence_score := pmin(-log10(pmax(endothelial_meta_FDR, 1e-300)), 20) + pmin(-log10(pmax(ED_p, 1e-300)), 20)]
  x[, evidence_tier := fifelse(endothelial_meta_FDR < 0.05 & endothelial_sign_concordant & ED_FDR < 0.05 & IH_ED_direction_concordant, "A", fifelse(endothelial_meta_FDR < 0.05 & endothelial_sign_concordant & ED_p < 0.05 & IH_ED_direction_concordant, "B", fifelse(endothelial_meta_FDR < 0.1 & endothelial_sign_concordant & ED_p < 0.1 & IH_ED_direction_concordant, "C", "unranked")))]
  x
}), fill = TRUE)
setorder(integrated, evidence_tier, -evidence_score)
fwrite(integrated, file.path(out_dir, "IH_ED_celltype_gene_evidence_doublet_sensitivity.tsv"), sep = "\t")
candidates <- integrated[evidence_tier %in% c("A", "B", "C")]
candidate_summary <- candidates[, .(best_tier = min(evidence_tier), best_score = max(evidence_score), supported_cell_types = paste(sort(unique(cell_type)), collapse = ";"), n_supported_cell_types = uniqueN(cell_type), endothelial_meta_z = first(endothelial_meta_z), endothelial_meta_FDR = first(endothelial_meta_FDR)), by = gene_symbol]
setorder(candidate_summary, best_tier, -n_supported_cell_types, -best_score)
fwrite(candidate_summary, file.path(out_dir, "IH_ED_prioritized_genes_doublet_sensitivity.tsv"), sep = "\t")

bulk_hallmark <- fread(file.path(integration_dir, "Hallmark_cross_dataset_evidence.tsv"))
gsea_files <- list.files(ed_dir, pattern = "^GSEA_Hallmark_.*_nonDM_ED_vs_normal\\.tsv$", full.names = TRUE)
hallmark_ed <- rbindlist(lapply(gsea_files, function(path) {
  cell_type <- sub("^GSEA_Hallmark_", "", basename(path)); cell_type <- sub("_nonDM_ED_vs_normal\\.tsv$", "", cell_type)
  x <- fread(path, select = c("pathway", "NES", "padj")); x[, cell_type := cell_type]; x
}))
hallmark_integrated <- merge(bulk_hallmark, hallmark_ed, by = "pathway", allow.cartesian = TRUE)
hallmark_integrated[, IH_ED_direction_concordant := sign(endothelial_NES_mean) == sign(NES)]
hallmark_integrated[, cross_context_significant := endothelial_direction_concordant & endothelial_max_FDR < 0.05 & padj < 0.05 & IH_ED_direction_concordant]
setorder(hallmark_integrated, -cross_context_significant, padj, endothelial_max_FDR)
fwrite(hallmark_integrated, file.path(out_dir, "IH_ED_Hallmark_celltype_evidence_doublet_sensitivity.tsv"), sep = "\t")

historical_candidates <- fread(file.path(integration_dir, "IH_ED_prioritized_genes.tsv"))
focus <- c("TYMS", "EFNB2", "LRRC17")
focus_table <- merge(historical_candidates[gene_symbol %in% focus, .(gene_symbol, historical_best_tier = best_tier, historical_supported_cell_types = supported_cell_types, historical_n_supported_cell_types = n_supported_cell_types)], candidate_summary[gene_symbol %in% focus, .(gene_symbol, sensitivity_best_tier = best_tier, sensitivity_supported_cell_types = supported_cell_types, sensitivity_n_supported_cell_types = n_supported_cell_types)], by = "gene_symbol", all = TRUE)
fwrite(focus_table, file.path(out_dir, "TYMS_EFNB2_LRRC17_integrated_sensitivity.tsv"), sep = "\t")
summary <- data.table(metric = c("tested_gene_celltype_rows", "ranked_candidate_genes", "tier_A_rows", "tier_A_candidate_genes", "cross_context_Hallmark_significant"), value = c(nrow(integrated), nrow(candidate_summary), sum(integrated$evidence_tier == "A"), sum(candidate_summary$best_tier == "A"), sum(hallmark_integrated$cross_context_significant)))
fwrite(summary, file.path(out_dir, "integration_summary.tsv"), sep = "\t")
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
message("Gate 03B candidate integration sensitivity complete")
