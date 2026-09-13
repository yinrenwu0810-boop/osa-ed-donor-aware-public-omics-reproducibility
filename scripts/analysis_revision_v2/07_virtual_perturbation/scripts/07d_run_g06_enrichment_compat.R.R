args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript 07d_run_g06_enrichment_compat.R <vp_root>")
source_path <- file.path(normalizePath(args[[1L]], winslash = "/", mustWork = TRUE), "scripts", "07_run_g06_enrichment.R")
source_lines <- readLines(source_path, warn = FALSE, encoding = "UTF-8")
patched_lines <- gsub("freeze$random_seed", "freeze$method$random_seed", source_lines, fixed = TRUE)
if (sum(source_lines != patched_lines) != 2L) stop("Random-seed path compatibility substitution count is not two")
replacements <- list(
  "  setorder(result, padj, -abs(NES), pathway, na.last = TRUE)" = c(
    "  result[, abs_NES_order__ := abs(NES)]",
    "  setorder(result, padj, -abs_NES_order__, pathway, na.last = TRUE)",
    "  result[, abs_NES_order__ := NULL]"
  ),
  "setorder(main_results, target_gene, database, padj, -abs(NES), pathway, na.last = TRUE)" = c(
    "main_results[, abs_NES_order__ := abs(NES)]",
    "setorder(main_results, target_gene, database, padj, -abs_NES_order__, pathway, na.last = TRUE)",
    "main_results[, abs_NES_order__ := NULL]"
  ),
  "setorder(donor_hallmark, target_gene, donor, padj, -abs(NES), pathway, na.last = TRUE)" = c(
    "donor_hallmark[, abs_NES_order__ := abs(NES)]",
    "setorder(donor_hallmark, target_gene, donor, padj, -abs_NES_order__, pathway, na.last = TRUE)",
    "donor_hallmark[, abs_NES_order__ := NULL]"
  ),
  "setorder(primary_hallmark, target_gene, role, selection_order, padj, -abs(NES), pathway, na.last = TRUE)" = c(
    "primary_hallmark[, abs_NES_order__ := abs(NES)]",
    "setorder(primary_hallmark, target_gene, role, selection_order, padj, -abs_NES_order__, pathway, na.last = TRUE)",
    "primary_hallmark[, abs_NES_order__ := NULL]"
  ),
  "setorder(ed_gsea, target_gene, padj, -abs(NES), pathway, na.last = TRUE)" = c(
    "ed_gsea[, abs_NES_order__ := abs(NES)]",
    "setorder(ed_gsea, target_gene, padj, -abs_NES_order__, pathway, na.last = TRUE)",
    "ed_gsea[, abs_NES_order__ := NULL]"
  )
)
for (old in names(replacements)) {
  index <- which(patched_lines == old)
  if (length(index) != 1L) stop("Expected exactly one ordering expression: ", old)
  patched_lines <- append(patched_lines[-index], replacements[[old]], after = index - 1L)
}
eval(parse(text = patched_lines, keep.source = TRUE), envir = new.env(parent = globalenv()))
