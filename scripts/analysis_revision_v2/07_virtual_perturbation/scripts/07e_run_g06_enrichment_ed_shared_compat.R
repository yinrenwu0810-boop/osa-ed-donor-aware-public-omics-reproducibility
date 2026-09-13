args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript 07e_run_g06_enrichment_ed_shared_compat.R <vp_root>")
source_path <- file.path(normalizePath(args[[1L]], winslash = "/", mustWork = TRUE), "scripts", "07_run_g06_enrichment.R")
source_lines <- readLines(source_path, warn = FALSE, encoding = "UTF-8")
patched_lines <- gsub("freeze$random_seed", "freeze$method$random_seed", source_lines, fixed = TRUE)
if (sum(source_lines != patched_lines) != 2L) stop("Random-seed path compatibility substitution count is not two")
replacements <- list(
  "  setorder(result, padj, -abs(NES), pathway, na.last = TRUE)" = c("  result[, abs_NES_order__ := abs(NES)]", "  setorder(result, padj, -abs_NES_order__, pathway, na.last = TRUE)", "  result[, abs_NES_order__ := NULL]"),
  "setorder(main_results, target_gene, database, padj, -abs(NES), pathway, na.last = TRUE)" = c("main_results[, abs_NES_order__ := abs(NES)]", "setorder(main_results, target_gene, database, padj, -abs_NES_order__, pathway, na.last = TRUE)", "main_results[, abs_NES_order__ := NULL]"),
  "setorder(donor_hallmark, target_gene, donor, padj, -abs(NES), pathway, na.last = TRUE)" = c("donor_hallmark[, abs_NES_order__ := abs(NES)]", "setorder(donor_hallmark, target_gene, donor, padj, -abs_NES_order__, pathway, na.last = TRUE)", "donor_hallmark[, abs_NES_order__ := NULL]"),
  "setorder(primary_hallmark, target_gene, role, selection_order, padj, -abs(NES), pathway, na.last = TRUE)" = c("primary_hallmark[, abs_NES_order__ := abs(NES)]", "setorder(primary_hallmark, target_gene, role, selection_order, padj, -abs_NES_order__, pathway, na.last = TRUE)", "primary_hallmark[, abs_NES_order__ := NULL]"),
  "setorder(ed_gsea, target_gene, padj, -abs(NES), pathway, na.last = TRUE)" = c("ed_gsea[, abs_NES_order__ := abs(NES)]", "setorder(ed_gsea, target_gene, padj, -abs_NES_order__, pathway, na.last = TRUE)", "ed_gsea[, abs_NES_order__ := NULL]"),
  "ed_universe <- ed[in_vp_universe == TRUE]" = c("ed_shared_genes <- intersect(gene_universe, ed$gene_symbol)", "ed_universe <- ed[gene_symbol %in% ed_shared_genes]"),
  "if (nrow(ed_universe) != 1388L || !setequal(ed_universe$gene_symbol, gene_universe)) stop(\"ED and VP universes do not match\")" = "if (nrow(ed_universe) != 1370L || !setequal(ed_universe$gene_symbol, ed_shared_genes)) stop(\"Frozen ED shared universe is not 1370 genes\")",
  "  ranks <- strict_rank(profile[, .(gene, median_standardized_rank)], \"median_standardized_rank\")" = c("  ranks <- strict_rank(profile[, .(gene, median_standardized_rank)], \"median_standardized_rank\")", "  ranks <- ranks[names(ranks) %in% ed_shared_genes]", "  ranks <- strict_rank(data.table(gene = names(ranks), shared_score = as.numeric(ranks)), \"shared_score\")"),
  "  top_genes <- profile[is_2of3_top5pct_consensus == TRUE, gene]" = "  top_genes <- intersect(profile[is_2of3_top5pct_consensus == TRUE, gene], ed_shared_genes)",
  "      universe_size = length(gene_universe), target_consensus_size = length(top_genes)," = "      universe_size = length(ed_shared_genes), target_consensus_size = length(top_genes),",
  "      p_value = phyper(length(common) - 1L, length(signature), length(gene_universe) - length(signature), length(top_genes), lower.tail = FALSE)," = "      p_value = phyper(length(common) - 1L, length(signature), length(ed_shared_genes) - length(signature), length(top_genes), lower.tail = FALSE),",
  "  gene_universe_size = length(gene_universe)," = c("  gene_universe_size = length(gene_universe),", "  ED_shared_universe_size = length(ed_shared_genes),")
)
for (old in names(replacements)) {
  index <- which(patched_lines == old)
  if (length(index) != 1L) stop("Expected exactly one compatibility expression: ", old)
  patched_lines <- append(patched_lines[-index], replacements[[old]], after = index - 1L)
}
eval(parse(text = patched_lines, keep.source = TRUE), envir = new.env(parent = globalenv()))
