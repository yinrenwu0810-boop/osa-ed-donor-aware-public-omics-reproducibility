args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript 07f_validate_g06_enrichment_ed_shared_compat.R <vp_root>")
source_path <- file.path(normalizePath(args[[1L]], winslash = "/", mustWork = TRUE), "scripts", "07a_validate_g06_enrichment.R")
source_lines <- readLines(source_path, warn = FALSE, encoding = "UTF-8")
old_source <- "ed_source[, in_vp_universe := gene_symbol %in% gene_universe]"
new_source <- c(old_source, "ed_shared_genes <- intersect(gene_universe, ed_source$gene_symbol)", "check(length(ed_shared_genes) == 1370L, \"ED shared universe is not 1370 genes\")")
old_rank <- "  rank_vector <- strict_rank(consensus[target_gene == spec$target_gene & condition == spec$condition, .(gene, median_standardized_rank)], \"median_standardized_rank\")"
new_rank <- c(old_rank, "  rank_vector <- rank_vector[names(rank_vector) %in% ed_shared_genes]", "  rank_vector <- strict_rank(data.table(gene = names(rank_vector), shared_score = as.numeric(rank_vector)), \"shared_score\")")
for (pair in list(list(old_source, new_source), list(old_rank, new_rank))) {
  index <- which(source_lines == pair[[1L]])
  if (length(index) != 1L) stop("Expected exactly one ED shared-universe audit expression")
  source_lines <- append(source_lines[-index], pair[[2L]], after = index - 1L)
}
eval(parse(text = source_lines, keep.source = TRUE), envir = new.env(parent = globalenv()))
