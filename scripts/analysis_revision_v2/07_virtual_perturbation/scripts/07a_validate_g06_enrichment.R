args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript 07a_validate_g06_enrichment.R <vp_root>")
vp_root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(vp_root, "../.."), winslash = "/", mustWork = TRUE)
.libPaths(c(file.path(vp_root, "02_env", "Rlib"), .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(fgsea); library(jsonlite); library(digest) })

freeze_path <- file.path(vp_root, "01_protocol", "VP_G06_enrichment_rule_freeze_v1.json")
result_dir <- file.path(vp_root, "07_enrichment", "v1")
pass_path <- file.path(vp_root, "validation", "VP_G06_enrichment_AUDIT01.json")
fail_path <- file.path(vp_root, "validation", "VP_G06_enrichment_AUDIT01_FAIL_RETAINED.json")
if (file.exists(pass_path) || file.exists(fail_path)) stop("Refusing to overwrite an existing G06 audit")
freeze <- fromJSON(freeze_path)
errors <- character()
add_error <- function(message) errors <<- c(errors, message)
check <- function(value, message) if (!isTRUE(value)) add_error(message)

read_gmt <- function(path) {
  fields <- strsplit(readLines(path, warn = FALSE, encoding = "UTF-8"), "\t", fixed = TRUE)
  out <- lapply(fields, function(x) unique(x[-c(1L, 2L)]))
  names(out) <- vapply(fields, `[[`, character(1), 1L)
  out
}
strict_rank <- function(dt, score_col) {
  x <- copy(dt)[is.finite(get(score_col))]
  setorderv(x, c(score_col, "gene"), c(-1L, 1L), na.last = TRUE)
  x[, gsea_score := (nrow(x) - seq_len(.N) + 1) / nrow(x)]
  setNames(x$gsea_score, x$gene)
}
run_fgsea <- function(pathways, ranks, seed_value, min_size = 15L, max_size = 500L) {
  set.seed(seed_value)
  x <- suppressWarnings(fgseaMultilevel(pathways, ranks, minSize = min_size, maxSize = max_size, eps = 0, scoreType = "std", nproc = 1L))
  x[, leadingEdge := vapply(leadingEdge, function(z) paste(z, collapse = ";"), character(1))]
  x[, pathway := as.character(pathway)]
  x
}
compare_fgsea <- function(observed, expected, label) {
  check(setequal(observed$pathway, expected$pathway) && nrow(observed) == nrow(expected), paste(label, "pathway membership mismatch"))
  joined <- merge(observed, expected, by = "pathway", suffixes = c(".obs", ".exp"), all = TRUE)
  if (anyNA(joined$pathway)) { add_error(paste(label, "merge mismatch")); return(invisible(NULL)) }
  for (column in c("pval", "padj", "log2err", "ES", "NES", "size")) {
    a <- joined[[paste0(column, ".obs")]]; b <- joined[[paste0(column, ".exp")]]
    same <- (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) <= 1e-12 * pmax(1, abs(a), abs(b)))
    if (!all(same)) add_error(paste(label, column, "recomputation mismatch"))
  }
  if (!all(joined$leadingEdge.obs == joined$leadingEdge.exp)) add_error(paste(label, "leadingEdge recomputation mismatch"))
}
validate_rank_profile <- function(x, label, universe) {
  check(nrow(x) == 1388L && !anyDuplicated(x$gene) && setequal(x$gene, universe), paste(label, "gene universe mismatch"))
  check(setequal(round(x$gsea_score * 1388), 1:1388), paste(label, "strict rank-score grid mismatch"))
}

expected_files <- c(
  "gene_set_eligibility.tsv", "target_consensus_rank_inputs.tsv", "target_consensus_enrichment.tsv",
  "donor_hallmark_enrichment.tsv", "hypoxia_donor_enrichment.tsv", "hypoxia_donor_consistency.tsv",
  "matched_profile_hallmark_enrichment.tsv", "matched_profile_rank_inputs.tsv",
  "hypoxia_matched_control_enrichment.tsv", "hypoxia_control_calibration.tsv", "ed_full_rank.tsv",
  "ed_signature_sets.tsv", "ed_signature_gsea.tsv", "ed_top5_overlap.tsv", "ed_rank_concordance.tsv",
  "hypoxia_primary_summary.tsv", "g06_report.json", "sessionInfo.txt"
)
manifest <- fread(file.path(result_dir, "artifact_manifest.sha256.tsv"))
check(setequal(manifest$file, expected_files) && nrow(manifest) == length(expected_files), "Artifact manifest membership mismatch")
for (i in seq_len(nrow(manifest))) {
  path <- file.path(result_dir, manifest$file[i])
  check(file.exists(path), paste("Missing artifact", manifest$file[i]))
  if (file.exists(path)) {
    check(file.info(path)$size == manifest$bytes[i], paste("Byte count mismatch", manifest$file[i]))
    check(identical(digest(file = path, algo = "sha256"), manifest$sha256[i]), paste("SHA-256 mismatch", manifest$file[i]))
  }
}

gene_universe <- fread(file.path(vp_root, "04_prepared", "gene_universe.tsv"))$gene
check(length(gene_universe) == 1388L && !anyDuplicated(gene_universe), "Frozen gene universe mismatch")
gmt_paths <- list(
  Hallmark = file.path(project_root, "analysis", "data", "processed", "gene_sets", "Hallmark.gmt"),
  Reactome = file.path(project_root, "analysis", "data", "processed", "gene_sets", "Reactome.gmt"),
  GO_BP = file.path(project_root, "analysis", "data", "processed", "gene_sets", "GO_BP.gmt")
)
sets <- lapply(gmt_paths, read_gmt)
eligibility <- fread(file.path(result_dir, "gene_set_eligibility.tsv"))
expected_eligibility <- rbindlist(lapply(names(sets), function(db) {
  overlaps <- lengths(lapply(sets[[db]], intersect, y = gene_universe))
  data.table(database = db, pathway = names(sets[[db]]), original_size = lengths(sets[[db]]), overlap_size = overlaps,
             status = fifelse(overlaps >= 15L & overlaps <= 500L, "TESTED", "EXCLUDED_SIZE"))
}))
setorder(eligibility, database, pathway); setorder(expected_eligibility, database, pathway)
check(identical(eligibility, expected_eligibility), "Gene-set eligibility recomputation mismatch")
eligible <- lapply(names(sets), function(db) {
  keep <- expected_eligibility[database == db & status == "TESTED", pathway]
  lapply(sets[[db]][keep], intersect, y = gene_universe)
}); names(eligible) <- names(sets)

consensus <- fread(file.path(vp_root, "06_consensus", "candidate_donor_consensus.tsv"))
rankings <- fread(file.path(vp_root, "06_consensus", "all_virtual_KO_rankings.tsv.gz"))
target_rank_inputs <- fread(file.path(result_dir, "target_consensus_rank_inputs.tsv"))
main_results <- fread(file.path(result_dir, "target_consensus_enrichment.tsv"))
target_specs <- unique(consensus[, .(target_gene, condition)]); setorder(target_specs, target_gene)
for (i in seq_len(nrow(target_specs))) {
  spec <- target_specs[i]
  observed_rank <- target_rank_inputs[target_gene == spec$target_gene & condition == spec$condition]
  validate_rank_profile(observed_rank, paste("target", spec$target_gene), gene_universe)
  expected_rank <- strict_rank(consensus[target_gene == spec$target_gene & condition == spec$condition, .(gene, median_standardized_rank)], "median_standardized_rank")
  joined_rank <- merge(observed_rank[, .(gene, observed = gsea_score)], data.table(gene = names(expected_rank), expected = as.numeric(expected_rank)), by = "gene")
  check(nrow(joined_rank) == 1388L && max(abs(joined_rank$observed - joined_rank$expected)) < 1e-15, paste("Target rank reconstruction mismatch", spec$target_gene))
  for (db in names(eligible)) {
    observed <- main_results[target_gene == spec$target_gene & condition == spec$condition & database == db]
    check(length(unique(observed$call_seed)) == 1L, paste("Main call seed mismatch", spec$target_gene, db))
    expected <- run_fgsea(eligible[[db]], expected_rank, unique(observed$call_seed))
    compare_fgsea(observed, expected, paste("main", spec$target_gene, db))
  }
}

donor_results <- fread(file.path(result_dir, "donor_hallmark_enrichment.tsv"))
donor_keys <- unique(donor_results[, .(target_gene, condition, donor)])
for (key_i in seq_len(nrow(donor_keys))) {
  key <- donor_keys[key_i]
  source <- rankings[target_gene == key$target_gene & condition == key$condition & donor == key$donor, .(gene, median_standardized_rank)]
  rank_vector <- strict_rank(source, "median_standardized_rank")
  observed <- donor_results[target_gene == key$target_gene & condition == key$condition & donor == key$donor]
  expected <- run_fgsea(eligible$Hallmark, rank_vector, unique(observed$call_seed))
  compare_fgsea(observed, expected, paste("donor", key$target_gene, key$donor))
}
hypoxia_donor <- fread(file.path(result_dir, "hypoxia_donor_enrichment.tsv"))
check(nrow(hypoxia_donor) == 9L && all(hypoxia_donor$pathway == "HALLMARK_HYPOXIA"), "Donor hypoxia row count mismatch")

rank_one_file <- function(path) {
  x <- fread(path, select = c("gene", "distance")); strict_rank(x, "distance")
}
aggregate_files <- function(paths) {
  ranks <- lapply(paths, rank_one_file); genes <- Reduce(intersect, lapply(ranks, names))
  matrix_values <- vapply(ranks, function(x) x[genes], numeric(length(genes)))
  strict_rank(data.table(gene = genes, score = apply(matrix_values, 1L, median)), "score")
}
selection <- fread(file.path(vp_root, "06_consensus", "negative_control_selection_v2", "negative_control_selection.tsv"))
tasks <- fread(file.path(vp_root, "01_protocol", "control_execution_v2", "control_task_matrix.tsv"))
batches <- fread(file.path(vp_root, "01_protocol", "control_execution_v2", "control_batch_matrix.tsv"))
registry <- fread(file.path(vp_root, "05_runs", "run_registry.tsv"))
matched_ranks <- fread(file.path(result_dir, "matched_profile_rank_inputs.tsv"))
matched_results <- fread(file.path(result_dir, "matched_profile_hallmark_enrichment.tsv"))
profile_keys <- unique(matched_ranks[, .(role, target_gene, condition, profile_gene, selection_order)])
check(nrow(profile_keys) == 33L, "Matched rank profile count mismatch")
for (key_i in seq_len(nrow(profile_keys))) {
  key <- profile_keys[key_i]
  observed_rank <- matched_ranks[role == key$role & target_gene == key$target_gene & condition == key$condition & profile_gene == key$profile_gene]
  validate_rank_profile(observed_rank, paste("matched", key$target_gene, key$profile_gene), gene_universe)
  donor_names <- sort(unique(selection[target_gene == key$target_gene & condition == key$condition, donor]))
  if (key$role == "target") {
    source_rows <- registry[gene == key$target_gene & condition == key$condition & donor %in% donor_names & seed == freeze$primary_seed & status == "PASS_TECHNICAL"]
    source_paths <- file.path(source_rows$output_attempt, "differential_regulation.tsv")
  } else {
    source_rows <- tasks[control_gene == key$profile_gene & condition == key$condition & donor %in% donor_names & seed == freeze$primary_seed]
    source_rows <- source_rows[vapply(strsplit(matched_targets, ";", fixed = TRUE), function(x) key$target_gene %in% x, logical(1))]
    source_rows <- merge(source_rows, batches[, .(batch_order, output_batch)], by = "batch_order", all.x = TRUE, sort = FALSE)
    source_paths <- file.path(vp_root, source_rows$output_batch, "controls", key$profile_gene, "attempt_01", "differential_regulation.tsv")
  }
  check(length(source_paths) == 3L, paste("Matched source count mismatch", key$target_gene, key$profile_gene))
  expected_rank <- aggregate_files(source_paths)
  joined <- merge(observed_rank[, .(gene, observed = gsea_score)], data.table(gene = names(expected_rank), expected = as.numeric(expected_rank)), by = "gene")
  check(nrow(joined) == 1388L && max(abs(joined$observed - joined$expected)) < 1e-15, paste("Matched rank reconstruction mismatch", key$target_gene, key$profile_gene))
  observed <- matched_results[role == key$role & target_gene == key$target_gene & condition == key$condition & profile_gene == key$profile_gene]
  expected <- run_fgsea(eligible$Hallmark, expected_rank, unique(observed$call_seed))
  compare_fgsea(observed, expected, paste("matched", key$target_gene, key$profile_gene))
}

ed_source <- fread(file.path(project_root, "revision_v2", "01_work", "GSE206528_doublet_sensitivity", "attempt_02", "edgeR_Hallmark", "differential_fibroblast_nonDM_ED_vs_normal.tsv"))
ed_source[, in_vp_universe := gene_symbol %in% gene_universe]
ed_sets <- list(
  ED_UP_FDR = sort(ed_source[in_vp_universe == TRUE & FDR < 0.05 & logFC > 0, gene_symbol]),
  ED_DOWN_FDR = sort(ed_source[in_vp_universe == TRUE & FDR < 0.05 & logFC < 0, gene_symbol])
)
observed_ed_sets <- fread(file.path(result_dir, "ed_signature_sets.tsv"))
for (name in names(ed_sets)) check(identical(sort(observed_ed_sets[signature == name, gene]), ed_sets[[name]]), paste("ED signature mismatch", name))
ed_gsea <- fread(file.path(result_dir, "ed_signature_gsea.tsv"))
for (i in seq_len(nrow(target_specs))) {
  spec <- target_specs[i]
  rank_vector <- strict_rank(consensus[target_gene == spec$target_gene & condition == spec$condition, .(gene, median_standardized_rank)], "median_standardized_rank")
  observed <- ed_gsea[target_gene == spec$target_gene & condition == spec$condition]
  expected <- run_fgsea(ed_sets[lengths(ed_sets) >= 5L & lengths(ed_sets) <= 500L], rank_vector, unique(observed$call_seed), 5L, 500L)
  compare_fgsea(observed, expected, paste("ED", spec$target_gene))
}

report <- fromJSON(file.path(result_dir, "g06_report.json"))
check(identical(report$status, "PASS_COMPUTED_PENDING_INDEPENDENT_AUDIT"), "G06 report status mismatch")
check(identical(report$freeze_sha256, digest(file = freeze_path, algo = "sha256")), "G06 report freeze hash mismatch")
audit <- list(
  gate = "VP-G06", stage = "ordered_enrichment_and_ED_signature_overlap_independent_audit",
  status = if (length(errors)) "FAIL_RETAINED" else "PASS",
  created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  freeze_sha256 = digest(file = freeze_path, algo = "sha256"),
  result_manifest_sha256 = digest(file = file.path(result_dir, "artifact_manifest.sha256.tsv"), algo = "sha256"),
  target_database_profiles_recomputed = 9L, donor_hallmark_profiles_recomputed = 9L,
  matched_hallmark_profiles_recomputed = 33L, ED_signature_profiles_recomputed = 3L,
  artifact_count_audited = nrow(manifest), errors = as.list(errors),
  boundary = "Independent source-to-output reconstruction and fgseaMultilevel recomputation. Biological direction, treatment effect, and causality are not established."
)
destination <- if (length(errors)) fail_path else pass_path
write_json(audit, destination, pretty = TRUE, auto_unbox = TRUE)
cat(toJSON(audit, pretty = TRUE, auto_unbox = TRUE), "\n")
if (length(errors)) quit(status = 1L)
