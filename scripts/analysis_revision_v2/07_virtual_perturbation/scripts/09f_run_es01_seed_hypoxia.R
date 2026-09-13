args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript 09f_run_es01_seed_hypoxia.R <vp_root>", call. = FALSE)
vp_root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(vp_root, "../.."), winslash = "/", mustWork = TRUE)
attempt <- file.path(vp_root, "10_exploratory_hypoxia_sensitivity", "attempt_20260910_02")
.libPaths(c(file.path(vp_root, "02_env", "Rlib"), .libPaths()))

suppressPackageStartupMessages({ library(data.table); library(fgsea); library(jsonlite); library(digest) })

preflight_path <- file.path(attempt, "00_preflight.json")
scope_path <- file.path(attempt, "01_scope_freeze.json")
qc_path <- file.path(attempt, "02_seed_rank_qc.tsv")
pair_path <- file.path(attempt, "02_seed_pair_stability.tsv")
recon_path <- file.path(attempt, "02_seed_aggregate_reconciliation.tsv")
outputs <- file.path(attempt, c("03_target_donor_seed_hypoxia.tsv", "03_target_donor_seed_hallmark_full.tsv.gz", "03_seed_hypoxia_summary.tsv"))
if (any(file.exists(outputs))) stop("Refusing to overwrite existing VP-ES01-03 output.", call. = FALSE)
if (!all(file.exists(c(preflight_path, scope_path, qc_path, pair_path, recon_path)))) stop("Required ES01-00--02 records are missing.", call. = FALSE)
preflight <- fromJSON(preflight_path, simplifyVector = FALSE)
scope <- fromJSON(scope_path, simplifyVector = FALSE)
if (!identical(preflight$status, "PASS_INPUTS_READ_ONLY") || length(preflight$errors) != 0L) stop("ES01 preflight is not a clean PASS.", call. = FALSE)
if (!identical(scope$status, "FROZEN_EXPLORATORY_SPECIFICATION") || !identical(scope$analysis_label, "EXPLORATORY_SENSITIVITY_NOT_CONFIRMATORY")) stop("ES01 scope is not frozen.", call. = FALSE)
if (!identical(scope$primary_ranking$metric, "distance") || !identical(scope$primary_ranking$tie_break, "gene symbol ascending")) stop("Frozen primary ranking rule mismatch.", call. = FALSE)
if (!identical(scope$primary_gene_set$formal_multiple_testing_family, "retain the original VP-G06 within-Hallmark FDR family; do not shrink it")) stop("Frozen FDR-family rule mismatch.", call. = FALSE)
qc <- fread(qc_path); pairs <- fread(pair_path); recon <- fread(recon_path)
if (nrow(qc) != 45L || any(qc$status != "PASS_SEED_RANK_QC") || nrow(pairs) != 90L || any(pairs$status != "PASS_RECOMPUTED") || nrow(recon) != 12492L || any(recon$comparison_status != "PASS_EXACT_RECONCILIATION")) stop("ES01-02 reconciliation is not a clean PASS.", call. = FALSE)

g06_freeze_path <- file.path(vp_root, "01_protocol", "VP_G06_enrichment_rule_freeze_v1.json")
g06 <- fromJSON(g06_freeze_path, simplifyVector = FALSE)
hallmark_path <- file.path(project_root, "analysis", "data", "processed", "gene_sets", "Hallmark.gmt")
expected_gmt_sha <- g06$input_sha256[["analysis/data/processed/gene_sets/Hallmark.gmt"]]
if (!identical(digest(file = hallmark_path, algo = "sha256"), expected_gmt_sha)) stop("HOLD_INPUT_DRIFT: frozen Hallmark GMT hash mismatch.", call. = FALSE)
if (!identical(as.character(packageVersion("fgsea")), as.character(g06$environment$fgsea))) stop("fgsea version differs from VP-G06 freeze.", call. = FALSE)

read_gmt <- function(path) {
  fields <- strsplit(readLines(path, warn = FALSE, encoding = "UTF-8"), "\t", fixed = TRUE)
  names_out <- vapply(fields, `[[`, character(1), 1L)
  if (length(names_out) != 50L || anyDuplicated(names_out)) stop("Frozen Hallmark GMT must contain 50 uniquely named sets.", call. = FALSE)
  sets <- lapply(fields, function(x) unique(x[-c(1L, 2L)])); names(sets) <- names_out; sets
}
strict_rank <- function(dt) {
  x <- copy(dt)[is.finite(distance)]
  if (nrow(x) != 1388L || anyDuplicated(x$gene)) stop("Seed ranking is not a unique 1,388-gene profile.", call. = FALSE)
  setorderv(x, c("distance", "gene"), c(-1L, 1L), na.last = TRUE)
  x[, gsea_score := (nrow(x) - seq_len(.N) + 1) / nrow(x)]
  setNames(x$gsea_score, x$gene)
}
run_fgsea <- function(pathways, ranks, call_seed) {
  set.seed(call_seed)
  result <- suppressWarnings(fgseaMultilevel(pathways = pathways, stats = ranks, minSize = 15L, maxSize = 500L, eps = 0, scoreType = "std", nproc = 1L))
  result[, leading_edge := vapply(leadingEdge, function(x) paste(x, collapse = ";"), character(1))]
  result[, leading_edge_size := lengths(leadingEdge)]
  result[, pathway := as.character(pathway)]
  setorder(result, padj, -abs(NES), pathway, na.last = TRUE)
  result[, hallmark_rank := seq_len(.N)]
  result
}

universe <- fread(file.path(vp_root, "04_prepared", "gene_universe.tsv"))$gene
if (length(universe) != 1388L || anyDuplicated(universe)) stop("Fixed gene universe mismatch.", call. = FALSE)
hallmark <- read_gmt(hallmark_path)
eligible <- lapply(hallmark, intersect, y = universe)
eligible <- eligible[lengths(eligible) >= 15L & lengths(eligible) <= 500L]
formal_eligibility <- fread(file.path(vp_root, "07_enrichment", "v1", "gene_set_eligibility.tsv"))[database == "Hallmark"]
formal_tested <- formal_eligibility[status == "TESTED", pathway]
if (!identical(sort(names(eligible)), sort(formal_tested)) || !"HALLMARK_HYPOXIA" %in% names(eligible)) stop("Frozen Hallmark eligibility differs from VP-G06.", call. = FALSE)

registry <- fread(file.path(vp_root, "05_runs", "run_registry.tsv"))
if (nrow(registry) != 45L || any(registry$status != "PASS_TECHNICAL")) stop("Main run registry mismatch.", call. = FALSE)
registry[, seed_order := as.integer(seed)]
setorder(registry, gene, condition, donor, seed_order)
registry[, seed_order := NULL]
expected_seeds <- as.character(2026082701:2026082705)
group_counts <- registry[, .(n = .N, seeds = paste(sort(as.character(seed)), collapse = ";")), by = .(gene, condition, donor)]
if (nrow(group_counts) != 9L || any(group_counts$n != 5L) || any(group_counts$seeds != paste(expected_seeds, collapse = ";"))) stop("Frozen 3x3x5 design mismatch.", call. = FALSE)

tmp_dir <- file.path(attempt, paste0("03_seed_hypoxia.tmp_", Sys.getpid()))
if (dir.exists(tmp_dir)) stop("Temporary ES01-03 directory already exists.", call. = FALSE)
dir.create(tmp_dir, recursive = FALSE)

seed_counter <- as.integer(g06$method$random_seed)
full <- list()
for (i in seq_len(nrow(registry))) {
  task <- registry[i]
  input_path <- file.path(task$output_attempt, "differential_regulation.tsv")
  input <- fread(input_path)
  if (!identical(names(input), c("gene", "distance", "Z", "FC", "p.value", "p.adj"))) stop("Differential-regulation schema mismatch: ", task$run_id, call. = FALSE)
  if (nrow(input) != 1388L || anyDuplicated(input$gene) || !setequal(input$gene, universe) || !(task$gene %in% input$gene) || any(!is.finite(input$distance))) stop("Seed input QC mismatch: ", task$run_id, call. = FALSE)
  ranks <- strict_rank(input[, .(gene, distance)])
  seed_counter <- seed_counter + 1L
  result <- run_fgsea(eligible, ranks, seed_counter)
  result[, `:=`(profile_type = "single_seed_network_rank", target_gene = task$gene, condition = task$condition, donor = task$donor, seed = as.character(task$seed), database = "Hallmark", call_seed = seed_counter)]
  setcolorder(result, c("profile_type", "target_gene", "condition", "donor", "seed", "database", "call_seed", "pathway", "pval", "padj", "log2err", "ES", "NES", "size", "leading_edge", "leading_edge_size", "hallmark_rank"))
  full[[length(full) + 1L]] <- result
}
full <- rbindlist(full, use.names = TRUE)
full[, seed_order := as.integer(seed)]
setorder(full, target_gene, condition, donor, seed_order, hallmark_rank)
full[, seed_order := NULL]
if (nrow(full) != 45L * length(eligible) || any(full[, .N, by = .(target_gene, condition, donor, seed)]$N != length(eligible))) stop("Full seed-level Hallmark output cardinality mismatch.", call. = FALSE)
hypoxia <- full[pathway == "HALLMARK_HYPOXIA", .(target_gene, condition, donor, seed, profile_type, database, call_seed, ES, NES, p_value = pval, FDR_within_Hallmark = padj, log2err, set_size = size, leading_edge, leading_edge_size, hallmark_rank)]
if (nrow(hypoxia) != 45L) stop("Expected exactly 45 seed-level HALLMARK_HYPOXIA rows.", call. = FALSE)
summary <- hypoxia[, .(seeds = .N, NES_median = median(NES), NES_min = min(NES), NES_max = max(NES), NES_range = max(NES) - min(NES), positive_NES_seed_count = sum(NES > 0), negative_NES_seed_count = sum(NES < 0), FDR_lt_0_05_seed_count = sum(FDR_within_Hallmark < 0.05, na.rm = TRUE), leading_edge_size_median = median(leading_edge_size)), by = .(target_gene, condition, donor)]
hypoxia[, seed_order := as.integer(seed)]
setorder(hypoxia, target_gene, condition, donor, seed_order)
hypoxia[, seed_order := NULL]
setorder(summary, target_gene, condition, donor)
fwrite(hypoxia, file.path(tmp_dir, "03_target_donor_seed_hypoxia.tsv"), sep = "\t")
fwrite(full, file.path(tmp_dir, "03_target_donor_seed_hallmark_full.tsv.gz"), sep = "\t", compress = "gzip")
fwrite(summary, file.path(tmp_dir, "03_seed_hypoxia_summary.tsv"), sep = "\t")
for (name in basename(outputs)) if (!file.rename(file.path(tmp_dir, name), file.path(attempt, name))) stop("Could not seal ES01-03 output: ", name, call. = FALSE)
unlink(tmp_dir, recursive = TRUE)
cat(toJSON(list(gate = "VP-ES01-03", status = "PASS_COMPUTED_PENDING_NEXT_GATE", seed_profiles = 45L, full_hallmark_rows = nrow(full), hypoxia_rows = nrow(hypoxia), call_seed_first = as.integer(g06$method$random_seed) + 1L, call_seed_last = seed_counter), auto_unbox = TRUE), "\n")
