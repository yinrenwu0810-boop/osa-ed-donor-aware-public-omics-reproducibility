args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript 09l_run_es01_mapping_multiverse_attempt03.R <vp_root>", call. = FALSE)

vp <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(vp, "../.."), winslash = "/", mustWork = TRUE)
attempt <- file.path(vp, "10_exploratory_hypoxia_sensitivity", "attempt_20260910_03")
.libPaths(c(file.path(vp, "02_env", "Rlib"), .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(fgsea); library(jsonlite); library(digest) })

out_names <- c("06_hypoxia_mapping_qc.tsv", "06_analysis_multiverse.tsv", "06_specification_concordance.tsv")
out_paths <- file.path(attempt, out_names)
if (any(file.exists(out_paths))) stop("Refusing to overwrite VP-ES01-06 attempt_20260910_03 outputs.", call. = FALSE)

scope_path <- file.path(attempt, "01_scope_freeze.json")
preflight_path <- file.path(attempt, "00_preflight.json")
scope <- fromJSON(scope_path, simplifyVector = FALSE)
if (!identical(scope$status, "FROZEN_CORRECTED_ES06_SPECIFICATION")) stop("Corrected ES06 scope is not frozen.", call. = FALSE)
if (!identical(digest(file = preflight_path, algo = "sha256"), scope$preflight_sha256)) stop("HOLD_INPUT_DRIFT: attempt preflight hash differs from scope freeze.", call. = FALSE)

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- sub("^--file=", "", script_arg[[1L]])
script_sha256 <- digest(file = script_path, algo = "sha256")
scope_sha256 <- digest(file = scope_path, algo = "sha256")

freeze_path <- file.path(vp, "01_protocol", "VP_G06_enrichment_rule_freeze_v1.json")
freeze <- fromJSON(freeze_path, simplifyVector = FALSE)
gmt_path <- file.path(root, "analysis", "data", "processed", "gene_sets", "Hallmark.gmt")
if (!identical(digest(file = gmt_path, algo = "sha256"), freeze$input_sha256[["analysis/data/processed/gene_sets/Hallmark.gmt"]])) stop("HOLD_INPUT_DRIFT: Hallmark GMT hash mismatch.", call. = FALSE)

universe <- fread(file.path(vp, "04_prepared", "gene_universe.tsv"))$gene
if (length(universe) != 1388L || anyDuplicated(universe)) stop("HOLD_POTENTIAL_IMPLEMENTATION_DEFECT: frozen gene universe is not 1,388 unique genes.", call. = FALSE)

read_gmt_raw <- function(path) {
  fields <- strsplit(readLines(path, warn = FALSE, encoding = "UTF-8"), "\t", fixed = TRUE)
  names_out <- vapply(fields, `[[`, character(1), 1L)
  if (anyDuplicated(names_out)) stop("HOLD_POTENTIAL_IMPLEMENTATION_DEFECT: duplicated Hallmark pathway names.", call. = FALSE)
  sets_raw <- lapply(fields, function(x) x[-c(1L, 2L)])
  names(sets_raw) <- names_out
  sets_raw
}

strict_rank <- function(dt, score_col) {
  x <- copy(dt)[is.finite(get(score_col))]
  setorderv(x, c(score_col, "gene"), c(-1L, 1L), na.last = TRUE)
  if (nrow(x) != 1388L || anyDuplicated(x$gene) || !setequal(x$gene, universe)) stop("HOLD_POTENTIAL_IMPLEMENTATION_DEFECT: rank input does not exactly match frozen universe.", call. = FALSE)
  setNames((nrow(x) - seq_len(nrow(x)) + 1) / nrow(x), x$gene)
}

aggregate_rank_vectors <- function(vectors, method) {
  genes <- sort(names(vectors[[1L]]))
  if (!all(vapply(vectors, function(x) setequal(names(x), genes), logical(1)))) stop("HOLD_POTENTIAL_IMPLEMENTATION_DEFECT: rank vector gene membership differs.", call. = FALSE)
  values <- do.call(cbind, lapply(vectors, function(x) unname(x[genes])))
  score <- if (identical(method, "median")) apply(values, 1L, median) else rowMeans(values)
  strict_rank(data.table(gene = genes, aggregate_score = score), "aggregate_score")
}

jaccard <- function(a, b) {
  aa <- if (is.na(a) || !nzchar(a)) character() else strsplit(a, ";", fixed = TRUE)[[1L]]
  bb <- if (is.na(b) || !nzchar(b)) character() else strsplit(b, ";", fixed = TRUE)[[1L]]
  den <- length(union(aa, bb))
  if (den == 0L) return(NA_real_)
  length(intersect(aa, bb)) / den
}

run_hallmark <- function(ranks, call_seed) {
  set.seed(call_seed)
  z <- suppressWarnings(fgseaMultilevel(pathways = hallmark_sets, stats = ranks, minSize = 15L, maxSize = 500L, eps = 0, scoreType = "std", nproc = 1L))
  z[, pathway := as.character(pathway)]
  z[, leading_edge := vapply(leadingEdge, paste, collapse = ";", character(1))]
  z[, leading_edge_size := lengths(leadingEdge)]
  z[, abs_NES_sort := -abs(NES)]
  setorder(z, padj, abs_NES_sort, pathway, na.last = TRUE)
  z[, abs_NES_sort := NULL]
  z[, hallmark_rank := seq_len(.N)]
  if (nrow(z) != 32L || !"HALLMARK_HYPOXIA" %in% z$pathway) stop("HOLD_POTENTIAL_IMPLEMENTATION_DEFECT: Hallmark tested family is not the frozen 32-set family.", call. = FALSE)
  z[pathway == "HALLMARK_HYPOXIA"]
}

raw_sets <- read_gmt_raw(gmt_path)
hypoxia_raw <- raw_sets[["HALLMARK_HYPOXIA"]]
hypoxia_unique <- unique(hypoxia_raw)
hallmark_sets <- lapply(raw_sets, unique)
hallmark_sets <- lapply(hallmark_sets, intersect, y = universe)
hallmark_sets <- hallmark_sets[lengths(hallmark_sets) >= 15L & lengths(hallmark_sets) <= 500L]
formal_eligible <- fread(file.path(vp, "07_enrichment", "v1", "gene_set_eligibility.tsv"))[database == "Hallmark" & status == "TESTED", pathway]
if (!identical(sort(names(hallmark_sets)), sort(formal_eligible))) stop("HOLD_POTENTIAL_IMPLEMENTATION_DEFECT: reconstructed Hallmark eligibility differs from VP-G06.", call. = FALSE)
formal_hypoxia_eligibility <- fread(file.path(vp, "07_enrichment", "v1", "gene_set_eligibility.tsv"))[database == "Hallmark" & pathway == "HALLMARK_HYPOXIA"]
if (nrow(formal_hypoxia_eligibility) != 1L || formal_hypoxia_eligibility$original_size != length(hypoxia_unique) || formal_hypoxia_eligibility$overlap_size != length(hallmark_sets[["HALLMARK_HYPOXIA"]])) stop("HOLD_POTENTIAL_IMPLEMENTATION_DEFECT: HALLMARK_HYPOXIA source/overlap sizes differ from VP-G06.", call. = FALSE)

registry <- fread(file.path(vp, "05_runs", "run_registry.tsv"))
setorder(registry, gene, condition, donor, seed)
if (nrow(registry) != 45L || any(registry$status != "PASS_TECHNICAL") || nrow(unique(registry[, .(gene, condition, donor, seed)])) != 45L) stop("HOLD_RUN_REGISTRY_MISMATCH: expected 45 unique PASS_TECHNICAL main runs.", call. = FALSE)
target_specs <- unique(registry[, .(target_gene = gene, condition)])
setorder(target_specs, target_gene)
if (nrow(target_specs) != 3L) stop("HOLD_RUN_REGISTRY_MISMATCH: target specification count is not three.", call. = FALSE)

consensus <- fread(file.path(vp, "06_consensus", "candidate_donor_consensus.tsv"))
all_rankings <- fread(file.path(vp, "06_consensus", "all_virtual_KO_rankings.tsv.gz"))
formal_full <- fread(file.path(vp, "07_enrichment", "v1", "target_consensus_enrichment.tsv"))[database == "Hallmark"]
formal_primary <- fread(file.path(vp, "07_enrichment", "v1", "hypoxia_primary_summary.tsv"))
formal_rank_inputs <- fread(file.path(vp, "07_enrichment", "v1", "target_consensus_rank_inputs.tsv"))
formal_donor <- fread(file.path(vp, "07_enrichment", "v1", "hypoxia_donor_enrichment.tsv"))[pathway == "HALLMARK_HYPOXIA"]
if (nrow(formal_primary) != 3L || nrow(formal_donor) != 9L) stop("HOLD_INPUT_DRIFT: formal VP-G06 hypoxia tables have unexpected cardinality.", call. = FALSE)

raw_seed_ranks <- list()
for (i in seq_len(nrow(registry))) {
  raw_path <- file.path(registry$output_attempt[[i]], "differential_regulation.tsv")
  if (!file.exists(raw_path) || file.info(raw_path)$size <= 0L) stop("HOLD_INPUT_DRIFT: missing raw run file.", call. = FALSE)
  d <- fread(raw_path, select = c("gene", "distance", "Z"))
  if (nrow(d) != 1388L || anyDuplicated(d$gene) || !setequal(d$gene, universe) || any(!is.finite(d$distance)) || any(!is.finite(d$Z))) stop("HOLD_POTENTIAL_IMPLEMENTATION_DEFECT: raw run schema/metrics invalid.", call. = FALSE)
  raw_seed_ranks[[registry$run_id[[i]]]] <- list(distance = strict_rank(d, "distance"), Z = strict_rank(d, "Z"))
}

mapping_rows <- list()
missing <- sort(setdiff(hypoxia_unique, universe))
case_only <- sort(hypoxia_unique[!(hypoxia_unique %in% universe) & toupper(hypoxia_unique) %in% toupper(universe)])
mapping_rows[[length(mapping_rows) + 1L]] <- data.table(
  row_type = "SOURCE_MAPPING", target_gene = NA_character_, condition = NA_character_, donor = NA_character_,
  source_member_count_raw = length(hypoxia_raw), source_member_count_unique = length(hypoxia_unique), duplicate_source_members = length(hypoxia_raw) - length(hypoxia_unique),
  universe_member_count = length(hallmark_sets[["HALLMARK_HYPOXIA"]]), universe_mapping_fraction = length(hallmark_sets[["HALLMARK_HYPOXIA"]]) / length(hypoxia_unique),
  missing_members = paste(missing, collapse = ";"), missing_member_reason = if (length(missing)) "not_in_frozen_1388_gene_universe" else "NONE",
  case_only_nonexact_members = paste(case_only, collapse = ";"), alias_or_retired_symbol_assessment = "NOT_ASSESSED_NO_FROZEN_ALIAS_REFERENCE_TABLE",
  detectable_member_count = NA_integer_, detectable_member_fraction = NA_real_, leading_edge_size = NA_integer_, degree_matched_members = NA_integer_,
  leading_edge_degree_median = NA_real_, leading_edge_degree_q1 = NA_real_, leading_edge_degree_q3 = NA_real_, leading_edge_degree_min = NA_real_, leading_edge_degree_max = NA_real_,
  leading_edge_rank_median = NA_real_, leading_edge_rank_q1 = NA_real_, leading_edge_rank_q3 = NA_real_, leading_edge_rank_min = NA_real_, leading_edge_rank_max = NA_real_,
  mapping_qc_status = if (length(case_only) || length(hypoxia_raw) != length(hypoxia_unique)) "HOLD_POTENTIAL_IMPLEMENTATION_DEFECT" else "PASS_EXACT_SYMBOL_MAPPING"
)
if (mapping_rows[[1L]]$mapping_qc_status != "PASS_EXACT_SYMBOL_MAPPING") stop("HOLD_POTENTIAL_IMPLEMENTATION_DEFECT: source duplicate or case-only mapping found.", call. = FALSE)

degree_ref <- fread(file.path(vp, "06_consensus", "negative_control_degree_reference", "degree_reference.tsv"))
if (nrow(degree_ref) != 6L * 1388L || anyDuplicated(degree_ref[, .(gene, condition, donor)])) stop("HOLD_INPUT_DRIFT: degree reference has unexpected structure.", call. = FALSE)
for (i in seq_len(nrow(target_specs))) {
  spec <- target_specs[i]
  donors <- sort(unique(registry[gene == spec$target_gene & condition == spec$condition, donor]))
  if (length(donors) != 3L) stop("HOLD_RUN_REGISTRY_MISMATCH: expected three donors per target.", call. = FALSE)
  for (donor_name in donors) {
    dref <- degree_ref[condition == spec$condition & donor == donor_name & gene %in% hallmark_sets[["HALLMARK_HYPOXIA"]]]
    if (nrow(dref) != length(hallmark_sets[["HALLMARK_HYPOXIA"]])) stop("HOLD_INPUT_DRIFT: degree reference fails HALLMARK_HYPOXIA coverage.", call. = FALSE)
    mapping_rows[[length(mapping_rows) + 1L]] <- data.table(
      row_type = "DONOR_DETECTABILITY", target_gene = spec$target_gene, condition = spec$condition, donor = donor_name,
      source_member_count_raw = NA_integer_, source_member_count_unique = NA_integer_, duplicate_source_members = NA_integer_, universe_member_count = length(hallmark_sets[["HALLMARK_HYPOXIA"]]), universe_mapping_fraction = NA_real_,
      missing_members = NA_character_, missing_member_reason = NA_character_, case_only_nonexact_members = NA_character_, alias_or_retired_symbol_assessment = NA_character_,
      detectable_member_count = sum(dref$detection_rate > 0), detectable_member_fraction = mean(dref$detection_rate > 0), leading_edge_size = NA_integer_, degree_matched_members = NA_integer_,
      leading_edge_degree_median = NA_real_, leading_edge_degree_q1 = NA_real_, leading_edge_degree_q3 = NA_real_, leading_edge_degree_min = NA_real_, leading_edge_degree_max = NA_real_,
      leading_edge_rank_median = NA_real_, leading_edge_rank_q1 = NA_real_, leading_edge_rank_q3 = NA_real_, leading_edge_rank_min = NA_real_, leading_edge_rank_max = NA_real_, mapping_qc_status = "PASS_FROZEN_DONOR_REFERENCE"
    )
    fd <- formal_donor[target_gene == spec$target_gene & condition == spec$condition & donor == donor_name]
    le <- if (nzchar(fd$leadingEdge[[1L]])) strsplit(fd$leadingEdge[[1L]], ";", fixed = TRUE)[[1L]] else character()
    ranks <- all_rankings[target_gene == spec$target_gene & condition == spec$condition & donor == donor_name, .(gene, median_standardized_rank)]
    joined <- merge(data.table(gene = le), dref[, .(gene, mean_outdegree_raw_10net)], by = "gene", all.x = TRUE)
    joined <- merge(joined, ranks, by = "gene", all.x = TRUE)
    if (anyNA(joined$mean_outdegree_raw_10net) || anyNA(joined$median_standardized_rank)) stop("HOLD_POTENTIAL_IMPLEMENTATION_DEFECT: formal donor leading edge cannot be mapped to degree/rank inputs.", call. = FALSE)
    mapping_rows[[length(mapping_rows) + 1L]] <- data.table(
      row_type = "LEADING_EDGE_DONOR_NETWORK_RANK_DISTRIBUTION", target_gene = spec$target_gene, condition = spec$condition, donor = donor_name,
      source_member_count_raw = NA_integer_, source_member_count_unique = NA_integer_, duplicate_source_members = NA_integer_, universe_member_count = NA_integer_, universe_mapping_fraction = NA_real_,
      missing_members = NA_character_, missing_member_reason = NA_character_, case_only_nonexact_members = NA_character_, alias_or_retired_symbol_assessment = NA_character_,
      detectable_member_count = NA_integer_, detectable_member_fraction = NA_real_, leading_edge_size = length(le), degree_matched_members = nrow(joined),
      leading_edge_degree_median = median(joined$mean_outdegree_raw_10net), leading_edge_degree_q1 = as.numeric(quantile(joined$mean_outdegree_raw_10net, 0.25, names = FALSE)), leading_edge_degree_q3 = as.numeric(quantile(joined$mean_outdegree_raw_10net, 0.75, names = FALSE)), leading_edge_degree_min = min(joined$mean_outdegree_raw_10net), leading_edge_degree_max = max(joined$mean_outdegree_raw_10net),
      leading_edge_rank_median = median(joined$median_standardized_rank), leading_edge_rank_q1 = as.numeric(quantile(joined$median_standardized_rank, 0.25, names = FALSE)), leading_edge_rank_q3 = as.numeric(quantile(joined$median_standardized_rank, 0.75, names = FALSE)), leading_edge_rank_min = min(joined$median_standardized_rank), leading_edge_rank_max = max(joined$median_standardized_rank), mapping_qc_status = "PASS_FROZEN_LEADING_EDGE_MAPPING"
    )
  }
}
mapping_qc <- rbindlist(mapping_rows, fill = TRUE)
if (nrow(mapping_qc) != 19L) stop("HOLD_POTENTIAL_IMPLEMENTATION_DEFECT: mapping QC row count is not 19.", call. = FALSE)

mv_rows <- list()
for (i in seq_len(nrow(target_specs))) {
  spec <- target_specs[i]
  target <- spec$target_gene
  condition <- spec$condition
  s0_input <- consensus[target_gene == target & condition == condition, .(gene, median_standardized_rank)]
  s0_ranks <- strict_rank(s0_input, "median_standardized_rank")
  saved_s0 <- formal_rank_inputs[target_gene == target & condition == condition, .(gene, gsea_score)]
  if (nrow(saved_s0) != 1388L || anyDuplicated(saved_s0$gene) || !setequal(saved_s0$gene, names(s0_ranks))) stop("FAIL_RECONCILIATION: formal S0 rank input differs in structure.", call. = FALSE)
  rank_input_max_abs_difference <- max(abs(s0_ranks[saved_s0$gene] - saved_s0$gsea_score))
  if (rank_input_max_abs_difference > 1e-15) stop("FAIL_RECONCILIATION: independently reconstructed S0 ranking differs from formal input.", call. = FALSE)
  formal_s0 <- formal_full[target_gene == target & condition == condition & pathway == "HALLMARK_HYPOXIA"]
  if (nrow(formal_s0) != 1L) stop("HOLD_INPUT_DRIFT: formal S0 Hallmark result missing.", call. = FALSE)
  s0 <- run_hallmark(s0_ranks, formal_s0$call_seed[[1L]])
  s0_es_diff <- abs(s0$ES - formal_s0$ES); s0_nes_diff <- abs(s0$NES - formal_s0$NES); s0_p_diff <- abs(s0$pval - formal_s0$pval); s0_fdr_diff <- abs(s0$padj - formal_s0$padj)
  s0_le_match <- identical(s0$leading_edge, formal_s0$leadingEdge)
  s0_status <- if (max(c(s0_es_diff, s0_nes_diff, s0_p_diff, s0_fdr_diff)) <= 1e-12 && s0_le_match) "PASS_NUMERIC_RECONCILIATION" else "FAIL_RECONCILIATION"
  if (s0_status != "PASS_NUMERIC_RECONCILIATION") stop("FAIL_RECONCILIATION: independently rerun S0 does not reproduce VP-G06.", call. = FALSE)
  mv_rows[[length(mv_rows) + 1L]] <- data.table(specification = "S0", evidence_label = "FORMAL_BASELINE", target_gene = target, condition = condition, aggregation_definition = "frozen distance plus five-seed median standardized rank plus three-donor consensus", ES = s0$ES, NES = s0$NES, p_value = s0$pval, FDR_within_Hallmark = s0$padj, leading_edge = s0$leading_edge, leading_edge_size = s0$leading_edge_size, leading_edge_definition = "S0 rank-level fgsea leading edge", Hallmark_tested_set_count = 32L, call_seed = formal_s0$call_seed, result_status = s0_status, rank_input_max_abs_difference = rank_input_max_abs_difference, S0_formal_ES = formal_s0$ES, S0_formal_NES = formal_s0$NES, S0_formal_p_value = formal_s0$pval, S0_formal_FDR = formal_s0$padj, S0_ES_abs_difference = s0_es_diff, S0_NES_abs_difference = s0_nes_diff, S0_p_value_abs_difference = s0_p_diff, S0_FDR_abs_difference = s0_fdr_diff, S0_leading_edge_exact_match = s0_le_match, scope_sha256 = scope_sha256, execution_script_sha256 = script_sha256)

  donors <- sort(unique(registry[gene == target & condition == condition, donor]))
  s1_donor <- list(); s2_donor <- list()
  for (donor_name in donors) {
    rr <- registry[gene == target & condition == condition & donor == donor_name]
    setorder(rr, seed)
    s1_donor[[donor_name]] <- aggregate_rank_vectors(lapply(rr$run_id, function(id) raw_seed_ranks[[id]]$distance), "mean")
    s2_donor[[donor_name]] <- aggregate_rank_vectors(lapply(rr$run_id, function(id) raw_seed_ranks[[id]]$Z), "median")
  }
  s1_ranks <- aggregate_rank_vectors(s1_donor, "mean")
  s2_ranks <- aggregate_rank_vectors(s2_donor, "median")
  for (item in list(list(name = "S1", ranks = s1_ranks, seed = as.integer(formal_s0$call_seed) + 1000L, definition = "distance plus five-seed mean standardized rank plus three-donor mean rank"), list(name = "S2", ranks = s2_ranks, seed = as.integer(formal_s0$call_seed) + 2000L, definition = "Z descending plus S0 median aggregation"))) {
    z <- run_hallmark(item$ranks, item$seed)
    mv_rows[[length(mv_rows) + 1L]] <- data.table(specification = item$name, evidence_label = "EXPLORATORY", target_gene = target, condition = condition, aggregation_definition = item$definition, ES = z$ES, NES = z$NES, p_value = z$pval, FDR_within_Hallmark = z$padj, leading_edge = z$leading_edge, leading_edge_size = z$leading_edge_size, leading_edge_definition = paste0(item$name, " rank-level fgsea leading edge"), Hallmark_tested_set_count = 32L, call_seed = item$seed, result_status = "EXPLORATORY_COMPLETE", rank_input_max_abs_difference = NA_real_, S0_formal_ES = formal_s0$ES, S0_formal_NES = formal_s0$NES, S0_formal_p_value = formal_s0$pval, S0_formal_FDR = formal_s0$padj, S0_ES_abs_difference = NA_real_, S0_NES_abs_difference = NA_real_, S0_p_value_abs_difference = NA_real_, S0_FDR_abs_difference = NA_real_, S0_leading_edge_exact_match = NA, scope_sha256 = scope_sha256, execution_script_sha256 = script_sha256)
  }

  donor_nes <- numeric(); donor_edges <- character()
  for (donor_name in donors) {
    d_input <- all_rankings[target_gene == target & condition == condition & donor == donor_name, .(gene, median_standardized_rank)]
    d_ranks <- strict_rank(d_input, "median_standardized_rank")
    fd <- formal_donor[target_gene == target & condition == condition & donor == donor_name]
    if (nrow(fd) != 1L) stop("HOLD_INPUT_DRIFT: formal donor result missing.", call. = FALSE)
    dz <- run_hallmark(d_ranks, fd$call_seed[[1L]])
    if (max(abs(dz$ES - fd$ES), abs(dz$NES - fd$NES), abs(dz$pval - fd$pval), abs(dz$padj - fd$padj)) > 1e-12 || !identical(dz$leading_edge, fd$leadingEdge)) stop("FAIL_RECONCILIATION: donor-first component does not reproduce VP-G06 donor result.", call. = FALSE)
    donor_nes[[donor_name]] <- dz$NES
    donor_edges[[donor_name]] <- dz$leading_edge
  }
  s3_edge_union <- paste(sort(unique(unlist(strsplit(donor_edges, ";", fixed = TRUE)))), collapse = ";")
  mv_rows[[length(mv_rows) + 1L]] <- data.table(specification = "S3", evidence_label = "EXPLORATORY", target_gene = target, condition = condition, aggregation_definition = "donor-first Hallmark fgsea followed by descriptive median donor NES", ES = NA_real_, NES = median(donor_nes), p_value = NA_real_, FDR_within_Hallmark = NA_real_, leading_edge = s3_edge_union, leading_edge_size = length(unique(unlist(strsplit(donor_edges, ";", fixed = TRUE)))), leading_edge_definition = "descriptive union of three donor-specific fgsea leading edges; no aggregate-rank inferential leading edge", Hallmark_tested_set_count = 32L, call_seed = NA_integer_, result_status = "EXPLORATORY_DESCRIPTIVE_NO_P_OR_FDR", rank_input_max_abs_difference = NA_real_, S0_formal_ES = formal_s0$ES, S0_formal_NES = formal_s0$NES, S0_formal_p_value = formal_s0$pval, S0_formal_FDR = formal_s0$padj, S0_ES_abs_difference = NA_real_, S0_NES_abs_difference = NA_real_, S0_p_value_abs_difference = NA_real_, S0_FDR_abs_difference = NA_real_, S0_leading_edge_exact_match = NA, scope_sha256 = scope_sha256, execution_script_sha256 = script_sha256)
}
multiverse <- rbindlist(mv_rows, fill = TRUE)
setorder(multiverse, target_gene, specification)
if (nrow(multiverse) != 12L || any(multiverse[specification %in% c("S0", "S1", "S2"), Hallmark_tested_set_count != 32L])) stop("HOLD_POTENTIAL_IMPLEMENTATION_DEFECT: multiverse output cardinality/FDR family failure.", call. = FALSE)

concordance <- rbindlist(lapply(unique(multiverse$target_gene), function(target) {
  z <- multiverse[target_gene == target]
  setorder(z, specification)
  pair_index <- combn(seq_len(nrow(z)), 2L)
  rbindlist(lapply(seq_len(ncol(pair_index)), function(i) {
    a <- z[pair_index[1L, i]]; b <- z[pair_index[2L, i]]
    data.table(target_gene = target, condition = a$condition, specification_a = a$specification, specification_b = b$specification, NES_a = a$NES, NES_b = b$NES, NES_b_minus_a = b$NES - a$NES, abs_NES_difference = abs(b$NES - a$NES), direction_concordant = sign(a$NES) == sign(b$NES), leading_edge_jaccard = jaccard(a$leading_edge, b$leading_edge), leading_edge_size_a = a$leading_edge_size, leading_edge_size_b = b$leading_edge_size, leading_edge_definition_a = a$leading_edge_definition, leading_edge_definition_b = b$leading_edge_definition, inference_comparability = if (a$specification == "S3" || b$specification == "S3") "S3_DESCRIPTIVE_NES_NO_P_OR_FDR" else "RANK_LEVEL_FGSEA_COMPARISON", scope_sha256 = scope_sha256, execution_script_sha256 = script_sha256)
  }))
}))
if (nrow(concordance) != 18L || any(!is.finite(concordance$leading_edge_jaccard))) stop("HOLD_POTENTIAL_IMPLEMENTATION_DEFECT: concordance rows/Jaccard values incomplete.", call. = FALSE)

tmp_paths <- file.path(attempt, paste0(out_names, ".tmp_", Sys.getpid()))
fwrite(mapping_qc, tmp_paths[[1L]], sep = "\t", na = "NA")
fwrite(multiverse, tmp_paths[[2L]], sep = "\t", na = "NA")
fwrite(concordance, tmp_paths[[3L]], sep = "\t", na = "NA")
if (!all(file.rename(tmp_paths, out_paths))) stop("Could not atomically publish all ES01-06 outputs.", call. = FALSE)
cat(toJSON(list(gate = "VP-ES01-06", attempt = "attempt_20260910_03", status = "PASS_COMPUTED_PENDING_INDEPENDENT_AUDIT", mapping_qc_rows = nrow(mapping_qc), multiverse_rows = nrow(multiverse), concordance_rows = nrow(concordance), S0_status = unique(multiverse[specification == "S0", result_status]), Hallmark_tested_set_count = 32L, execution_script_sha256 = script_sha256), auto_unbox = TRUE), "\n")
