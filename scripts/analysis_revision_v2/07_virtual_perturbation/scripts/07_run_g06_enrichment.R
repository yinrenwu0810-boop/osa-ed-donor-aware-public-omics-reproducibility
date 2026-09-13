args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript 07_run_g06_enrichment.R <vp_root>")
vp_root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(vp_root, "../.."), winslash = "/", mustWork = TRUE)
.libPaths(c(file.path(vp_root, "02_env", "Rlib"), .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(jsonlite)
  library(digest)
})

freeze_path <- file.path(vp_root, "01_protocol", "VP_G06_enrichment_rule_freeze_v1.json")
control_audit_path <- file.path(vp_root, "validation", "VP_G05_control_runs_v2_AUDIT01.json")
out_dir <- file.path(vp_root, "07_enrichment", "v1")
if (dir.exists(out_dir)) stop("Refusing to overwrite sealed output directory: ", out_dir)
freeze <- fromJSON(freeze_path, simplifyVector = FALSE)
if (!identical(freeze$status, "FROZEN_BEFORE_ENRICHMENT")) stop("G06 rule is not frozen")
if (!identical(digest(file = control_audit_path, algo = "sha256"), freeze$control_run_audit_sha256)) stop("Control audit hash drift")
for (script_name in names(freeze$implementation_sha256)) {
  observed <- digest(file = file.path(vp_root, script_name), algo = "sha256")
  if (!identical(observed, freeze$implementation_sha256[[script_name]])) stop("Implementation hash drift: ", script_name)
}

tmp_dir <- paste0(out_dir, ".tmp_", Sys.getpid())
if (dir.exists(tmp_dir)) stop("Temporary output already exists: ", tmp_dir)
dir.create(tmp_dir, recursive = TRUE)

read_gmt <- function(path) {
  fields <- strsplit(readLines(path, warn = FALSE, encoding = "UTF-8"), "\t", fixed = TRUE)
  names_out <- vapply(fields, `[[`, character(1), 1L)
  if (anyDuplicated(names_out)) stop("Duplicate pathway names in ", path)
  sets <- lapply(fields, function(x) unique(x[-c(1L, 2L)]))
  names(sets) <- names_out
  sets
}

strict_rank <- function(dt, score_col) {
  x <- copy(dt)[is.finite(get(score_col))]
  setorderv(x, c(score_col, "gene"), c(-1L, 1L), na.last = TRUE)
  if (anyDuplicated(x$gene)) stop("Duplicate genes in rank input")
  x[, gsea_score := (nrow(x) - seq_len(.N) + 1) / nrow(x)]
  setNames(x$gsea_score, x$gene)
}

run_fgsea <- function(pathways, ranks, seed_value, min_size = 15L, max_size = 500L) {
  set.seed(seed_value)
  result <- suppressWarnings(fgseaMultilevel(
    pathways = pathways, stats = ranks, minSize = min_size, maxSize = max_size,
    eps = 0, scoreType = "std", nproc = 1L
  ))
  result[, leadingEdge := vapply(leadingEdge, function(x) paste(x, collapse = ";"), character(1))]
  result[, pathway := as.character(pathway)]
  setorder(result, padj, -abs(NES), pathway, na.last = TRUE)
  result
}

add_meta <- function(result, fields) {
  for (name in rev(names(fields))) result[, (name) := fields[[name]]]
  setcolorder(result, c(names(fields), setdiff(names(result), names(fields))))
  result
}

gene_universe <- fread(file.path(vp_root, "04_prepared", "gene_universe.tsv"))$gene
if (length(gene_universe) != 1388L || anyDuplicated(gene_universe)) stop("Frozen gene universe is not 1388 unique genes")
rankings <- fread(file.path(vp_root, "06_consensus", "all_virtual_KO_rankings.tsv.gz"))
consensus <- fread(file.path(vp_root, "06_consensus", "candidate_donor_consensus.tsv"))
selection <- fread(file.path(vp_root, "06_consensus", "negative_control_selection_v2", "negative_control_selection.tsv"))
tasks <- fread(file.path(vp_root, "01_protocol", "control_execution_v2", "control_task_matrix.tsv"))
batches <- fread(file.path(vp_root, "01_protocol", "control_execution_v2", "control_batch_matrix.tsv"))
registry <- fread(file.path(vp_root, "05_runs", "run_registry.tsv"))

gmt_paths <- list(
  Hallmark = file.path(project_root, "analysis", "data", "processed", "gene_sets", "Hallmark.gmt"),
  Reactome = file.path(project_root, "analysis", "data", "processed", "gene_sets", "Reactome.gmt"),
  GO_BP = file.path(project_root, "analysis", "data", "processed", "gene_sets", "GO_BP.gmt")
)
gene_sets <- lapply(gmt_paths, read_gmt)
eligibility <- rbindlist(lapply(names(gene_sets), function(database_name) {
  sets <- gene_sets[[database_name]]
  data.table(
    database = database_name,
    pathway = names(sets),
    original_size = lengths(sets),
    overlap_size = lengths(lapply(sets, intersect, y = gene_universe)),
    status = fifelse(lengths(lapply(sets, intersect, y = gene_universe)) >= 15L & lengths(lapply(sets, intersect, y = gene_universe)) <= 500L, "TESTED", "EXCLUDED_SIZE")
  )
}))
setorder(eligibility, database, pathway)
fwrite(eligibility, file.path(tmp_dir, "gene_set_eligibility.tsv"), sep = "\t")
eligible_sets <- lapply(names(gene_sets), function(database_name) {
  keep <- eligibility[database == database_name & status == "TESTED", pathway]
  lapply(gene_sets[[database_name]][keep], intersect, y = gene_universe)
})
names(eligible_sets) <- names(gene_sets)

target_specs <- unique(consensus[, .(target_gene, condition)])
setorder(target_specs, target_gene)
seed_counter <- as.integer(freeze$random_seed)
main_results <- list()
target_rank_inputs <- list()
for (i in seq_len(nrow(target_specs))) {
  spec <- target_specs[i]
  profile <- consensus[target_gene == spec$target_gene & condition == spec$condition, .(gene, median_standardized_rank)]
  ranks <- strict_rank(profile, "median_standardized_rank")
  if (!setequal(names(ranks), gene_universe)) stop("Consensus rank universe mismatch: ", spec$target_gene)
  target_rank_inputs[[length(target_rank_inputs) + 1L]] <- data.table(
    target_gene = spec$target_gene, condition = spec$condition, gene = names(ranks), gsea_score = as.numeric(ranks)
  )
  for (database_name in names(eligible_sets)) {
    seed_counter <- seed_counter + 1L
    result <- run_fgsea(eligible_sets[[database_name]], ranks, seed_counter)
    main_results[[length(main_results) + 1L]] <- add_meta(result, list(
      profile_type = "five_seed_three_donor_consensus", target_gene = spec$target_gene,
      condition = spec$condition, database = database_name, call_seed = seed_counter
    ))
  }
}
target_rank_inputs <- rbindlist(target_rank_inputs)
setorder(target_rank_inputs, target_gene, -gsea_score, gene)
fwrite(target_rank_inputs, file.path(tmp_dir, "target_consensus_rank_inputs.tsv"), sep = "\t")
main_results <- rbindlist(main_results, fill = TRUE)
setorder(main_results, target_gene, database, padj, -abs(NES), pathway, na.last = TRUE)
fwrite(main_results, file.path(tmp_dir, "target_consensus_enrichment.tsv"), sep = "\t")

donor_hallmark <- list()
for (i in seq_len(nrow(target_specs))) {
  spec <- target_specs[i]
  donor_names <- sort(unique(rankings[target_gene == spec$target_gene & condition == spec$condition, donor]))
  for (donor_name in donor_names) {
    profile <- rankings[target_gene == spec$target_gene & condition == spec$condition & donor == donor_name, .(gene, median_standardized_rank)]
    ranks <- strict_rank(profile, "median_standardized_rank")
    seed_counter <- seed_counter + 1L
    result <- run_fgsea(eligible_sets$Hallmark, ranks, seed_counter)
    donor_hallmark[[length(donor_hallmark) + 1L]] <- add_meta(result, list(
      profile_type = "five_seed_donor_consensus", target_gene = spec$target_gene,
      condition = spec$condition, donor = donor_name, database = "Hallmark", call_seed = seed_counter
    ))
  }
}
donor_hallmark <- rbindlist(donor_hallmark, fill = TRUE)
setorder(donor_hallmark, target_gene, donor, padj, -abs(NES), pathway, na.last = TRUE)
fwrite(donor_hallmark, file.path(tmp_dir, "donor_hallmark_enrichment.tsv"), sep = "\t")
hypoxia_donor <- donor_hallmark[pathway == "HALLMARK_HYPOXIA"]
fwrite(hypoxia_donor, file.path(tmp_dir, "hypoxia_donor_enrichment.tsv"), sep = "\t")
hypoxia_consistency <- hypoxia_donor[, .(
  donors = .N,
  positive_NES_donors = sum(NES > 0),
  negative_NES_donors = sum(NES < 0),
  FDR_lt_0_05_donors = sum(padj < 0.05, na.rm = TRUE),
  median_NES = median(NES),
  NES_min = min(NES),
  NES_max = max(NES)
), by = .(target_gene, condition)]
setorder(hypoxia_consistency, target_gene)
fwrite(hypoxia_consistency, file.path(tmp_dir, "hypoxia_donor_consistency.tsv"), sep = "\t")

rank_one_file <- function(path) {
  dt <- fread(path, select = c("gene", "distance"))
  strict_rank(dt, "distance")
}
aggregate_files <- function(paths) {
  donor_ranks <- lapply(paths, rank_one_file)
  genes <- Reduce(intersect, lapply(donor_ranks, names))
  if (length(genes) != 1388L) stop("Primary-seed donor rank universe mismatch")
  values <- vapply(donor_ranks, function(x) x[genes], numeric(length(genes)))
  aggregated <- data.table(gene = genes, donor_median_rank = apply(values, 1L, median))
  strict_rank(aggregated, "donor_median_rank")
}

primary_hallmark <- list()
matched_rank_inputs <- list()
controls_by_target <- unique(selection[, .(target_gene, condition, control_gene, selection_order)])
setorder(controls_by_target, target_gene, selection_order)
for (i in seq_len(nrow(target_specs))) {
  spec <- target_specs[i]
  donor_names <- sort(unique(selection[target_gene == spec$target_gene & condition == spec$condition, donor]))
  direct <- registry[gene == spec$target_gene & condition == spec$condition & donor %in% donor_names & seed == freeze$primary_seed & status == "PASS_TECHNICAL"]
  if (nrow(direct) != 3L) stop("Expected three primary target runs: ", spec$target_gene)
  direct_paths <- file.path(direct$output_attempt, "differential_regulation.tsv")
  ranks <- aggregate_files(direct_paths)
  matched_rank_inputs[[length(matched_rank_inputs) + 1L]] <- data.table(
    role = "target", target_gene = spec$target_gene, condition = spec$condition,
    profile_gene = spec$target_gene, selection_order = 0L, gene = names(ranks), gsea_score = as.numeric(ranks)
  )
  seed_counter <- seed_counter + 1L
  result <- run_fgsea(eligible_sets$Hallmark, ranks, seed_counter)
  primary_hallmark[[length(primary_hallmark) + 1L]] <- add_meta(result, list(
    profile_type = "primary_seed_three_donor_consensus", role = "target", target_gene = spec$target_gene,
    condition = spec$condition, profile_gene = spec$target_gene, selection_order = 0L,
    database = "Hallmark", call_seed = seed_counter
  ))
  selected_controls <- controls_by_target[target_gene == spec$target_gene & condition == spec$condition]
  for (j in seq_len(nrow(selected_controls))) {
    control <- selected_controls[j]
    task_rows <- tasks[control_gene == control$control_gene & condition == spec$condition & donor %in% donor_names & seed == freeze$primary_seed]
    task_rows <- task_rows[vapply(strsplit(matched_targets, ";", fixed = TRUE), function(x) spec$target_gene %in% x, logical(1))]
    if (nrow(task_rows) != 3L) stop("Expected three matched-control runs: ", spec$target_gene, "/", control$control_gene)
    task_rows <- merge(task_rows, batches[, .(batch_order, output_batch)], by = "batch_order", all.x = TRUE, sort = FALSE)
    control_paths <- file.path(vp_root, task_rows$output_batch, "controls", control$control_gene, "attempt_01", "differential_regulation.tsv")
    control_ranks <- aggregate_files(control_paths)
    matched_rank_inputs[[length(matched_rank_inputs) + 1L]] <- data.table(
      role = "control", target_gene = spec$target_gene, condition = spec$condition,
      profile_gene = control$control_gene, selection_order = control$selection_order,
      gene = names(control_ranks), gsea_score = as.numeric(control_ranks)
    )
    seed_counter <- seed_counter + 1L
    control_result <- run_fgsea(eligible_sets$Hallmark, control_ranks, seed_counter)
    primary_hallmark[[length(primary_hallmark) + 1L]] <- add_meta(control_result, list(
      profile_type = "primary_seed_three_donor_consensus", role = "control", target_gene = spec$target_gene,
      condition = spec$condition, profile_gene = control$control_gene, selection_order = control$selection_order,
      database = "Hallmark", call_seed = seed_counter
    ))
  }
}
matched_rank_inputs <- rbindlist(matched_rank_inputs)
setorder(matched_rank_inputs, target_gene, role, selection_order, -gsea_score, gene)
fwrite(matched_rank_inputs, file.path(tmp_dir, "matched_profile_rank_inputs.tsv"), sep = "\t")
primary_hallmark <- rbindlist(primary_hallmark, fill = TRUE)
setorder(primary_hallmark, target_gene, role, selection_order, padj, -abs(NES), pathway, na.last = TRUE)
fwrite(primary_hallmark, file.path(tmp_dir, "matched_profile_hallmark_enrichment.tsv"), sep = "\t")
hypoxia_profiles <- primary_hallmark[pathway == "HALLMARK_HYPOXIA"]
fwrite(hypoxia_profiles, file.path(tmp_dir, "hypoxia_matched_control_enrichment.tsv"), sep = "\t")
hypoxia_calibration <- rbindlist(lapply(sort(unique(hypoxia_profiles$target_gene)), function(target_name) {
  x <- hypoxia_profiles[target_gene == target_name]
  target_row <- x[role == "target"]
  control_rows <- x[role == "control"]
  data.table(
    target_gene = target_name,
    condition = target_row$condition,
    target_NES = target_row$NES,
    target_padj_within_Hallmark = target_row$padj,
    matched_controls = nrow(control_rows),
    signed_NES_empirical_percentile = mean(control_rows$NES <= target_row$NES),
    absolute_NES_empirical_percentile = mean(abs(control_rows$NES) <= abs(target_row$NES)),
    conservative_absolute_upper_tail_value = (1 + sum(abs(control_rows$NES) >= abs(target_row$NES))) / (nrow(control_rows) + 1L)
  )
}))
setorder(hypoxia_calibration, target_gene)
fwrite(hypoxia_calibration, file.path(tmp_dir, "hypoxia_control_calibration.tsv"), sep = "\t")

ed_path <- file.path(project_root, "revision_v2", "01_work", "GSE206528_doublet_sensitivity", "attempt_02", "edgeR_Hallmark", "differential_fibroblast_nonDM_ED_vs_normal.tsv")
ed <- fread(ed_path)
ed[, signed_stat := sqrt(pmax(F, 0)) * sign(logFC)]
ed[, absolute_stat := abs(signed_stat)]
ed[, in_vp_universe := gene_symbol %in% gene_universe]
setorder(ed, -absolute_stat, gene_symbol)
fwrite(ed, file.path(tmp_dir, "ed_full_rank.tsv"), sep = "\t")
ed_universe <- ed[in_vp_universe == TRUE]
if (nrow(ed_universe) != 1388L || !setequal(ed_universe$gene_symbol, gene_universe)) stop("ED and VP universes do not match")
ed_sets <- list(
  ED_UP_FDR = sort(ed_universe[FDR < 0.05 & logFC > 0, gene_symbol]),
  ED_DOWN_FDR = sort(ed_universe[FDR < 0.05 & logFC < 0, gene_symbol])
)
ed_set_rows <- rbindlist(lapply(names(ed_sets), function(set_name) data.table(signature = set_name, gene = ed_sets[[set_name]])))
fwrite(ed_set_rows, file.path(tmp_dir, "ed_signature_sets.tsv"), sep = "\t")

ed_gsea <- list(); ed_overlap <- list(); ed_concordance <- list()
for (i in seq_len(nrow(target_specs))) {
  spec <- target_specs[i]
  profile <- consensus[target_gene == spec$target_gene & condition == spec$condition, .(gene, median_standardized_rank, is_2of3_top5pct_consensus)]
  ranks <- strict_rank(profile[, .(gene, median_standardized_rank)], "median_standardized_rank")
  eligible_ed <- ed_sets[lengths(ed_sets) >= 5L & lengths(ed_sets) <= 500L]
  seed_counter <- seed_counter + 1L
  if (length(eligible_ed)) {
    result <- run_fgsea(eligible_ed, ranks, seed_counter, min_size = 5L, max_size = 500L)
    ed_gsea[[length(ed_gsea) + 1L]] <- add_meta(result, list(
      target_gene = spec$target_gene, condition = spec$condition, database = "ED_FDR_SIGNATURES", call_seed = seed_counter
    ))
  }
  top_genes <- profile[is_2of3_top5pct_consensus == TRUE, gene]
  overlap_part <- rbindlist(lapply(names(ed_sets), function(set_name) {
    signature <- ed_sets[[set_name]]
    common <- sort(intersect(top_genes, signature))
    data.table(
      target_gene = spec$target_gene, condition = spec$condition, signature = set_name,
      universe_size = length(gene_universe), target_consensus_size = length(top_genes),
      signature_size = length(signature), overlap_size = length(common),
      p_value = phyper(length(common) - 1L, length(signature), length(gene_universe) - length(signature), length(top_genes), lower.tail = FALSE),
      overlap_genes = paste(common, collapse = ";")
    )
  }))
  overlap_part[, FDR_within_target := p.adjust(p_value, method = "BH")]
  ed_overlap[[length(ed_overlap) + 1L]] <- overlap_part
  target_scores <- consensus[target_gene == spec$target_gene & condition == spec$condition, .(gene, perturbation_score = median_standardized_rank)]
  joined <- merge(target_scores, ed_universe[, .(gene = gene_symbol, ed_absolute_stat = absolute_stat)], by = "gene")
  ed_concordance[[length(ed_concordance) + 1L]] <- data.table(
    target_gene = spec$target_gene, condition = spec$condition, shared_genes = nrow(joined),
    spearman_rho_perturbation_vs_ED_absolute_stat = cor(joined$perturbation_score, joined$ed_absolute_stat, method = "spearman")
  )
}
ed_gsea <- rbindlist(ed_gsea, fill = TRUE)
setorder(ed_gsea, target_gene, padj, -abs(NES), pathway, na.last = TRUE)
fwrite(ed_gsea, file.path(tmp_dir, "ed_signature_gsea.tsv"), sep = "\t")
ed_overlap <- rbindlist(ed_overlap, fill = TRUE)
setorder(ed_overlap, target_gene, signature)
fwrite(ed_overlap, file.path(tmp_dir, "ed_top5_overlap.tsv"), sep = "\t")
ed_concordance <- rbindlist(ed_concordance, fill = TRUE)
setorder(ed_concordance, target_gene)
fwrite(ed_concordance, file.path(tmp_dir, "ed_rank_concordance.tsv"), sep = "\t")

main_hypoxia <- main_results[database == "Hallmark" & pathway == "HALLMARK_HYPOXIA", .(
  target_gene, condition, ES, NES, pval, padj, size, leadingEdge
)]
setorder(main_hypoxia, target_gene)
fwrite(main_hypoxia, file.path(tmp_dir, "hypoxia_primary_summary.tsv"), sep = "\t")

report <- list(
  gate = "VP-G06",
  stage = "ordered_enrichment_and_ED_signature_overlap",
  status = "PASS_COMPUTED_PENDING_INDEPENDENT_AUDIT",
  created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  freeze_sha256 = digest(file = freeze_path, algo = "sha256"),
  gene_universe_size = length(gene_universe),
  target_count = nrow(target_specs),
  database_counts = lapply(names(gene_sets), function(x) list(database = x, total = length(gene_sets[[x]]), tested = sum(eligibility$database == x & eligibility$status == "TESTED"))),
  target_consensus_enrichment_rows = nrow(main_results),
  donor_hallmark_rows = nrow(donor_hallmark),
  matched_profile_hallmark_rows = nrow(primary_hallmark),
  ED_signature_sizes_in_frozen_universe = as.list(lengths(ed_sets)),
  random_seed_start = freeze$random_seed,
  random_seed_last_call = seed_counter,
  interpretation_boundary = "Positive NES means concentration toward the high network-perturbation end; negative NES means concentration toward the low end. Neither denotes pathway up/down regulation or ED-expression reversal.",
  downstream_gate_status = "VP-G07_NOT_STARTED"
)
write_json(report, file.path(tmp_dir, "g06_report.json"), pretty = TRUE, auto_unbox = TRUE)
capture.output(sessionInfo(), file = file.path(tmp_dir, "sessionInfo.txt"))

artifact_names <- c(
  "gene_set_eligibility.tsv", "target_consensus_rank_inputs.tsv", "target_consensus_enrichment.tsv", "donor_hallmark_enrichment.tsv",
  "hypoxia_donor_enrichment.tsv", "hypoxia_donor_consistency.tsv", "matched_profile_hallmark_enrichment.tsv",
  "matched_profile_rank_inputs.tsv", "hypoxia_matched_control_enrichment.tsv", "hypoxia_control_calibration.tsv", "ed_full_rank.tsv",
  "ed_signature_sets.tsv", "ed_signature_gsea.tsv", "ed_top5_overlap.tsv", "ed_rank_concordance.tsv",
  "hypoxia_primary_summary.tsv", "g06_report.json", "sessionInfo.txt"
)
manifest <- data.table(
  file = artifact_names,
  bytes = file.info(file.path(tmp_dir, artifact_names))$size,
  sha256 = vapply(file.path(tmp_dir, artifact_names), digest, character(1), algo = "sha256", file = TRUE)
)
fwrite(manifest, file.path(tmp_dir, "artifact_manifest.sha256.tsv"), sep = "\t")
if (!dir.create(dirname(out_dir), recursive = TRUE, showWarnings = FALSE) && !dir.exists(dirname(out_dir))) stop("Cannot create output parent")
if (!file.rename(tmp_dir, out_dir)) stop("Could not seal output directory")
cat(toJSON(report, pretty = TRUE, auto_unbox = TRUE), "\n")
