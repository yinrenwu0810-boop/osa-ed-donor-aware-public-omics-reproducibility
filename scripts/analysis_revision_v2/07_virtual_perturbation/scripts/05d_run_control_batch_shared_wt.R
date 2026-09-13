args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag) {
  index <- match(flag, args)
  if (is.na(index) || index == length(args)) stop(sprintf("Missing value for %s", flag), call. = FALSE)
  args[[index + 1L]]
}

matrix_dir <- normalizePath(value_after("--matrix-dir"), winslash = "/", mustWork = TRUE)
controls_file <- normalizePath(value_after("--controls-file"), winslash = "/", mustWork = TRUE)
benchmark_file <- normalizePath(value_after("--benchmark-file"), winslash = "/", mustWork = TRUE)
output_batch <- normalizePath(value_after("--output-batch"), winslash = "/", mustWork = FALSE)
condition <- value_after("--condition")
donor <- value_after("--donor")
benchmark_gene <- value_after("--benchmark-gene")
seed <- as.integer(value_after("--seed"))
n_cores <- as.integer(value_after("--n-cores"))
freeze_sha256 <- value_after("--freeze-sha256")
if (seed != 2026082701L || n_cores != 4L) stop("Frozen v2 control execution requires seed 2026082701 and nCores=4.", call. = FALSE)

vp <- "C:/Users/22394/Documents/Codex/2026-07-12/acad/revision_v2/07_virtual_perturbation"
rlib <- file.path(vp, "02_env/Rlib")
.libPaths(c(rlib, .libPaths()))
if (as.character(utils::packageVersion("scTenifoldKnk")) != "1.0.3") stop("scTenifoldKnk 1.0.3 is required.", call. = FALSE)

sha256_file <- function(path) {
  if (file.info(path)$size == 0) return("e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855")
  result <- suppressWarnings(system2("certutil", c("-hashfile", shQuote(normalizePath(path, winslash = "\\", mustWork = TRUE)), "SHA256"), stdout = TRUE, stderr = TRUE))
  candidates <- gsub("\\s", "", result[grepl("^[0-9A-Fa-f ]+$", result)])
  candidates <- candidates[nchar(candidates) == 64L]
  if (length(candidates) != 1L) stop(sprintf("Unable to hash %s", path), call. = FALSE)
  tolower(candidates[[1L]])
}
write_tsv <- function(x, path) utils::write.table(x, file = path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
write_json <- function(x, path) writeLines(rjson::toJSON(x), path, useBytes = TRUE)
compute_dr <- function(WT, gene, rng_state) {
  assign(".Random.seed", rng_state, envir = .GlobalEnv)
  KO <- WT; KO[gene, ] <- 0
  MA <- scTenifoldNet::manifoldAlignment(WT, KO, d = 2L, nCores = n_cores)
  getFromNamespace("dRegulation", "scTenifoldKnk")(MA, gene)
}

suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(scTenifoldNet))
suppressPackageStartupMessages(library(scTenifoldKnk))
controls <- utils::read.delim(controls_file, check.names = FALSE, stringsAsFactors = FALSE)
if (!identical(names(controls), c("control_gene", "matched_targets")) || anyDuplicated(controls$control_gene)) stop("Invalid batch control list.", call. = FALSE)
genes <- utils::read.delim(file.path(matrix_dir, "genes.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
cells <- utils::read.delim(file.path(matrix_dir, "cells.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
connection <- gzfile(file.path(matrix_dir, "matrix.mtx.gz"), open = "rt")
counts <- Matrix::readMM(connection); close(connection)
if (!inherits(counts, "sparseMatrix")) counts <- Matrix::Matrix(counts, sparse = TRUE)
rownames(counts) <- genes$gene; colnames(counts) <- cells$cell_id
if (any(!controls$control_gene %in% genes$gene) || !benchmark_gene %in% genes$gene) stop("Control or benchmark gene absent from frozen matrix.", call. = FALSE)

dir.create(output_batch, recursive = TRUE, showWarnings = FALSE)
shared_path <- file.path(output_batch, "shared_wt_reference.rds")
batch_parameters_path <- file.path(output_batch, "batch_parameters.json")
benchmark_path <- file.path(output_batch, "direct_equivalence_benchmark.json")
if (!file.exists(shared_path)) {
  set.seed(seed)
  started <- proc.time()
  networks <- scTenifoldNet::makeNetworks(X = counts, q = 0.9, nNet = 10L, nCells = 500L, scaleScores = TRUE, symmetric = FALSE, nComp = 3L, nCores = n_cores, label = paste(condition, donor, "control WT"))
  tensor <- scTenifoldNet::tensorDecomposition(xList = networks, K = 3L, maxError = 1e-5, maxIter = 1000L, nDecimal = 3L)
  WT <- as.matrix(tensor$X); diag(WT) <- 0; WT <- t(WT)
  rng_state <- get(".Random.seed", envir = .GlobalEnv)
  saveRDS(list(WT = WT, rng_state = rng_state), shared_path, compress = TRUE)
  parameters <- list(gate = "VP-G05", stage = "shared_WT_for_matched_control_KO", condition = condition, donor = donor, seed = seed, n_cores = n_cores, controls_file = controls_file, controls_file_sha256 = sha256_file(controls_file), freeze_sha256 = freeze_sha256, model_parameters = list(qc = FALSE, nc_lambda = 0, nc_nNet = 10, nc_nCells = 500, nc_nComp = 3, nc_scaleScores = TRUE, nc_symmetric = FALSE, nc_q = 0.9, td_K = 3, td_maxIter = 1000, td_maxError = 1e-5, td_nDecimal = 3, ma_nDim = 2, nCores = 4), input_sha256 = list(matrix_mtx_gz = sha256_file(file.path(matrix_dir, "matrix.mtx.gz")), genes_tsv = sha256_file(file.path(matrix_dir, "genes.tsv")), cells_tsv = sha256_file(file.path(matrix_dir, "cells.tsv"))), elapsed_seconds = unname((proc.time() - started)[["elapsed"]]))
  write_json(parameters, batch_parameters_path)
} else {
  shared <- readRDS(shared_path); WT <- shared$WT; rng_state <- shared$rng_state
}
if (!exists("WT")) { shared <- readRDS(shared_path); WT <- shared$WT; rng_state <- shared$rng_state }

if (!file.exists(benchmark_path)) {
  observed <- compute_dr(WT, benchmark_gene, rng_state)
  expected <- utils::read.delim(benchmark_file, check.names = FALSE, stringsAsFactors = FALSE)
  observed <- observed[match(expected$gene, observed$gene), names(expected)]
  numeric_fields <- setdiff(names(expected), "gene")
  max_abs <- vapply(numeric_fields, function(field) max(abs(as.numeric(observed[[field]]) - as.numeric(expected[[field]]))), numeric(1))
  benchmark <- list(gate = "VP-G05", status = if (all(max_abs <= 1e-7)) "PASS_NUMERICAL_EQUIVALENCE" else "FAIL_RETAINED", benchmark_gene = benchmark_gene, direct_result = benchmark_file, direct_result_sha256 = sha256_file(benchmark_file), shared_wt_sha256 = sha256_file(shared_path), tolerance_absolute = 1e-7, maximum_absolute_difference = as.list(max_abs))
  write_json(benchmark, benchmark_path)
  if (benchmark$status != "PASS_NUMERICAL_EQUIVALENCE") stop("Shared-WT implementation failed direct-run numerical equivalence.", call. = FALSE)
}
benchmark <- rjson::fromJSON(file = benchmark_path)
if (benchmark$status != "PASS_NUMERICAL_EQUIVALENCE") stop("Batch benchmark is not PASS.", call. = FALSE)

results_root <- file.path(output_batch, "controls")
dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
progress_path <- file.path(output_batch, "progress.tsv")
for (index in seq_len(nrow(controls))) {
  control <- controls$control_gene[[index]]
  final <- file.path(results_root, control, "attempt_01")
  if (file.exists(final)) next
  temp <- paste0(final, ".tmp_", Sys.getpid())
  dir.create(temp, recursive = TRUE, showWarnings = FALSE)
  started_at <- format(Sys.time(), tz = "UTC", usetz = TRUE); started <- proc.time()
  output <- tryCatch({
    differential <- compute_dr(WT, control, rng_state)
    write_tsv(differential, file.path(temp, "differential_regulation.tsv"))
    parameters <- list(gate = "VP-G05", purpose = "Relaxed-v2 matched-control virtual KO; not a biologically inert negative control.", condition = condition, donor = donor, control_gene = control, matched_targets = controls$matched_targets[[index]], seed = seed, n_cores = n_cores, shared_wt_reference = shared_path, shared_wt_sha256 = sha256_file(shared_path), direct_equivalence_benchmark_sha256 = sha256_file(benchmark_path), freeze_sha256 = freeze_sha256, started_at_utc = started_at, finished_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), elapsed_seconds = unname((proc.time() - started)[["elapsed"]]))
    write_json(parameters, file.path(temp, "run_parameters.json"))
    writeLines(capture.output(sessionInfo()), file.path(temp, "session_info.txt"), useBytes = TRUE)
    writeLines(character(0), file.path(temp, "stdout.log"), useBytes = TRUE)
    writeLines(character(0), file.path(temp, "stderr.log"), useBytes = TRUE)
    artifacts <- list.files(temp, full.names = TRUE, recursive = FALSE)
    manifest <- data.frame(file = basename(artifacts), bytes = as.numeric(file.info(artifacts)$size), sha256 = vapply(artifacts, sha256_file, character(1)), stringsAsFactors = FALSE)
    write_tsv(manifest, file.path(temp, "artifact_manifest.sha256.tsv"))
    dir.create(dirname(final), recursive = TRUE, showWarnings = FALSE)
    if (!file.rename(temp, final)) stop("Atomic control-result rename failed.", call. = FALSE)
    data.frame(control_gene = control, matched_targets = controls$matched_targets[[index]], status = "PASS_TECHNICAL", finished_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), output_attempt = final, stringsAsFactors = FALSE)
  }, error = function(e) {
    failure <- paste0(final, "_FAIL")
    if (dir.exists(temp)) file.rename(temp, failure)
    data.frame(control_gene = control, matched_targets = controls$matched_targets[[index]], status = "FAIL_RETAINED", finished_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), output_attempt = failure, stringsAsFactors = FALSE)
  })
  write.table(output, progress_path, sep = "\t", row.names = FALSE, col.names = !file.exists(progress_path), quote = FALSE, append = file.exists(progress_path))
  if (output$status != "PASS_TECHNICAL") stop(sprintf("Control %s failed; retained and stopping batch.", control), call. = FALSE)
}
complete <- list(gate = "VP-G05", stage = "matched_control_KO_batch", status = "PASS_TECHNICAL", condition = condition, donor = donor, seed = seed, controls_expected = nrow(controls), controls_completed = length(list.files(results_root, pattern = "differential_regulation.tsv", recursive = TRUE)), shared_wt_sha256 = sha256_file(shared_path), benchmark_sha256 = sha256_file(benchmark_path), completed_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
write_json(complete, file.path(output_batch, "batch_complete.json"))
cat(rjson::toJSON(complete), "\n")
