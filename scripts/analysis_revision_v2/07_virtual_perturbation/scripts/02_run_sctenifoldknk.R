args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag) {
  index <- match(flag, args)
  if (is.na(index) || index == length(args)) stop(sprintf("Missing value for %s", flag), call. = FALSE)
  args[[index + 1L]]
}

matrix_dir <- normalizePath(value_after("--matrix-dir"), winslash = "/", mustWork = TRUE)
target <- value_after("--target")
seed <- as.integer(value_after("--seed"))
output_attempt <- normalizePath(value_after("--output-attempt"), winslash = "/", mustWork = FALSE)
n_cores <- as.integer(value_after("--n-cores"))
rlib <- Sys.getenv("R_LIBS_USER", unset = "")
if (!nzchar(rlib)) stop("R_LIBS_USER must point to the task-specific R library.", call. = FALSE)
rlib <- normalizePath(rlib, winslash = "/", mustWork = TRUE)
.libPaths(c(rlib, .libPaths()))

if (n_cores != 4L) stop("Frozen VP-G01 protocol requires --n-cores 4.", call. = FALSE)
if (file.exists(output_attempt)) stop("Output attempt already exists; refusing to overwrite retained evidence.", call. = FALSE)
if (as.character(utils::packageVersion("scTenifoldKnk")) != "1.0.3") stop("Frozen protocol requires scTenifoldKnk 1.0.3.", call. = FALSE)

sha256_file <- function(path) {
  if (file.info(path)$size == 0) {
    # certutil returns ERROR_FILE_INVALID for empty files on this Windows host.
    # SHA-256 of an empty byte stream is fixed by the standard.
    return("e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855")
  }
  result <- suppressWarnings(system2("certutil", c("-hashfile", shQuote(normalizePath(path, winslash = "\\", mustWork = TRUE)), "SHA256"), stdout = TRUE, stderr = TRUE))
  candidates <- gsub("\\s", "", result[grepl("^[0-9A-Fa-f ]+$", result)])
  candidates <- candidates[nchar(candidates) == 64L]
  if (length(candidates) != 1L) stop(sprintf("Unable to calculate SHA-256 for %s", path), call. = FALSE)
  tolower(candidates[[1L]])
}

write_tsv <- function(data, path) {
  utils::write.table(data, file = path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
}

temporary_attempt <- paste0(output_attempt, ".tmp_", Sys.getpid())
dir.create(temporary_attempt, recursive = TRUE, showWarnings = FALSE)
stdout_connection <- file.path(temporary_attempt, "stdout.log")
stderr_connection <- file.path(temporary_attempt, "stderr.log")
stdout_con <- file(stdout_connection, open = "wt", encoding = "UTF-8")
stderr_con <- file(stderr_connection, open = "wt", encoding = "UTF-8")
sink(stdout_con, split = TRUE)
sink(stderr_con, type = "message")

restore_sinks <- function() {
  # Exactly one output and one message diversion are created above.  Do not
  # loop on sink.number(): R's reported sink depth can remain nonzero after
  # package-level message handling and previously caused a post-run CPU hold.
  if (sink.number(type = "message") > 0L) sink(type = "message")
  if (sink.number() > 0L) sink()
  if (isOpen(stderr_con)) close(stderr_con)
  if (isOpen(stdout_con)) close(stdout_con)
}

started_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)
started_clock <- proc.time()
result <- tryCatch({
  suppressPackageStartupMessages(library("Matrix", character.only = TRUE))
  suppressPackageStartupMessages(library("scTenifoldKnk", lib.loc = rlib, character.only = TRUE))
  genes <- utils::read.delim(file.path(matrix_dir, "genes.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
  cells <- utils::read.delim(file.path(matrix_dir, "cells.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
  if (!identical(genes$gene_order, seq_len(nrow(genes)))) stop("genes.tsv order is invalid.", call. = FALSE)
  if (!identical(cells$cell_order, seq_len(nrow(cells)))) stop("cells.tsv order is invalid.", call. = FALSE)
  if (sum(genes$gene == target) != 1L) stop(sprintf("Target %s is not uniquely present in genes.tsv.", target), call. = FALSE)
  target_qc <- utils::read.delim(file.path(matrix_dir, "target_detection.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
  target_qc_row <- target_qc[target_qc$gene == target, , drop = FALSE]
  if (nrow(target_qc_row) != 1L || target_qc_row$status[[1L]] != "PASS") stop("Target detection gate is not PASS for this matrix.", call. = FALSE)
  connection <- gzfile(file.path(matrix_dir, "matrix.mtx.gz"), open = "rt")
  gene_by_cell <- Matrix::readMM(connection)
  close(connection)
  if (!inherits(gene_by_cell, "sparseMatrix")) gene_by_cell <- Matrix::Matrix(gene_by_cell, sparse = TRUE)
  if (nrow(gene_by_cell) != nrow(genes) || ncol(gene_by_cell) != nrow(cells)) stop("Matrix dimensions do not match genes.tsv/cells.tsv.", call. = FALSE)
  rownames(gene_by_cell) <- genes$gene
  colnames(gene_by_cell) <- cells$cell_id
  if (any(!is.finite(gene_by_cell@x)) || any(gene_by_cell@x < 0)) stop("Input matrix contains invalid values.", call. = FALSE)
  set.seed(seed)
  output <- scTenifoldKnk::scTenifoldKnk(
    countMatrix = gene_by_cell,
    qc = FALSE,
    gKO = target,
    nc_lambda = 0,
    nc_nNet = 10,
    nc_nCells = 500,
    nc_nComp = 3,
    nc_scaleScores = TRUE,
    nc_symmetric = FALSE,
    nc_q = 0.90,
    td_K = 3,
    td_maxIter = 1000,
    td_maxError = 1e-5,
    td_nDecimal = 3,
    ma_nDim = 2,
    nCores = n_cores
  )
  differential <- as.data.frame(output$diffRegulation, stringsAsFactors = FALSE)
  required_columns <- c("gene", "distance", "Z", "FC", "p.value", "p.adj")
  if (!identical(colnames(differential), required_columns)) stop(sprintf("Unexpected differential-regulation columns: %s", paste(colnames(differential), collapse = ",")), call. = FALSE)
  if (nrow(differential) != nrow(genes) || anyDuplicated(differential$gene) || !setequal(differential$gene, genes$gene)) stop("Differential-regulation gene set does not exactly match input genes.", call. = FALSE)
  numeric_columns <- setdiff(required_columns, "gene")
  if (any(!is.finite(as.matrix(differential[numeric_columns])))) stop("Differential-regulation output contains NA or Inf.", call. = FALSE)
  if (any(differential$p.value < 0 | differential$p.value > 1 | differential$p.adj < 0 | differential$p.adj > 1)) stop("p.value or p.adj lies outside [0,1].", call. = FALSE)
  write_tsv(differential, file.path(temporary_attempt, "differential_regulation.tsv"))
  parameters <- list(
    gate = "VP-G04",
    purpose = "Main frozen virtual-KO; biological interpretation prohibited.",
    matrix_dir = matrix_dir,
    target = target,
    seed = seed,
    n_cores = n_cores,
    package = "scTenifoldKnk",
    package_version = as.character(utils::packageVersion("scTenifoldKnk")),
    input_dimensions_gene_by_cell = c(nrow(gene_by_cell), ncol(gene_by_cell)),
    frozen_parameters = list(qc = FALSE, nc_lambda = 0, nc_nNet = 10, nc_nCells = 500, nc_nComp = 3, nc_scaleScores = TRUE, nc_symmetric = FALSE, nc_q = 0.90, td_K = 3, td_maxIter = 1000, td_maxError = 1e-5, td_nDecimal = 3, ma_nDim = 2, nCores = 4),
    input_sha256 = list(matrix_mtx_gz = sha256_file(file.path(matrix_dir, "matrix.mtx.gz")), genes_tsv = sha256_file(file.path(matrix_dir, "genes.tsv")), cells_tsv = sha256_file(file.path(matrix_dir, "cells.tsv")), provenance_json = sha256_file(file.path(matrix_dir, "provenance.json")))
  )
  writeLines(rjson::toJSON(parameters), file.path(temporary_attempt, "run_parameters.json"), useBytes = TRUE)
  writeLines(capture.output(sessionInfo()), file.path(temporary_attempt, "session_info.txt"), useBytes = TRUE)
  gc_summary <- gc()
  timing <- proc.time() - started_clock
  resource_usage <- data.frame(
    metric = c("elapsed_seconds", "r_gc_max_vcells", "r_gc_max_ncells", "input_genes", "input_cells"),
    value = c(unname(timing[["elapsed"]]), gc_summary["Vcells", "max used"], gc_summary["Ncells", "max used"], nrow(gene_by_cell), ncol(gene_by_cell))
  )
  write_tsv(resource_usage, file.path(temporary_attempt, "resource_usage.tsv"))
  list(success = TRUE, started_at = started_at, finished_at = format(Sys.time(), tz = "UTC", usetz = TRUE), output_rows = nrow(differential))
}, error = function(error) {
  list(success = FALSE, error_type = class(error)[[1L]], error_message = conditionMessage(error), started_at = started_at, failed_at = format(Sys.time(), tz = "UTC", usetz = TRUE))
})
restore_sinks()

if (!isTRUE(result$success)) {
  writeLines(rjson::toJSON(c(list(gate = "VP-G04", status = "FAIL_RETAINED"), result)), file.path(temporary_attempt, "failure.json"), useBytes = TRUE)
  failure_attempt <- paste0(output_attempt, "_FAIL")
  if (file.exists(failure_attempt)) failure_attempt <- paste0(failure_attempt, "_", format(Sys.time(), "%Y%m%d%H%M%S", tz = "UTC"))
  file.rename(temporary_attempt, failure_attempt)
  stop(result$error_message, call. = FALSE)
}

artifact_paths <- list.files(temporary_attempt, full.names = TRUE, recursive = FALSE)
artifact_rows <- data.frame(
  file = basename(artifact_paths),
  bytes = as.numeric(file.info(artifact_paths)$size),
  sha256 = vapply(artifact_paths, sha256_file, character(1)),
  stringsAsFactors = FALSE
)
write_tsv(artifact_rows, file.path(temporary_attempt, "artifact_manifest.sha256.tsv"))
if (!file.rename(temporary_attempt, output_attempt)) stop("Atomic output-attempt rename failed.", call. = FALSE)
