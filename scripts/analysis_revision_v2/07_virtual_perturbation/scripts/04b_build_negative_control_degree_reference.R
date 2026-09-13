args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 0L) stop("This frozen script takes no arguments.", call. = FALSE)

root <- "C:/Users/22394/Documents/Codex/2026-07-12/acad"
vp <- file.path(root, "revision_v2/07_virtual_perturbation")
rlib <- file.path(vp, "02_env/Rlib")
.libPaths(c(rlib, .libPaths()))
rule_path <- file.path(vp, "01_protocol/VP_G05_negative_control_rule_freeze_v1.json")
out <- file.path(vp, "06_consensus/negative_control_degree_reference")
if (!file.exists(rule_path)) stop("Negative-control rule is not frozen.", call. = FALSE)
if (file.exists(out)) stop("Degree-reference output already exists; refusing overwrite.", call. = FALSE)
if (as.character(utils::packageVersion("scTenifoldKnk")) != "1.0.3") stop("Frozen package version scTenifoldKnk 1.0.3 is required.", call. = FALSE)

sha256_file <- function(path) {
  output <- suppressWarnings(system2("certutil", c("-hashfile", shQuote(normalizePath(path, winslash = "\\", mustWork = TRUE)), "SHA256"), stdout = TRUE, stderr = TRUE))
  candidates <- gsub("\\s", "", output[grepl("^[0-9A-Fa-f ]+$", output)])
  candidates <- candidates[nchar(candidates) == 64L]
  if (length(candidates) != 1L) stop(sprintf("Unable to hash %s", path), call. = FALSE)
  tolower(candidates[[1L]])
}
write_tsv <- function(x, path) utils::write.table(x, file = path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
decile <- function(metric, genes) {
  order_index <- order(metric, genes, na.last = NA)
  ranks <- integer(length(metric)); ranks[order_index] <- seq_along(order_index)
  floor(10 * (ranks - 1L) / length(metric)) + 1L
}

temp <- paste0(out, ".tmp_", Sys.getpid())
dir.create(temp, recursive = TRUE, showWarnings = FALSE)
on_error <- function() {
  if (dir.exists(temp)) file.rename(temp, paste0(out, ".attempt_01_FAIL"))
}
tryCatch({
  suppressPackageStartupMessages(library(Matrix))
  suppressPackageStartupMessages(library(scTenifoldNet))
  definitions <- list(
    list(condition = "normal", donors = c("Normal_1", "Normal_2", "Normal_3")),
    list(condition = "organic_ED_nonDM", donors = c("non-DM_1", "non-DM_2", "non-DM_3"))
  )
  rows <- list(); qc_rows <- list(); index <- 1L; qc_index <- 1L
  for (definition in definitions) for (donor in definition$donors) {
    matrix_dir <- file.path(vp, "04_prepared/fibroblast", definition$condition, donor)
    genes <- utils::read.delim(file.path(matrix_dir, "genes.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
    cells <- utils::read.delim(file.path(matrix_dir, "cells.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
    connection <- gzfile(file.path(matrix_dir, "matrix.mtx.gz"), open = "rt")
    counts <- Matrix::readMM(connection); close(connection)
    if (!inherits(counts, "sparseMatrix")) counts <- Matrix::Matrix(counts, sparse = TRUE)
    if (nrow(counts) != nrow(genes) || ncol(counts) != nrow(cells) || !identical(genes$gene_order, seq_len(nrow(genes)))) stop(sprintf("Frozen matrix schema mismatch: %s", matrix_dir), call. = FALSE)
    rownames(counts) <- genes$gene
    set.seed(2026082701L)
    networks <- scTenifoldNet::makeNetworks(X = counts, nNet = 10L, nCells = 500L, nComp = 3L, scaleScores = TRUE, symmetric = FALSE, q = 0.9, nCores = 4L, label = paste(definition$condition, donor))
    if (length(networks) != 10L) stop(sprintf("Expected 10 networks for %s", donor), call. = FALSE)
    degrees <- matrix(0, nrow = nrow(counts), ncol = 10L, dimnames = list(genes$gene, NULL))
    for (network_index in seq_along(networks)) {
      network <- networks[[network_index]]
      if (!identical(rownames(network), genes$gene) || !identical(colnames(network), genes$gene)) stop(sprintf("Network gene order mismatch: %s net %d", donor, network_index), call. = FALSE)
      degrees[, network_index] <- Matrix::colSums(network != 0)
      qc_rows[[qc_index]] <- data.frame(condition = definition$condition, donor = donor, network_index = network_index, genes = nrow(network), nonzero_edges = length(network@x), stringsAsFactors = FALSE); qc_index <- qc_index + 1L
    }
    mean_degree <- rowMeans(degrees)
    detection <- Matrix::rowSums(counts != 0) / ncol(counts)
    rows[[index]] <- data.frame(gene = genes$gene, condition = definition$condition, donor = donor, fibroblast_cells = ncol(counts), detection_rate = detection, mean_outdegree_raw_10net = mean_degree, detection_decile = decile(detection, genes$gene), outdegree_decile = decile(mean_degree, genes$gene), stringsAsFactors = FALSE); index <- index + 1L
  }
  degree_rows <- do.call(rbind, rows)
  qc <- do.call(rbind, qc_rows)
  write_tsv(degree_rows, file.path(temp, "degree_reference.tsv"))
  write_tsv(qc, file.path(temp, "network_qc.tsv"))
  metadata <- list(gate = "VP-G05", stage = "negative_control_outdegree_reference", status = "PASS_TECHNICAL", created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), rule_sha256 = sha256_file(rule_path), package = "scTenifoldKnk/scTenifoldNet", package_version = as.character(utils::packageVersion("scTenifoldKnk")), definition = "Mean non-zero raw-PCNet regulator-column count across 10 seeded makeNetworks networks; no cross-donor aggregation.")
  writeLines(rjson::toJSON(metadata), file.path(temp, "metadata.json"), useBytes = TRUE)
  writeLines(capture.output(sessionInfo()), file.path(temp, "session_info.txt"), useBytes = TRUE)
  artifact_paths <- list.files(temp, full.names = TRUE, recursive = FALSE)
  manifest <- data.frame(file = basename(artifact_paths), bytes = as.numeric(file.info(artifact_paths)$size), sha256 = vapply(artifact_paths, sha256_file, character(1)), stringsAsFactors = FALSE)
  write_tsv(manifest, file.path(temp, "artifact_manifest.sha256.tsv"))
  if (!file.rename(temp, out)) stop("Atomic degree-reference rename failed.", call. = FALSE)
  cat(sprintf('{"status":"PASS_TECHNICAL","output":"%s"}\n', out))
}, error = function(e) { on_error(); stop(e) })
