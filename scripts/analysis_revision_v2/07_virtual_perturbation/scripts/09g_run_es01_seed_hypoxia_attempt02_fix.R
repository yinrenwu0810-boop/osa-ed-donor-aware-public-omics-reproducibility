args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript 09g_run_es01_seed_hypoxia_attempt02_fix.R <vp_root>", call. = FALSE)
vp_root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
attempt <- file.path(vp_root, "10_exploratory_hypoxia_sensitivity", "attempt_20260910_02")
source_script <- file.path(vp_root, "scripts", "09f_run_es01_seed_hypoxia.R")
failure_path <- file.path(attempt, "03_seed_enrichment_ATTEMPT01_IMPLEMENTATION_FAIL_RETAINED.json")
suppressPackageStartupMessages({ library(jsonlite); library(digest) })
if (!file.exists(failure_path)) {
  write_json(list(
    gate = "VP-ES01-03", stage = "single_seed_full_hallmark_enrichment", status = "IMPLEMENTATION_FAIL_RETAINED",
    created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    failed_script = "scripts/09f_run_es01_seed_hypoxia.R", failed_script_sha256 = digest(file = source_script, algo = "sha256"),
    failure = "data.table::setorder does not accept the inline -abs(NES) expression after the first fgsea result; no 03_* formal output was written.",
    correction = "The append-only runner 09g replaces only that sort expression with an explicit temporary numeric sort column, retaining the full frozen input, fgsea, FDR, and output rules.",
    downstream_action = "Retry VP-ES01-03 only; do not enter VP-ES01-04 until its formal outputs pass structural checks."
  ), failure_path, pretty = TRUE, auto_unbox = TRUE)
}
src <- paste(readLines(source_script, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
needle <- "  setorder(result, padj, -abs(NES), pathway, na.last = TRUE)"
replacement <- paste(c(
  "  result[, abs_NES_sort := -abs(NES)]",
  "  setorder(result, padj, abs_NES_sort, pathway, na.last = TRUE)",
  "  result[, abs_NES_sort := NULL]"
), collapse = "\n")
if (!grepl(needle, src, fixed = TRUE)) stop("Retained source runner no longer has the expected sortable failure line.", call. = FALSE)
src <- sub(needle, replacement, src, fixed = TRUE)
eval(parse(text = src), envir = new.env(parent = globalenv()))
