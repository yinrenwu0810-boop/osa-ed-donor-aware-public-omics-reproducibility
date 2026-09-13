args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag) {
  index <- match(flag, args)
  if (is.na(index) || index == length(args)) stop(sprintf("Missing value for %s", flag), call. = FALSE)
  args[[index + 1L]]
}

rlib <- normalizePath(value_after("--rlib"), winslash = "/", mustWork = FALSE)
out_dir <- normalizePath(value_after("--out-dir"), winslash = "/", mustWork = FALSE)
dir.create(rlib, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(rlib, .libPaths()))

stdout_log <- file.path(out_dir, "install_stdout.log")
stderr_log <- file.path(out_dir, "install_stderr.log")
stdout_connection <- file(stdout_log, open = "wt", encoding = "UTF-8")
stderr_connection <- file(stderr_log, open = "wt", encoding = "UTF-8")
sink(stdout_connection, split = TRUE)
sink(stderr_connection, type = "message")
on.exit({ sink(type = "message"); sink(); close(stderr_connection); close(stdout_connection) }, add = TRUE)

required_version <- "1.0.3"
installed_version <- tryCatch(as.character(utils::packageVersion("scTenifoldKnk")), error = function(e) NA_character_)
if (is.na(installed_version) || installed_version != required_version) {
  utils::install.packages(
    "scTenifoldKnk",
    lib = rlib,
    repos = "https://cloud.r-project.org",
    dependencies = NA,
    quiet = FALSE
  )
}
installed_version <- as.character(utils::packageVersion("scTenifoldKnk"))
if (installed_version != required_version) {
  stop(sprintf("scTenifoldKnk version drift: installed %s, required %s", installed_version, required_version), call. = FALSE)
}
suppressPackageStartupMessages(library("scTenifoldKnk", lib.loc = rlib, character.only = TRUE))
writeLines(installed_version, file.path(out_dir, "scTenifoldKnk_version.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"), useBytes = TRUE)
writeLines(normalizePath(.libPaths(), winslash = "/", mustWork = FALSE), file.path(out_dir, "library_paths.txt"), useBytes = TRUE)
cat(sprintf("INSTALL_PASS scTenifoldKnk=%s R=%s\n", installed_version, R.version.string))
