#!/usr/bin/env Rscript

main <- function(args) {
  usage <- paste(
    "Usage: Rscript --vanilla scripts/register_gse176031.R",
    "  --data-root DIRECTORY --output-dir EMPTY_DIRECTORY",
    "  [--series-matrix LOCAL_FILE_OR_HTTPS_URL]",
    "  [--confirmed-genome-build hg19] [--r-artifact EXISTING_RDS]",
    "Builds local metadata/header records only; does not load count matrices or run analysis.",
    sep = "\n"
  )
  if (identical(args, "--help")) {
    cat(usage, "\n")
    return(invisible(NULL))
  }
  options <- list("series-matrix" = paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176031/matrix/",
    "GSE176031_series_matrix.txt.gz"), "r-artifact" = character())
  allowed <- c("data-root", "output-dir", "series-matrix", "confirmed-genome-build", "r-artifact")
  while (length(args)) {
    key <- sub("^--", "", args[1])
    if (!startsWith(args[1], "--") || !key %in% allowed || length(args) < 2L ||
        startsWith(args[2], "--") || !nzchar(args[2])) stop(usage)
    if (key == "r-artifact") options[[key]] <- c(options[[key]], args[2]) else options[[key]] <- args[2]
    args <- args[-c(1L, 2L)]
  }
  if (is.null(options[["data-root"]]) || is.null(options[["output-dir"]])) stop(usage)
  output <- options[["output-dir"]]
  if (file.exists(output) && (!dir.exists(output) ||
      length(list.files(output, all.files = TRUE, no.. = TRUE)))) {
    stop("Output directory is not empty; existing outputs will not be overwritten.")
  }
  script <- sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])
  root <- dirname(dirname(normalizePath(script)))
  source(file.path(root, "analysis", "R", "gse176031_registry.R"))
  registry <- build_gse176031_registry(
    options[["data-root"]], options[["series-matrix"]],
    r_artifacts = options[["r-artifact"]],
    confirmed_genome_build = options[["confirmed-genome-build"]]
  )
  if (!dir.exists(output) && !dir.create(output, recursive = TRUE)) stop("Could not create output directory.")
  tables <- list(samples = registry$samples, source_metadata = registry$metadata,
                 inputs = registry$inputs, issues = registry$issues, artifacts = registry$artifacts)
  for (name in names(tables)) {
    utils::write.csv(tables[[name]], file.path(output, paste0(name, ".csv")), row.names = FALSE, na = "")
  }
  saveRDS(registry, file.path(output, "registry.rds"))
  cat("Registered", nrow(registry$samples), "GEO sample records and", nrow(registry$inputs), "DGE files.\n")
  cat("Unambiguous available inputs:", sum(!is.na(registry$samples$dge_file)), "\n")
  cat("Errors:", sum(registry$issues$severity == "error"),
      "Warnings:", sum(registry$issues$severity == "warning"), "\n")
  cat("Local outputs:", normalizePath(output), "\n")
  if (any(registry$issues$severity == "error")) stop("Registry contains input errors; inspect issues.csv before reuse.")
}

tryCatch(main(commandArgs(TRUE)), error = function(e) {
  message(conditionMessage(e))
  quit(status = 1L)
})
