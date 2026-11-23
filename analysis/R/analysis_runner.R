# Standalone decision runner; persistent records never contain expression objects.
.analysis_metadata <- function(object) {
  if (inherits(object, "SingleCellExperiment")) return(S4Vectors::metadata(object))
  if (inherits(object, "Seurat")) return(object@misc)
  stop("Analysis steps require an existing SingleCellExperiment or Seurat object.")
}

.analysis_mark <- function(object, state) {
  if (inherits(object, "SingleCellExperiment")) {
    S4Vectors::metadata(object)$analysis_runner <- state
  } else if (inherits(object, "Seurat")) {
    object@misc$analysis_runner <- state
  } else stop("Analysis steps must return a SingleCellExperiment or Seurat object.")
  object
}

.analysis_preserve_geo <- function(input, candidate) {
  original <- .analysis_metadata(input)$geo
  returned <- .analysis_metadata(candidate)$geo
  cells <- colnames(candidate)
  if (is.null(cells) || anyNA(cells) || any(!nzchar(cells)) || anyDuplicated(cells) ||
      anyDuplicated(colnames(input)) || !all(cells %in% colnames(input))) {
    stop("Step output must retain unique input cell identifiers (subsets are allowed).")
  }
  if (!is.null(original)) {
    if (!is.null(returned) && !identical(original, returned)) stop("Step changed original GEO provenance.")
    if (inherits(candidate, "SingleCellExperiment")) S4Vectors::metadata(candidate)$geo <- original else
      candidate@misc$geo <- original
  }
  before <- if (inherits(input, "SingleCellExperiment")) as.data.frame(SummarizedExperiment::colData(input)) else input[[]]
  after <- if (inherits(candidate, "SingleCellExperiment")) as.data.frame(SummarizedExperiment::colData(candidate)) else candidate[[]]
  for (field in names(before)[startsWith(names(before), "geo_")]) {
    values <- before[[field]][match(cells, colnames(input))]
    if (field %in% names(after) && !identical(unname(after[[field]]), unname(values))) {
      stop("Step changed original GEO cell metadata: ", field)
    }
    if (inherits(candidate, "SingleCellExperiment")) SummarizedExperiment::colData(candidate)[[field]] <- values else
      candidate[[field]] <- values
  }
  candidate
}

.analysis_data_only <- function(value) {
  if (isS4(value) || is.environment(value) || is.function(value) ||
      (!is.null(dim(value)) && !is.data.frame(value))) {
    stop("Parameters and diagnostics must be compact data, not objects, matrices, functions, or environments.")
  }
  if (is.list(value)) {
    for (part in value) .analysis_data_only(part)
  } else if (!is.null(value) && !is.atomic(value)) {
    stop("Parameters and diagnostics contain unsupported values.")
  }
  invisible(TRUE)
}

.analysis_fingerprint <- function(object) digest::digest(object, algo = "sha256")

.analysis_artifacts <- function(paths = character()) {
  if (is.null(paths)) paths <- character()
  if (!is.character(paths) || anyNA(paths) || any(!file.exists(paths)) || any(dir.exists(paths))) {
    stop("Diagnostic evidence must reference existing files.")
  }
  data.frame(path = normalizePath(paths, mustWork = TRUE),
             checksum = unname(tools::md5sum(paths)), stringsAsFactors = FALSE)
}

.analysis_write_once <- function(record, path) {
  if (file.exists(path)) stop("An analysis record already exists and will not be overwritten.")
  temporary <- tempfile("record-", tmpdir = dirname(path))
  on.exit(unlink(temporary))
  saveRDS(record, temporary)
  if (!file.rename(temporary, path)) stop("Could not save analysis record.")
}

.analysis_seeded <- function(seed, execute) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  previous <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(kind))
    if (had_seed) assign(".Random.seed", previous, envir = .GlobalEnv) else
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  })
  set.seed(seed)
  execute()
}

.analysis_code <- function(execute, files) {
  files <- lapply(files, function(path) {
    if (!file.exists(path)) stop("Declared source file does not exist: ", path)
    list(path = normalizePath(path), checksum = unname(tools::md5sum(path)),
         text = readLines(path, warn = FALSE))
  })
  git <- Sys.which("git")
  revision <- if (nzchar(git)) suppressWarnings(system2(git, c("rev-parse", "HEAD"),
                                                       stdout = TRUE, stderr = FALSE)) else character()
  status <- if (nzchar(git)) suppressWarnings(system2(git, c("status", "--porcelain"),
                                                     stdout = TRUE, stderr = FALSE)) else character()
  list(function_body = paste(deparse(execute), collapse = "\n"), files = files,
       git_revision = if (length(revision) && is.null(attr(revision, "status"))) revision[1] else NA_character_,
       git_dirty = if (is.null(attr(status, "status"))) length(status) > 0L else NA,
       R_version = R.version.string, rng_kind = RNGkind())
}

.analysis_run_path <- function(directory, run_id) {
  if (!is.character(run_id) || length(run_id) != 1L || is.na(run_id) ||
      !grepl("^run-[A-Za-z0-9-]+$", run_id)) stop("Invalid analysis run identifier.")
  file.path(directory, run_id)
}

.analysis_decisions <- function(directory) {
  paths <- list.files(directory, pattern = "^decision\\.rds$", recursive = TRUE, full.names = TRUE)
  lapply(paths, readRDS)
}

.analysis_head <- function(directory) {
  accepted <- Filter(function(x) x$decision == "accept", .analysis_decisions(directory))
  if (!length(accepted)) return(NULL)
  accepted[[which.max(vapply(accepted, function(x) x$sequence, integer(1)))]]$run_id
}

.analysis_chain <- function(directory, head = .analysis_head(directory)) {
  chain <- character()
  while (!is.null(head)) {
    if (head %in% chain) stop("Analysis history contains a cycle.")
    chain <- c(chain, head)
    proposal <- readRDS(file.path(.analysis_run_path(directory, head), "proposal.rds"))
    head <- proposal$parent_run
  }
  chain
}

.analysis_locked <- function(directory, execute) {
  lock <- file.path(directory, ".analysis-lock")
  if (!dir.create(lock, showWarnings = FALSE)) stop("Analysis log is locked; inspect the active writer before retrying.")
  on.exit(unlink(lock, recursive = TRUE))
  execute()
}

decide_analysis_step <- function(proposal, decision, reason) {
  if (!is.character(decision) || length(decision) != 1L || is.na(decision) ||
      !decision %in% c("accept", "reject")) stop("Decision must be explicitly accept or reject.")
  if (!is.character(reason) || length(reason) != 1L || is.na(reason) || !nzchar(trimws(reason))) {
    stop("A nonempty decision reason is required.")
  }
  directory <- normalizePath(proposal$directory, mustWork = TRUE)
  run <- .analysis_run_path(directory, proposal$run_id)
  .analysis_locked(directory, function() {
    if (file.exists(file.path(run, "decision.rds"))) stop("This proposal already has an immutable decision.")
    if (!file.exists(file.path(run, "result.rds")) || file.exists(file.path(run, "failure.rds"))) {
      stop("Only a successful proposal can be reviewed.")
    }
    record <- readRDS(file.path(run, "proposal.rds"))
    result <- readRDS(file.path(run, "result.rds"))
    object <- NULL
    if (decision == "accept") {
      if (!identical(record$head_at_proposal, .analysis_head(directory))) {
        stop("Analysis advanced since this proposal; rerun before acceptance.")
      }
      if (!identical(.analysis_fingerprint(proposal$object), result$object_fingerprint) ||
          !identical(proposal$diagnostics, result$diagnostics) ||
          !identical(proposal$artifacts, result$artifacts)) stop("Proposal changed after execution.")
      if (any(!file.exists(result$artifacts$path)) ||
          !identical(unname(tools::md5sum(result$artifacts$path)), result$artifacts$checksum)) {
        stop("Diagnostic evidence files changed after execution.")
      }
      if (!identical(unname(tools::md5sum(record$input$path)), record$input$checksum)) {
        stop("Input checkpoint changed after execution.")
      }
      object <- .analysis_mark(proposal$object, list(run_id = proposal$run_id,
                                                    directory = directory, status = "accepted"))
    }
    decisions <- .analysis_decisions(directory)
    sequence <- if (length(decisions)) max(vapply(decisions, function(x) x$sequence, integer(1))) + 1L else 1L
    .analysis_write_once(list(run_id = proposal$run_id, decision = decision, reason = reason,
                              sequence = sequence, decided_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
                              object_fingerprint = if (is.null(object)) NULL else .analysis_fingerprint(object)),
                         file.path(run, "decision.rds"))
    object
  })
}

read_analysis_history <- function(directory) {
  if (!file.exists(file.path(directory, "analysis.rds"))) stop("No analysis log exists in this directory.")
  runs <- list.dirs(directory, recursive = FALSE, full.names = TRUE)
  runs <- runs[grepl("^run-", basename(runs))]
  history <- data.frame(run_id = character(), step = character(), status = character(),
                        parent_run = character(), seed = integer(), reason = character(),
                        decision = character(), current = logical(),
                        stringsAsFactors = FALSE)
  history$parameters <- I(list())
  history$diagnostics <- I(list())
  history$code <- I(list())
  history$input <- I(list())
  history$artifacts <- I(list())
  head <- .analysis_head(directory)
  chain <- .analysis_chain(directory, head)
  for (run in sort(runs)) {
    path <- file.path(run, "proposal.rds")
    if (!file.exists(path)) next
    record <- readRDS(path)
    result <- if (file.exists(file.path(run, "result.rds"))) readRDS(file.path(run, "result.rds")) else NULL
    failed <- file.exists(file.path(run, "failure.rds"))
    decision <- if (file.exists(file.path(run, "decision.rds"))) readRDS(file.path(run, "decision.rds")) else NULL
    status <- if (failed) "failed" else if (is.null(result)) "running" else "proposed"
    if (!is.null(decision)) {
      status <- if (decision$decision == "reject") "rejected" else
        if (basename(run) %in% chain) "accepted" else "superseded"
    }
    row <- data.frame(run_id = basename(run), step = record$step,
                      status = status,
                      parent_run = if (is.null(record$parent_run)) NA_character_ else record$parent_run,
                      seed = record$seed,
                      reason = if (!is.null(decision)) decision$reason else
                        if (failed) readRDS(file.path(run, "failure.rds"))$message else "",
                      decision = if (is.null(decision)) "pending" else decision$decision,
                      current = identical(basename(run), head), stringsAsFactors = FALSE)
    row$parameters <- I(list(record$parameters))
    row$diagnostics <- I(list(if (is.null(result)) list() else result$diagnostics))
    code <- record$code
    if (!is.null(result)) code$packages_after <- result$packages
    row$code <- I(list(code))
    row$input <- I(list(record$input))
    row$artifacts <- I(list(if (is.null(result)) .analysis_artifacts() else result$artifacts))
    history <- rbind(history, row)
  }
  rownames(history) <- NULL
  history
}

propose_analysis_step <- function(object, directory, step, execute, parameters = list(),
                                   seed = 1L, input_path = NULL, parent_run = NULL,
                                   source_files = character()) {
  if (!requireNamespace("digest", quietly = TRUE)) stop("The runner requires an existing digest installation.")
  .analysis_metadata(object)
  if (!is.function(execute) || !is.character(step) || length(step) != 1L || is.na(step) || !nzchar(step)) {
    stop("Supply a step name and an R function accepting object and parameters.")
  }
  if (!is.list(parameters) || (length(parameters) &&
      (is.null(names(parameters)) || any(!nzchar(names(parameters))) || anyDuplicated(names(parameters))))) {
    stop("Parameters must be a named list with unique names.")
  }
  .analysis_data_only(parameters)
  if (length(seed) != 1L || !is.numeric(seed) || is.na(seed) || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != floor(seed)) stop("Seed must be a nonnegative integer.")
  fingerprint <- .analysis_fingerprint(object)
  manifest <- file.path(directory, "analysis.rds")
  if (!file.exists(manifest)) {
    if (!is.null(parent_run)) stop("A parent result must be explicitly accepted before downstream execution.")
    if (is.null(input_path) || !file.exists(input_path)) stop("Supply the existing input checkpoint path.")
    input <- list(path = normalizePath(input_path), checksum = unname(tools::md5sum(input_path)))
    if (dir.exists(directory) && length(list.files(directory, all.files = TRUE, no.. = TRUE))) {
      stop("Analysis directory is not empty; use a dedicated directory.")
    }
    if (!dir.exists(directory) && !dir.create(directory, recursive = TRUE)) stop("Cannot create analysis directory.")
    .analysis_locked(directory, function() {
      .analysis_write_once(list(schema_version = 1L, input = input, object_fingerprint = fingerprint), manifest)
    })
  }
  directory <- normalizePath(directory)
  baseline <- readRDS(manifest)
  if (!identical(baseline$schema_version, 1L)) stop("Unsupported analysis log schema.")
  input <- baseline$input
  if (!is.null(input_path) && !identical(normalizePath(input_path, mustWork = TRUE), input$path)) {
    stop("Input differs from this analysis baseline.")
  }
  if (!file.exists(input$path) || !identical(unname(tools::md5sum(input$path)), input$checksum)) {
    stop("Input checkpoint changed from the analysis baseline.")
  }
  head <- .analysis_head(directory)
  if (is.null(parent_run)) {
    if (!identical(baseline$object_fingerprint, fingerprint)) stop("Input differs from this analysis baseline.")
  } else {
    parent <- .analysis_run_path(directory, parent_run)
    if (!file.exists(file.path(parent, "decision.rds")) || !parent_run %in% .analysis_chain(directory, head)) {
      stop("Parent must be accepted and on the current analysis lineage.")
    }
    decision <- readRDS(file.path(parent, "decision.rds"))
    if (decision$decision != "accept" || !identical(decision$object_fingerprint, fingerprint)) {
      stop("Input does not match the accepted parent result.")
    }
  }
  run <- tempfile("run-", tmpdir = directory)
  if (!dir.create(run)) stop("Cannot create step record directory.")
  run_id <- basename(run)
  record <- list(run_id = run_id, step = step, parameters = parameters, seed = as.integer(seed),
                 parent_run = parent_run, head_at_proposal = head,
                 input = input, code = .analysis_code(execute, source_files),
                 created_at = format(Sys.time(), tz = "UTC", usetz = TRUE))
  .analysis_write_once(record, file.path(run, "proposal.rds"))
  outcome <- tryCatch({
    candidate <- .analysis_seeded(seed, function() execute(object, parameters))
    if (inherits(candidate, "SingleCellExperiment") || inherits(candidate, "Seurat")) candidate <- list(object = candidate)
    if (!is.list(candidate) || is.null(candidate$object)) stop("Step must return an object or a list containing object.")
    if (any(!names(candidate) %in% c("object", "diagnostics", "artifacts")) ||
        is.null(names(candidate)) || anyDuplicated(names(candidate))) stop("Unsupported step return fields.")
    if (is.null(candidate$diagnostics)) candidate$diagnostics <- list()
    .analysis_data_only(candidate$diagnostics)
    artifacts <- .analysis_artifacts(candidate$artifacts)
    candidate$object <- .analysis_preserve_geo(object, candidate$object)
    candidate$object <- .analysis_mark(candidate$object, list(run_id = run_id, directory = directory, status = "proposed"))
    result <- list(diagnostics = candidate$diagnostics, artifacts = artifacts,
                   object_fingerprint = .analysis_fingerprint(candidate$object),
                   object_class = class(candidate$object), dimensions = dim(candidate$object),
                   packages = vapply(sort(loadedNamespaces()), function(x) as.character(utils::packageVersion(x)), character(1)))
    .analysis_write_once(result, file.path(run, "result.rds"))
    list(object = candidate$object, diagnostics = candidate$diagnostics, artifacts = artifacts,
         run_id = run_id, directory = directory)
  }, error = function(e) {
    .analysis_write_once(list(message = conditionMessage(e)), file.path(run, "failure.rds"))
    stop("Analysis step failed: ", conditionMessage(e), call. = FALSE)
  })
  outcome
}
