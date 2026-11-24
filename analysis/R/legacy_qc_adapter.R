# Load the unchanged app helpers into a private environment. No package installs.
make_legacy_qc_steps <- function(app_root = "site_data") {
  app_root <- normalizePath(app_root, mustWork = TRUE)
  for (package in c("dplyr", "SingleCellExperiment")) {
    if (!requireNamespace(package, quietly = TRUE)) stop("Existing package required: ", package)
  }
  legacy <- new.env(parent = baseenv())
  legacy[["%>%"]] <- dplyr::`%>%`
  legacy$pull <- dplyr::pull
  # createMetadata() loads a relative annotation file. Redirect only that load,
  # without changing the caller's working directory or editing the original.
  legacy$load <- function(file, envir = parent.frame(), ...) {
    if (!file %in% c("annotations_hsap.Rda", "annotations_mmus.Rda")) stop("Unexpected annotation file.")
    base::load(file.path(app_root, file), envir = envir, ...)
  }
  source_files <- file.path(app_root, "R", c("createMetadata.R", "filter_cells.R", "createQCPlots.R", "annotation_reference_cache.R", "gene_annotations.R"))
  for (path in source_files) sys.source(path, envir = legacy)

  metadata_step <- function(object, parameters) {
    if (!inherits(object, "SingleCellExperiment")) stop("Existing QC metadata requires an SCE object.")
    if (!parameters$organism %in% c("Homo sapiens", "Mus musculus")) stop("Specify human or mouse organism.")
    sample_column <- parameters$sample_column
    if (is.null(sample_column) || !sample_column %in% names(SummarizedExperiment::colData(object))) {
      stop("Specify an existing sample_column; sample/capture identities are not inferred.")
    }
    groups <- as.character(SummarizedExperiment::colData(object)[[sample_column]])
    if (anyNA(groups) || any(!nzchar(groups)) || !length(groups)) stop("Sample identifiers must be nonempty.")
    counts <- SummarizedExperiment::assay(object, "counts")
    indices <- split(seq_len(ncol(object)), groups)
    inputs <- lapply(indices, function(idx) counts[, idx, drop = FALSE])
    message("Calling existing createMetadata() for ", length(inputs), " sample(s).")
    calculated <- legacy$createMetadata(inputs, parameters$organism, TRUE)
    for (name in names(indices)) {
      idx <- indices[[name]]
      result <- calculated[[name]]
      if (!identical(result$counts, inputs[[name]])) {
        stop("Existing createMetadata changed features/counts; review gene mapping before proceeding.")
      }
      fields <- result[setdiff(names(result), "counts")]
      for (field in names(fields)) {
        if (field %in% c("cells","sample") && field %in% names(SummarizedExperiment::colData(object)) &&
            !identical(unname(as.character(object[[field]][idx])),unname(as.character(fields[[field]])))) {
          stop("Requested QC grouping conflicts with existing ",field," metadata; choose the existing mapping explicitly.")
        }
        if (!field %in% names(SummarizedExperiment::colData(object))) object[[field]] <- rep(NA, ncol(object))
        values <- object[[field]]
        values[idx] <- fields[[field]]
        object[[field]] <- values
      }
    }
    annotation <- file.path(app_root, if (parameters$organism == "Homo sapiens") "annotations_hsap.Rda" else "annotations_mmus.Rda")
    evidence <- if (is.null(parameters$plot_dir)) list(artifacts = character(), plots = data.frame()) else
      .legacy_qc_report(legacy, object, parameters$plot_dir)
    list(object = object, diagnostics = list(cells = ncol(object), features = nrow(object),
      excluded_cells = 0L, samples = length(indices), reused = "createMetadata",
      annotation_path = annotation, annotation_checksum = unname(tools::md5sum(annotation)),
      gene_audit = audit_legacy_qc_metrics(object, annotation), plots = evidence$plots),
      artifacts = evidence$artifacts)
  }
  filter_step <- function(object, parameters) {
    if (!inherits(object, "SingleCellExperiment")) stop("Existing QC filtering requires an SCE object.")
    required <- names(formals(legacy$filter_cells))[-c(1,2)]
    thresholds <- parameters$thresholds
    if (!is.list(thresholds) || !setequal(names(thresholds), required) || anyDuplicated(names(thresholds)) ||
        !all(vapply(thresholds, function(x) is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x), logical(1)))) {
      stop("Supply all nine explicit finite filtering thresholds.")
    }
    message("Calling existing filter_cells() with the recorded thresholds.")
    result <- do.call(legacy$filter_cells, c(list(raw_data = object,
      metrics = as.data.frame(SummarizedExperiment::colData(object))), thresholds))
    artifacts <- character()
    plots <- data.frame()
    if (!is.null(parameters$plot_dir)) {
      before <- .legacy_qc_report(legacy, object, file.path(parameters$plot_dir, "before"), thresholds)
      after <- .legacy_qc_report(legacy, result$se, file.path(parameters$plot_dir, "after"), thresholds)
      artifacts <- c(before$artifacts, after$artifacts)
      plots <- rbind(transform(before$plots, stage = "before"), transform(after$plots, stage = "after"))
    }
    list(object = result$se, diagnostics = list(cells_before = ncol(object), cells_after = ncol(result$se),
      excluded_cells = ncol(object) - ncol(result$se), reused = "filter_cells", message = result$mes,
      plots = plots), artifacts = artifacts)
  }
  list(metadata = metadata_step, filter = filter_step, source_files = source_files)
}

# Inspect the exact legacy definitions; do not alter gene sets or count values.
audit_legacy_qc_metrics <- function(object, annotation_path) {
  annotation <- new.env(parent = baseenv())
  load(annotation_path, envir = annotation)
  a <- annotation$annotations
  counts <- SummarizedExperiment::assay(object, "counts")
  sets <- list(mitochondrial = unique(a$gene_name[!is.na(a$seq_name) & a$seq_name == "MT"]),
    hemoglobin = unique(a$gene_name[a$description %in% grep("^hemoglobin", a$description,
                                                         value = TRUE, ignore.case = TRUE)]))
  genes <- do.call(rbind, lapply(names(sets), function(set) {
    names <- sets[[set]]
    rows <- if (set == "mitochondrial") a[!is.na(a$seq_name) & a$seq_name == "MT",] else
      a[a$description %in% grep("^hemoglobin", a$description, value = TRUE, ignore.case = TRUE),]
    details <- function(field) vapply(names, function(gene) paste(unique(rows[[field]][rows$gene_name == gene]), collapse = "; "), character(1))
    present <- names %in% rownames(counts)
    total <- rep(NA_real_, length(names))
    total[present] <- Matrix::rowSums(counts[names[present],,drop=FALSE])
    data.frame(set = set, gene_name = names, gene_biotype = details("gene_biotype"),
      gene_id = details("gene_id"), description = details("description"), present = present, total_umis = total)
  }))
  rownames(genes) <- NULL
  umi <- Matrix::colSums(counts)
  mt <- Matrix::colSums(counts[rownames(counts) %in% sets$mitochondrial,,drop=FALSE])
  hb <- Matrix::colSums(counts[rownames(counts) %in% sets$hemoglobin,,drop=FALSE])
  expected <- list(nUMI = umi, nGene = Matrix::colSums(counts > 0), mtUMI = mt, rbcUMI = hb,
    mitoRatio = mt/umi, rbcRatio = hb/umi, log10GenesPerUMI = log10(Matrix::colSums(counts > 0))/log10(umi))
  checks <- do.call(rbind, lapply(names(expected), function(field) {
    observed <- object[[field]]
    if (is.null(observed)) stop("QC field missing for audit: ", field)
    wanted <- expected[[field]]
    finite <- is.finite(observed) & is.finite(wanted)
    data.frame(metric = field, max_absolute_difference = if (any(finite)) max(abs(observed[finite]-wanted[finite])) else NA_real_,
      nonfinite_cells = sum(!is.finite(wanted)),
      nonfinite_agree = identical(unname(is.na(observed)),unname(is.na(wanted))) &&
        identical(unname(is.infinite(observed)),unname(is.infinite(wanted))))
  }))
  list(genes = genes, checks = checks, annotation_path = normalizePath(annotation_path),
    annotation_checksum = unname(tools::md5sum(annotation_path)),
    mitochondrial_definition = "Annotation seq_name == MT; includes coding genes, rRNAs and tRNAs",
    hemoglobin_definition = "Annotation descriptions starting with hemoglobin, case-insensitive; pseudogenes are not excluded",
    ratio_definition = "Per-cell selected-gene counts divided by all-gene counts; fraction 0..1, not percentage",
    annotation_build = "Release/genome build not recorded in bundled annotation table; hg19 correspondence unverified")
}

.legacy_qc_report <- function(legacy, object, directory, thresholds = NULL) {
  for (package in c("ggplot2", "scater")) {
    if (!requireNamespace(package, quietly = TRUE)) stop("Existing plotting package required: ", package)
  }
  if (dir.exists(directory)) stop("Use a new plot directory; existing evidence is preserved.")
  if (!dir.create(directory, recursive = TRUE)) stop("Cannot create plot directory.")
  directory <- normalizePath(directory)
  aliases <- c("ggplot","aes","geom_density","scale_x_log10","ylab","geom_vline","geom_boxplot",
    "ggtitle","geom_point","stat_smooth","scale_y_log10","geom_hline","facet_wrap")
  for (name in aliases) legacy[[name]] <- getExportedValue("ggplot2", name)
  legacy$lm <- stats::lm
  legacy$plotScater <- scater::plotScater
  legacy$plotHighestExprs <- scater::plotHighestExprs
  metrics <- as.data.frame(SummarizedExperiment::colData(object))
  bounds <- function(x) if (any(is.finite(x))) range(x[is.finite(x)]) else c(NA_real_, NA_real_)
  u <- bounds(metrics$nUMI); g <- bounds(metrics$nGene); m <- bounds(metrics$mitoRatio)
  r <- bounds(metrics$rbcRatio); n <- bounds(metrics$log10GenesPerUMI)
  if (!is.null(thresholds)) {
    u <- unlist(thresholds[c("minUMIsPerCell","maxUMIsPerCell")])
    g <- unlist(thresholds[c("minGenesPerCell","maxGenesPerCell")])
    m <- unlist(thresholds[c("minMitoPerCell","maxMitoPerCell")])
    r <- unlist(thresholds[c("minRbcPerCell","maxRbcPerCell")])
    n[1] <- thresholds$noveltyPerCell
  }
  calls <- list(createUMIsPlot = list(metrics,u[1],u[2]), createGenesHistPlot = list(metrics,g[1],g[2]),
    createGenesBoxPlot = list(metrics), createUMIsGenesPlot = list(metrics,g[1],g[2]),
    createUMIsGenesPlotRBC = list(metrics,g[1],g[2]), createmitoRatioPlot = list(metrics,m[1],m[2]),
    createrbcRatioPlot = list(metrics,r[1],r[2]), createnoveltyPlot = list(metrics,n[1],n[2]),
    createlibraryPlot = list(object), createhighExprPlot = list(object))
  metrics_path <- file.path(directory, "metrics.csv")
  utils::write.csv(metrics[, c("cells","sample","nUMI","nGene","log10GenesPerUMI","mtUMI","mitoRatio","rbcUMI","rbcRatio")],
                   metrics_path, row.names = FALSE)
  pdf_path <- file.path(directory, "existing-qc-plots.pdf")
  grDevices::pdf(pdf_path, width = 9, height = 6)
  on.exit(grDevices::dev.off())
  artifacts <- c(metrics_path, pdf_path)
  statuses <- list()
  for (name in names(calls)) {
    warnings <- character()
    path <- file.path(directory, paste0(name, ".png"))
    error <- tryCatch(withCallingHandlers({
      plot <- do.call(legacy[[name]], calls[[name]])
      ggplot2::ggsave(path, plot = plot, width = 9, height = 6, dpi = 110)
      print(plot)
      NULL
    }, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) conditionMessage(e))
    success <- is.null(error)
    message("  ", name, ": ", if (success) "rendered" else paste("FAILED", error))
    if (success) artifacts <- c(artifacts, path)
    statuses[[name]] <- data.frame(function_name = name, status = if (success) "rendered" else "failed",
      warnings = paste(unique(warnings), collapse = " | "), error = if (success) "" else error)
  }
  plots <- do.call(rbind, statuses)
  status_path <- file.path(directory, "plot-status.csv")
  utils::write.csv(plots, status_path, row.names = FALSE)
  list(artifacts = c(artifacts, status_path), plots = plots)
}

# A local view of the pending metadata proposal; decisions still go through R.
write_legacy_qc_review <- function(input, proposal, path) {
  if (file.exists(path)) stop("Use a new review report path.")
  escape <- function(x) {
    x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    gsub('"', "&quot;", x, fixed = TRUE)
  }
  object <- proposal$object
  before <- as.data.frame(SummarizedExperiment::colData(input))
  metrics <- as.data.frame(SummarizedExperiment::colData(object))
  geo_fields <- names(before)[startsWith(names(before), "geo_")]
  checks <- c("Counts unchanged" = identical(SummarizedExperiment::assay(input, "counts"), SummarizedExperiment::assay(object, "counts")),
    "Cell order unchanged" = identical(colnames(input), colnames(object)),
    "Gene identifiers unchanged" = identical(rownames(input), rownames(object)),
    "Full GEO provenance unchanged" = identical(S4Vectors::metadata(input)$geo, S4Vectors::metadata(object)$geo),
    "All GEO cell fields unchanged" = identical(before[,geo_fields,drop=FALSE], metrics[,geo_fields,drop=FALSE]),
    "No cells removed" = ncol(input) == ncol(object))
  plots <- proposal$diagnostics$plots
  plot_rows <- paste0("<tr><td>", escape(plots$function_name), "</td><td>", escape(plots$status),
    "</td><td>", escape(paste(plots$warnings, plots$error)), "</td></tr>", collapse = "")
  check_rows <- paste0("<li>", ifelse(checks, "PASS", "FAIL"), " · ", escape(names(checks)), "</li>", collapse = "")
  images <- proposal$artifacts$path[grepl("\\.png$", proposal$artifacts$path)]
  titles <- c(createUMIsPlot = "UMI distribution", createGenesHistPlot = "Detected-gene distribution",
    createGenesBoxPlot = "Detected genes per sample", createUMIsGenesPlot = "UMIs vs genes · mitochondrial fraction",
    createUMIsGenesPlotRBC = "UMIs vs genes · hemoglobin fraction", createmitoRatioPlot = "Mitochondrial fraction",
    createrbcRatioPlot = "Hemoglobin fraction", createnoveltyPlot = "Genes per UMI complexity",
    createlibraryPlot = "Library overview", createhighExprPlot = "Highest expressed genes",
    "cell-cycle-scores" = "Cell-cycle scores and predicted phases", "doublet-scores" = "Doublet predictions",
    "ambient-decontx-fraction" = "Ambient RNA · DecontX model estimates")
  panels <- vapply(images, function(image) {
    name <- sub("\\.png$", "", basename(image))
    title <- if (name %in% names(titles)) titles[[name]] else name
    expanded <- name %in% if (any(grepl("ambient-decontx-fraction\\.png$",images))) {
      c("ambient-decontx-fraction","cell-cycle-scores")
    } else if (any(grepl("cell-cycle-scores\\.png$",images))) {
      c("createUMIsGenesPlot","cell-cycle-scores")
    } else c("createUMIsGenesPlot","createmitoRatioPlot")
    paste0("<details", if (expanded) " open" else "", "><summary>", escape(title), "</summary><img src=\"",
      escape(image), "\" alt=\"", escape(title), "\"></details>")
  }, character(1))
  optional <- proposal$diagnostics$optional_checks
  optional_html <- if (is.null(optional)) "" else paste0("<table><tr><th>Optional QC check</th><th>Status</th><th>Notes</th></tr>",
    paste0("<tr><td>",escape(optional$check),"</td><td>",escape(optional$status),"</td><td>",escape(optional$notes),"</td></tr>",collapse=""),"</table>")
  html <- paste0('<!doctype html><html><head><meta charset="utf-8"><style>
body{font-family:var(--font-sans);color:var(--foreground);margin:0}p{line-height:1.5}
.cards{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:10px}.card,details{background:var(--card);border:1px solid var(--border);border-radius:var(--radius);padding:16px}
.value{font-size:24px;font-weight:650}.muted{color:var(--muted-foreground)}summary{cursor:pointer;font-weight:600}details{margin:10px 0}img{display:block;width:100%;height:auto;margin-top:10px}table{border-collapse:collapse;width:100%;font-size:12px}td,th{padding:8px;border-bottom:1px solid var(--border);text-align:left;overflow-wrap:anywhere}li{margin:7px 0}
</style></head><body><p><strong>', escape(paste(unique(metrics$sample), collapse = ", ")),
    '</strong> · Existing app QC · <strong>awaiting your decision</strong></p><div class="cards">',
    '<div class="card"><div class="muted">Cells retained</div><div class="value">', ncol(object),
    '</div></div><div class="card"><div class="muted">Median UMIs</div><div class="value">', stats::median(metrics$nUMI),
    '</div></div><div class="card"><div class="muted">Median mitochondrial fraction</div><div class="value">',
    sprintf('%.2f%%',100*stats::median(metrics$mitoRatio)), '</div></div></div>',
    '<p>No filtering has run. The initial <code>nUMI &gt;= 100</code> rule would remove <strong>',
    sum(metrics$nUMI < 100), ' cells</strong>. Plot reference lines mark observed ranges, not recommended cutoffs.</p>',
    optional_html,'<details><summary>Input integrity checks</summary><ul>', check_rows, '</ul></details>',
    '<details><summary>Exact functions reused and plotting warnings</summary><table><thead><tr><th>Function</th><th>Result</th><th>Warnings / errors</th></tr></thead><tbody>',
    plot_rows, '</tbody></table></details>', paste(panels, collapse = "\n"),
    '<p class="muted">Review the metric calculations and plots first. Acceptance records a reason in the R history; it does not approve any filtering thresholds. Run: ',
    escape(proposal$run_id), '</p></body></html>')
  writeLines(html, path)
  invisible(normalizePath(path))
}

write_legacy_gene_review <- function(audit, directory) {
  if (dir.exists(directory)) stop("Use a new gene-review directory.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Existing ggplot2 installation required.")
  dir.create(directory, recursive = TRUE)
  directory <- normalizePath(directory)
  genes <- audit$genes
  utils::write.csv(genes, file.path(directory, "gene-sets.csv"), row.names = FALSE)
  utils::write.csv(audit$checks, file.path(directory, "metric-checks.csv"), row.names = FALSE)
  pooled <- stats::aggregate(total_umis ~ gene_biotype, genes[genes$set == "mitochondrial" & genes$present,], sum)
  mt_plot <- ggplot2::ggplot(pooled, ggplot2::aes(x = gene_biotype, y = total_umis)) +
    ggplot2::geom_col(fill = "#397a97") + ggplot2::coord_flip() + ggplot2::theme_bw() +
    ggplot2::labs(x = NULL, y = "Pooled mitochondrial UMIs", title = "Contributions to the existing MT gene set")
  hb_plot <- ggplot2::ggplot(genes[genes$set == "hemoglobin" & genes$present,],
    ggplot2::aes(x = gene_name, y = total_umis)) + ggplot2::geom_col(fill = "#9c5b62") +
    ggplot2::theme_bw() + ggplot2::labs(x = NULL, y = "Pooled hemoglobin UMIs", title = "Present hemoglobin genes")
  mt_path <- file.path(directory, "mitochondrial-contributions.png")
  hb_path <- file.path(directory, "hemoglobin-contributions.png")
  ggplot2::ggsave(mt_path, plot = mt_plot, width = 8, height = 3.3, dpi = 110)
  ggplot2::ggsave(hb_path, plot = hb_plot, width = 8, height = 3.3, dpi = 110)
  escape <- function(x) {
    x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    gsub('"', "&quot;", x, fixed = TRUE)
  }
  tables <- vapply(c("mitochondrial","hemoglobin"), function(set) {
    g <- genes[genes$set == set,]
    rows <- paste0("<tr><td>",escape(g$gene_name),"</td><td>",escape(g$gene_biotype),"</td><td>",
      ifelse(g$present,"present","absent"),"</td><td>",ifelse(g$present,g$total_umis,"—"),"</td></tr>",collapse="")
    paste0("<details><summary>",escape(set)," · ",sum(g$present)," present / ",nrow(g)," defined</summary>",
      "<table><tr><th>Gene</th><th>Biotype</th><th>In matrix?</th><th>Pooled UMIs</th></tr>",rows,"</table></details>")
  },character(1))
  checks <- audit$checks
  check_rows <- paste0("<tr><td>",escape(checks$metric),"</td><td>",checks$max_absolute_difference,
    "</td><td>",checks$nonfinite_cells,"</td></tr>",collapse="")
  html <- paste0('<!doctype html><html><head><meta charset="utf-8"><style>
body{font-family:var(--font-sans);color:var(--foreground);margin:0}p{line-height:1.5}img{width:100%;height:auto;display:block}details{background:var(--card);border:1px solid var(--border);border-radius:var(--radius);padding:16px;margin:10px 0}summary{cursor:pointer;font-weight:600}table{width:100%;border-collapse:collapse;font-size:13px}td,th{text-align:left;padding:7px;border-bottom:1px solid var(--border)}code{overflow-wrap:anywhere}.note{color:var(--muted-foreground)}</style></head><body>',
    '<p><strong>Exact gene-set audit · definitions unchanged · proposal still pending</strong></p>',
    '<p>Mitochondrial genes come from <code>seq_name == "MT"</code>. Hemoglobin genes come from descriptions matching <code>^hemoglobin</code>, case-insensitively. The latter definition does not exclude pseudogenes.</p>',
    '<p><code>mtUMI = sum(MT counts per cell)</code><br><code>mitoRatio = mtUMI / nUMI</code><br>',
    '<code>rbcUMI = sum(hemoglobin counts per cell)</code><br><code>rbcRatio = rbcUMI / nUMI</code></p>',
    '<p>Ratios are fractions (0–1), not percentages. <code>rbcRatio</code> measures hemoglobin transcripts, not the fraction of red blood cells.</p>',
    '<details><summary>Independent per-cell metric checks</summary><table><tr><th>Metric</th><th>Largest difference</th><th>Nonfinite cells</th></tr>',
    check_rows,'</table></details>',paste(tables,collapse="\n"),'<img src="',escape(mt_path),'" alt="Mitochondrial contributions">',
    '<img src="',escape(hb_path),'" alt="Hemoglobin contributions"><p class="note">Bars show pooled counts across this pilot, not average per-cell fractions. Absent genes have no row in the supplied matrix; absence is not a measured zero.</p>',
    '<p><strong>Annotation provenance needs review:</strong> ',escape(audit$annotation_build),'</p>',
    '<details><summary>Annotation file identity</summary><p><code>',escape(audit$annotation_path),'</code></p><p>MD5: <code>',
    escape(audit$annotation_checksum),'</code></p></details></body></html>')
  path <- file.path(directory,"gene-review.html")
  writeLines(html,path)
  invisible(path)
}
