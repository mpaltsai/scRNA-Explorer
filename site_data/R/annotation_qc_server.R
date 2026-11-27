# Annotation/metric preparation is separate from baseline construction. Prepared
# tables contain no expression matrix; the baseline SCE is made only on request.
annotation_qc_server <- function(input,output,session,records,app_root=".",reference_services=list()) {
  # Test seams replace catalogue/resource access, never the QC implementation.
  references <- do.call(annotation_reference_server,c(list(input=input,output=output,session=session,app_root=app_root),reference_services))
  prepared <- shiny::reactiveVal(NULL)
  baseline <- shiny::reactiveVal(NULL)
  status <- shiny::reactiveVal("Read data, choose the species and identifier types, then add annotations.")
  baseline_status <- shiny::reactiveVal("")
  busy <- shiny::reactiveVal(FALSE)
  revision <- 0L
  options <- shiny::reactive({
    x <- records(); if (is.null(x)) return(NULL)
    keys <- vapply(seq_along(x),function(i) {
      value <- input[[paste0("annotationKey_",i)]]
      # 10X column 1 declares primary IDs. For tables/RDS the user must choose;
      # never guess species or identifier type from an ENSG prefix.
      if (is.null(value) && identical(x[[i]]$source$format,"10X")) value <- "gene_id"
      if (is.null(value) || length(value)!=1L) "" else value
    },character(1))
    names(keys) <- names(x)
    c(list(organism=input$typeOrganism,keys=keys,strip_ensembl_version=isTRUE(input$annotationStripVersion)),references$configuration())
  })
  output$annotationIdControls <- shiny::renderUI({
    x <- records(); if (is.null(x)) return(NULL)
    shiny::tagList(lapply(seq_along(x),function(i) shiny::selectInput(paste0("annotationKey_",i),
      paste("Feature identifiers:",names(x)[i]),
      choices=c("Select identifier type"="","Gene symbols"="gene_name","Ensembl gene IDs"="gene_id"),
      selected=if (identical(x[[i]]$source$format,"10X")) "gene_id" else "")))
  })
  output$paramsAddMetadata <- shiny::renderText(status())
  output$initDataFilt <- shiny::renderText(baseline_status())
  output$annotationBusy <- shiny::renderUI({
    if (busy()) shiny::tags$div(role="status",shiny::icon("spinner",class="fa-spin")," Preparing annotations and QC metrics…")
  })
  output$annotationSummary <- shiny::renderTable({
    p <- prepared(); if (is.null(p)) return(NULL)
    do.call(rbind,lapply(p,`[[`,"summary"))
  },rownames=FALSE)
  output$ribosomalSummary <- shiny::renderTable({
    p <- prepared(); if (is.null(p)) return(NULL)
    do.call(rbind,lapply(p,`[[`,"ribosomal_summary"))
  },rownames=FALSE)
  output$annotationProvenance <- shiny::renderTable({
    p <- prepared(); if (is.null(p)) return(NULL)
    r <- p[[1]]$config$reference
    value <- function(field) {x<-r[[field]];if(is.null(x) || !length(x) || is.na(x)) "Unknown" else x}
    data.frame(species=p[[1]]$config$organism,cache=basename(r$path),md5=r$md5,
      resource_id=value("resource_id"),snapshot=value("snapshot_date"),
      annotation_release=value("annotation_release"),genome_build=value("genome_build"),
      source=r$source,reference_id=value("reference_id"),data_genome_build=p[[1]]$config$data_genome_build,
      ribosomal_source=if (identical(p[[1]]$config$ribosomal_reference$status,"available")) {
        paste(p[[1]]$config$ribosomal_reference$package,p[[1]]$config$ribosomal_reference$package_version)
      } else "Protein reference unavailable; RNA biotypes only")
  },rownames=FALSE)
  output$annotationIssues <- shiny::renderTable({
    p <- prepared(); if (is.null(p)) return(NULL)
    rows <- lapply(names(p),function(n) {
      d <- p[[n]]$feature_annotation
      d <- d[d$annotation_status!="matched" | d$ribosomal_annotation_status=="unresolved",,drop=FALSE]
      cbind(sample=rep(n,nrow(d)),d[,c("source_feature_id","annotation_status","annotation_candidates","ribosomal_annotation_status"),drop=FALSE])
    })
    utils::head(do.call(rbind,rows),20L)
  },rownames=FALSE)
  output$annotationMappingDownload <- shiny::downloadHandler(filename=function() "gene-annotation-mapping.csv",content=function(path) {
    p <- shiny::req(prepared()); x <- shiny::req(records())
    d <- do.call(rbind,lapply(names(p),function(n) data.frame(sample=n,
      internal_feature_id=rownames(x[[n]]$counts),p[[n]]$feature_annotation,check.names=FALSE)))
    utils::write.csv(d,path,row.names=FALSE)
  })
  output$qcMetricsDownload <- shiny::downloadHandler(filename=function() "qc-metrics-before-initial-filter.csv",content=function(path) {
    p <- shiny::req(prepared())
    utils::write.csv(do.call(rbind,lapply(p,`[[`,"metrics")),path,row.names=FALSE)
  })
  output$baselineSummary <- shiny::renderTable({
    b <- baseline(); if (is.null(b)) return(NULL)
    do.call(rbind,lapply(names(b),function(n) {
      info <- S4Vectors::metadata(b[[n]])$scrna_explorer_qc$initial_filter
      data.frame(sample=n,input_cells=info$input_cells,baseline_cells=info$retained_cells,
        excluded_below_100=length(info$excluded_cells))
    }))
  },rownames=FALSE)
  shiny::observeEvent(list(records(),options()),{
    revision <<- revision+1L
    prepared(NULL); baseline(NULL); busy(FALSE); baseline_status("")
    status(if (is.null(records())) "Read data before adding annotations." else "Annotation settings not yet applied. Choose species and identifier types, then add annotations.")
  },ignoreNULL=FALSE,priority=50)
  shiny::observeEvent(input$addMetadata,{
    revision <<- revision+1L
    token <- revision
    x <- records(); cfg <- options()
    prepared(NULL); baseline(NULL); baseline_status("")
    if (is.null(x)) {status("Read and review data first.");return()}
    if (length(cfg$organism)!=1L || !cfg$organism %in% c("Homo sapiens","Mus musculus") || any(!cfg$keys %in% c("gene_id","gene_name"))) {
      status("Select a species and an identifier type for every sample.");return()
    }
    busy(TRUE); status("Preparing cached gene annotations and QC metrics…")
    session$onFlushed(function() shiny::isolate({
      if (session$isClosed() || token!=revision) return(invisible(NULL))
      tryCatch({
        reference <- load_gene_annotation_cache(cfg$organism,app_root,reference_id=cfg$reference_id,
          cache_dir=references$cache_dir)
        known_build <- reference$provenance$genome_build
        if (nzchar(cfg$data_genome_build) && !is.na(known_build) &&
            !identical(known_build,cfg$data_genome_build)) stop("Selected reference build differs from the declared counts build. Choose a compatible reference.")
        if (cfg$reference_id!="bundled" && !nzchar(cfg$data_genome_build))
          stop("Declare the counts genome build before using a downloaded reference.")
        result <- stats::setNames(lapply(names(x),function(n)
          prepare_sample_qc(x[[n]],reference,cfg$keys[[n]],cfg$strip_ensembl_version)),names(x))
        for (n in names(result)) result[[n]]$config$data_genome_build <- if (nzchar(cfg$data_genome_build)) cfg$data_genome_build else "Unknown"
        if (token!=revision) return(invisible(NULL))
        prepared(result)
        missing_hemoglobin <- any(vapply(result,function(p) identical(p$summary$hemoglobin_status,"missing_reference_definition"),logical(1)))
        status(paste("Annotations and QC metrics prepared. Review mapping coverage and reference provenance, then run Minimal filtering (UMIs >= 100). Original feature IDs and counts are preserved; unmapped/ambiguous features were not removed.",
          if (missing_hemoglobin) "Hemoglobin QC is unavailable: this reference lacks the description-based gene definition. QC filtering remains blocked until that evidence is supplied; unavailable values are not zero contamination." else ""))
      },error=function(e) {if(token==revision) status(paste("Annotation stopped:",conditionMessage(e)))})
      if (token==revision) busy(FALSE)
    }),once=TRUE)
  },ignoreInit=TRUE)
  shiny::observeEvent(input$initFiltering,{
    x <- records(); p <- prepared()
    if (is.null(x) || is.null(p)) {baseline_status("Prepare and review gene annotations first.");return()}
    tryCatch({
      b <- stats::setNames(lapply(names(x),function(n) make_baseline_sce(x[[n]],p[[n]])),names(x))
      baseline(b)
      empty <- names(b)[vapply(b,ncol,integer(1))==0L]
      baseline_status(paste("Minimal filtering completed: retained UMIs >= 100.",
        if(length(empty)) paste("No baseline cells in:",paste(empty,collapse=", ")) else ""))
    },error=function(e) {baseline(NULL);baseline_status(paste("Baseline construction stopped:",conditionMessage(e)))})
  },ignoreInit=TRUE,priority=10)
  list(prepared=shiny::reactive(prepared()),baseline=shiny::reactive(baseline()),
    status=shiny::reactive(status()),baseline_status=shiny::reactive(baseline_status()))
}
