# Maintenance has its own state: checking/downloading never invalidates QC.
# Only choosing a different reference/build changes the annotation configuration.
annotation_reference_server <- function(input,output,session,app_root=".",cache_dir=NULL,
    catalog_fetch=NULL,resource_fetch=NULL) {
  api <- .annotation_cache_api(app_root)
  if (is.null(cache_dir)) cache_dir <- api$annotation_cache_dir(app_root)
  if (is.null(catalog_fetch)) catalog_fetch <- api$annotationhub_catalog
  if (is.null(resource_fetch)) resource_fetch <- api$annotationhub_fetch
  checked <- shiny::reactiveVal(NULL)
  message <- shiny::reactiveVal("")
  busy <- shiny::reactiveVal(FALSE)
  library_revision <- shiny::reactiveVal(0L)
  revision <- 0L
  target <- shiny::reactive(list(organism=input$typeOrganism,build=input$annotationDataBuild))
  configuration <- shiny::reactive(list(reference_id=if (is.null(input$annotationReference) ||
    !nzchar(input$annotationReference)) "bundled" else input$annotationReference,
    data_genome_build=if (is.null(input$annotationDataBuild)) "" else input$annotationDataBuild))

  output$annotationReferenceControls <- shiny::renderUI({
    organism <- input$typeOrganism
    if (length(organism)!=1L || !organism %in% c("Homo sapiens","Mus musculus")) return(NULL)
    library_revision()
    entries <- api$list_annotation_references(organism,app_root,cache_dir)
    selected <- shiny::isolate(input$annotationReference)
    if (is.null(selected) || !selected %in% names(entries)) selected <- "bundled"
    shiny::selectInput("annotationReference","Annotation reference",
      choices=stats::setNames(names(entries),vapply(entries,`[[`,character(1),"label")),selected=selected)
  })
  output$annotationGenomeControls <- shiny::renderUI({
    organism <- input$typeOrganism
    if (length(organism)!=1L || !organism %in% c("Homo sapiens","Mus musculus")) return(NULL)
    builds <- if (organism=="Homo sapiens") c("GRCh37 / hg19"="GRCh37","GRCh38 / hg38"="GRCh38") else
      c("GRCm38 / mm10"="GRCm38","GRCm39 / mm39"="GRCm39")
    shiny::selectInput("annotationDataBuild","Counts reference genome build",
      choices=c("Unknown"="",builds),selected="")
  })
  output$annotationUpdateStatus <- shiny::renderText(message())
  output$annotationUpdateBusy <- shiny::renderUI({
    if (busy()) shiny::tags$div(role="status",shiny::icon("spinner",class="fa-spin")," Checking or downloading annotation reference…")
  })
  output$annotationCandidates <- shiny::renderTable({
    state <- checked(); if (is.null(state) || is.null(state$candidates)) return(NULL)
    state$candidates[,intersect(c("resource_id","title","species","genome","rdatadateadded"),names(state$candidates)),drop=FALSE]
  },rownames=FALSE)
  output$annotationDownloadControls <- shiny::renderUI({
    state <- checked()
    if (is.null(state) || !is.null(state$error) || is.null(state$candidates) || nrow(state$candidates)==0L) return(NULL)
    d <- state$candidates
    shiny::tagList(
      shiny::selectInput("annotationCandidate","Reference to download",
        choices=c("Select a resource"="",stats::setNames(d$resource_id,paste(d$resource_id,d$title))),selected=""),
      shiny::textInput("annotationUpdateReason","Reason for adding this reference"),
      shiny::checkboxInput("annotationDownloadConsent","Download this resource into the project cache (may require hundreds of MB)",FALSE),
      shiny::actionButton("annotationDownload","Download as a separate reference"))
  })
  shiny::observeEvent(target(),{
    revision <<- revision+1L
    checked(NULL); busy(FALSE)
    t <- target()
    state <- tryCatch(api$annotation_update_state(t$organism,t$build,cache_dir),error=function(e) NULL)
    if (is.null(state)) {
      message("Choose the counts' species and genome build to check for compatible references. The bundled export remains usable with its build marked unknown.")
    } else {
      checked(state)
      message(paste(if (state$due) "Monthly annotation check is due." else paste("Next monthly check:",state$next_check),
        if (!is.na(state$last_success)) paste("Last successful catalogue check:",state$last_success) else "No successful check recorded.",
        if (!is.null(state$error)) paste("Last attempt failed:",state$error) else "",
        "Use Check for reference updates; normal annotation loading stays offline."))
    }
  },ignoreNULL=FALSE,priority=60)
  shiny::observeEvent(input$annotationCheckUpdates,{
    t <- target()
    valid <- tryCatch({api$.annotation_target(t$organism,t$build);TRUE},error=function(e) {message(conditionMessage(e));FALSE})
    if (!valid) return()
    revision <<- revision+1L; token <- revision
    busy(TRUE);message("Checking AnnotationHub catalogue metadata; existing QC results and references remain active…")
    session$onFlushed(function() shiny::isolate({
      if (session$isClosed() || token!=revision) return(invisible(NULL))
      tryCatch({
        state <- api$check_annotation_updates(t$organism,t$build,app_root,cache_dir,force=TRUE,catalog_fetch=catalog_fetch)
        if (token!=revision) return(invisible(NULL))
        checked(state)
        message(if (!is.null(state$error)) paste("Catalogue check failed:",state$error,
          "Previous successful results, if shown, are stale. No reference or QC result changed.") else
          paste("Checked snapshot",state$snapshot_date,"on",state$last_success,"—",nrow(state$candidates),
            "compatible Ensembl EnsDb candidates. The bundled export's release is unknown, so these are candidates for review, not verified newer versions. No reference selected or QC result changed."))
      },error=function(e) message(paste("Catalogue check stopped:",conditionMessage(e))))
      if (token==revision) busy(FALSE)
    }),once=TRUE)
  },ignoreInit=TRUE)
  shiny::observeEvent(input$annotationDownload,{
    state <- checked(); id <- input$annotationCandidate
    if (is.null(state) || !is.null(state$error) || length(id)!=1L || !id %in% state$candidates$resource_id) {
      message("Check the catalogue and select exactly one compatible reference first.");return()
    }
    reason <- input$annotationUpdateReason
    if (!isTRUE(input$annotationDownloadConsent) || length(reason)!=1L || is.na(reason) || !nzchar(trimws(reason))) {
      message("Confirm the reference download and provide a reason. No download started.");return()
    }
    candidate <- state$candidates[state$candidates$resource_id==id,,drop=FALSE]
    revision <<- revision+1L; token <- revision
    busy(TRUE);message(paste("Downloading",id,"as a separate reference. Prepared QC results keep their existing reference…"))
    session$onFlushed(function() shiny::isolate({
      if (session$isClosed() || token!=revision) return(invisible(NULL))
      tryCatch({
        entry <- api$download_annotation_reference(candidate,state$snapshot_date,reason,app_root,cache_dir,resource_fetch)
        library_revision(library_revision()+1L)
        message(paste("Saved",entry$label,"in the project cache. To use it, choose it under Annotation reference, then Add gene annotations. Existing QC results still use their recorded reference."))
      },error=function(e) message(paste("Reference download stopped:",conditionMessage(e),"Existing references and QC results are unchanged.")))
      if (token==revision) busy(FALSE)
    }),once=TRUE)
  },ignoreInit=TRUE)
  list(configuration=configuration,checked=shiny::reactive(checked()),status=shiny::reactive(message()),cache_dir=cache_dir)
}
