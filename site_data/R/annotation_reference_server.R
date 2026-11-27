# Maintenance has its own state: checking/downloading never invalidates QC.
# Only choosing a different reference/build changes the annotation configuration.
annotation_reference_server <- function(input,output,session,app_root=".",cache_dir=NULL,
    catalog_fetch=NULL,resource_fetch=NULL,supplement_fetch=NULL) {
  api <- .annotation_cache_api(app_root)
  if (is.null(cache_dir)) cache_dir <- api$annotation_cache_dir(app_root)
  if (is.null(catalog_fetch)) catalog_fetch <- api$annotationhub_catalog
  if (is.null(resource_fetch)) resource_fetch <- api$annotationhub_fetch
  if (is.null(supplement_fetch)) supplement_fetch <- api$grch37_biomart_fetch
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
    shiny::tagList(shiny::selectInput("annotationReference","Annotation reference",
      choices=stats::setNames(names(entries),vapply(entries,`[[`,character(1),"label")),selected=selected),
      shiny::helpText("Choose a cached reference matching the counts' genome build. To add another build, choose it above and use Annotation reference maintenance. The bundled export has an unverified build."))
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
    state$candidates[,intersect(c("resource_id","title","species","genome","rdataclass","rdatadateadded"),names(state$candidates)),drop=FALSE]
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
      shiny::helpText("GTF references provide gene IDs, names and biotypes. Descriptions and Entrez links are unavailable in the GTF itself. For human GRCh37 Ensembl 75, use the separate BioMart supplement controls to restore available missing fields. Without these fields, hemoglobin QC and some ribosomal protein flags remain unavailable."),
      shiny::actionButton("annotationDownload","Download as a separate reference"))
  })
  output$annotationSupplementControls <- shiny::renderUI({
    library_revision()
    if (!identical(input$typeOrganism,"Homo sapiens") || !identical(input$annotationDataBuild,"GRCh37")) return(NULL)
    entries <- api$list_annotation_references("Homo sapiens",app_root,cache_dir)
    eligible <- entries[vapply(entries,function(e) {
      p <- e$provenance
      identical(p$source,"AnnotationHub Ensembl GTF") && identical(p$genome_build,"GRCh37") &&
        identical(p$annotation_release,"75") && is.null(p$supplement)
    },logical(1))]
    if (!length(eligible)) return(NULL)
    shiny::tagList(shiny::h5("Complete missing GRCh37 annotation fields"),
      shiny::helpText("Retrieve descriptions and Entrez links from Ensembl GRCh37 BioMart, verified as Ensembl 75. Save a separate reference; counts and prepared QC are unchanged. Ordinary annotation loading stays offline."),
      shiny::selectInput("annotationSupplementParent","GRCh37 reference to supplement",
        choices=stats::setNames(names(eligible),vapply(eligible,`[[`,character(1),"label"))),
      shiny::textInput("annotationSupplementReason","Reason for supplementing this reference"),
      shiny::checkboxInput("annotationSupplementConsent","Retrieve BioMart annotations and save a separate local reference",FALSE),
      shiny::actionButton("annotationSupplement","Add missing fields from GRCh37 BioMart"))
  })
  selected_supplement <- shiny::reactive({
    library_revision()
    if (!identical(input$typeOrganism,"Homo sapiens")) return(NULL)
    e <- api$list_annotation_references("Homo sapiens",app_root,cache_dir)[[configuration()$reference_id]]
    if (is.null(e$provenance$supplement)) NULL else e
  })
  output$annotationSupplementCoverage <- shiny::renderTable({
    e <- selected_supplement(); if (is.null(e)) return(NULL)
    s <- e$provenance$supplement$summary
    data.frame(measure=names(s),reference_rows=as.integer(unlist(s)),row.names=NULL)
  },rownames=FALSE)
  output$annotationSupplementReviewControls <- shiny::renderUI({
    if (is.null(selected_supplement())) return(NULL)
    shiny::tagList(shiny::helpText("BioMart filled missing fields only. Multiple Entrez links are retained; conflicting descriptions remain unresolved. Download the complete per-row supplement review."),
      shiny::downloadButton("annotationSupplementMappingDownload","Download BioMart supplement review"))
  })
  output$annotationSupplementMappingDownload <- shiny::downloadHandler(
    filename=function() "grch37-biomart-supplement-review.csv",content=function(path) {
      e <- shiny::req(selected_supplement())
      d <- api$load_annotation_supplement_review(e$reference_id,app_root,cache_dir)
      utils::write.csv(d,path,row.names=FALSE,na="")
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
            "compatible Ensembl gene-annotation candidates (EnsDb or GTF). The bundled export's release is unknown, so these are candidates for review, not verified newer versions. No reference selected or QC result changed."))
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
  shiny::observeEvent(input$annotationSupplement,{
    reason <- input$annotationSupplementReason; id <- input$annotationSupplementParent
    if (!identical(input$typeOrganism,"Homo sapiens") || !identical(input$annotationDataBuild,"GRCh37") ||
        length(id)!=1L || is.na(id) || !isTRUE(input$annotationSupplementConsent) ||
        length(reason)!=1L || is.na(reason) || !nzchar(trimws(reason))) {
      message("Choose a human GRCh37 reference, confirm BioMart retrieval and provide a reason. No supplement requested.");return()
    }
    revision <<- revision+1L; token <- revision
    busy(TRUE);message("Retrieving verified GRCh37 BioMart fields; existing references and prepared QC remain active…")
    session$onFlushed(function() shiny::isolate({
      if (session$isClosed() || token!=revision) return(invisible(NULL))
      tryCatch({
        e <- api$supplement_grch37_reference(id,reason,app_root,cache_dir,supplement_fetch)
        library_revision(library_revision()+1L)
        message(paste("Saved",e$label,"as a separate reference. Select it under Annotation reference, then Add gene annotations and review QC again. Existing results remain pinned to their original reference."))
      },error=function(e) message(paste("BioMart supplement stopped:",conditionMessage(e),"Existing references and QC results are unchanged.")))
      if (token==revision) busy(FALSE)
    }),once=TRUE)
  },ignoreInit=TRUE)
  list(configuration=configuration,checked=shiny::reactive(checked()),status=shiny::reactive(message()),cache_dir=cache_dir)
}
