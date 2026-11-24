# Session-local orchestration for the input stage. request(directory, selection)
# returns import records, not filtered data. selection is a snapshot taken on
# Read data; later changes cannot silently change the meaning of a queued read.
# This small controller can be tested without
# starting the clustering/correlation sections of the application.
input_import_server <- function(input,output,session,request) {
  directory <- tempfile("scrna-import-session-")
  dir.create(directory)
  # Only this session's generated working directory is removed on session end.
  session$onSessionEnded(function() unlink(directory,recursive=TRUE))
  records <- shiny::reactiveVal(NULL)
  status <- shiny::reactiveVal("Read data to begin.")
  busy <- shiny::reactiveVal(FALSE)
  revision <- 0L
  selection <- shiny::reactive({
    # Subscribe only to controls relevant to the selected input. Hidden upload
    # controls, annotation settings and decision reasons must not reset import.
    mode <- input$testData
    if (identical(as.character(mode),"1")) return(list(testData=mode))
    type <- input$typeData
    key <- if (length(type)==1L) as.character(type) else ""
    fields <- switch(key,"1"="csvFile","2"=c("mtxInput","genes10XInput","barcodes10XInput"),
      "3"="rdsFile","4"="zipFile",character())
    result <- list(testData=mode,typeData=type)
    for (field in fields) result[field] <- list(input[[field]])
    result
  })
  ready <- shiny::reactive({
    x <- records()
    if (is.null(x) || !length(x) || !all(vapply(x,function(r) identical(r$status,"ready"),logical(1)))) return(NULL)
    x
  })
  decisions <- shiny::reactive({
    rows <- list()
    for (r in records()) for (d in r$decisions) rows[[length(rows)+1L]] <- data.frame(
      sample=r$sample_id,action=d$action,reason=d$reason,
      affected_rows=paste(d$affected_rows,collapse=","),time=d$time)
    if (length(rows)) do.call(rbind,rows) else data.frame(sample=character(),action=character(),
      reason=character(),affected_rows=character(),time=character())
  })
  output$readData <- shiny::renderText(status())
  output$importBusy <- shiny::renderUI({
    if (!busy()) return(NULL)
    shiny::tags$div(role="status",`aria-live`="polite",
      shiny::icon("spinner",class="fa-spin")," Reading data…")
  })
  output$importSummary <- shiny::renderTable({
    x <- records(); if (is.null(x)) return(NULL)
    do.call(rbind,lapply(x,function(r) data.frame(sample=r$sample_id,features=nrow(r$counts),
      cells=ncol(r$counts),generated_cell_ids=sum(r$cell_identity$generated),status=r$status)))
  },rownames=FALSE)
  output$importDuplicateTable <- shiny::renderTable({
    x <- records(); if (is.null(x)) return(NULL)
    pending <- x[vapply(x,function(r) identical(r$status,"needs_feature_decision"),logical(1))]
    if (!length(pending)) return(NULL)
    do.call(rbind,lapply(pending,import_duplicate_rows))
  },rownames=FALSE)
  output$importDecisionControls <- shiny::renderUI({
    x <- records(); if (is.null(x)) return(NULL)
    pending <- names(x)[vapply(x,function(r) identical(r$status,"needs_feature_decision"),logical(1))]
    if (!length(pending)) return(NULL)
    shiny::tagList(shiny::helpText("Duplicate primary feature IDs need your decision. Keep rows separately to preserve every count, or cancel and upload corrected data. This version does not merge or discard rows."),
      shiny::selectInput("importDecisionSample","Sample",choices=pending),
      shiny::selectInput("importDecisionAction","Decision",choices=c("Keep all rows separately"="keep_separate","Cancel import; correct the input"="cancel")),
      shiny::textInput("importDecisionReason","Reason"),
      shiny::actionButton("resolveImport","Record import decision"))
  })
  output$importDecisionHistory <- shiny::renderTable(decisions(),rownames=FALSE)
  output$importDecisionDownload <- shiny::downloadHandler(
    filename=function() "import-decisions.csv",
    content=function(file) utils::write.csv(decisions(),file,row.names=FALSE))
  shiny::observeEvent(selection(),{
    revision <<- revision+1L
    records(NULL)
    busy(FALSE)
    status("Not yet read. Input selection changed; click Read data.")
  },ignoreInit=TRUE,priority=100)
  shiny::observeEvent(input$startQC,{
    revision <<- revision+1L
    token <- revision
    selected_input <- selection()
    records(NULL)
    busy(TRUE)
    status("Reading data… Please wait.")
    # Send the cleared results and busy state first. Merely setting a label and
    # then doing synchronous I/O in this observer would delay the visible update.
    # Reading itself remains synchronous; this does not claim faster file I/O.
    session$onFlushed(function() shiny::isolate({
      if (session$isClosed() || token!=revision) return(invisible(NULL))
      tryCatch({
        x <- request(directory,selected_input)
        if (token!=revision || session$isClosed()) return(invisible(NULL))
        if (!length(x) || is.null(names(x)) || any(!nzchar(names(x))) || anyDuplicated(names(x))) stop("Each imported sample needs a unique name.")
        records(x)
        if (is.null(ready())) status("Import paused: review duplicate feature IDs below.") else {
          generated <- sum(vapply(x,function(r) sum(r$cell_identity$generated),numeric(1)))
          status(paste("Reading data completed:",length(x),"sample(s).",generated,"cell identifiers generated. Source IDs and metadata retained for annotation review."))
        }
      },error=function(e) {
        if (token==revision && !session$isClosed()) status(paste("Import stopped:",conditionMessage(e)))
      })
      if (token==revision) busy(FALSE)
    }),once=TRUE)
  },ignoreInit=TRUE)
  shiny::observeEvent(input$resolveImport,{
    tryCatch({
      x <- records(); sample <- input$importDecisionSample
      if (is.null(x) || length(sample)!=1L || !sample %in% names(x)) stop("Select a pending sample.")
      x[[sample]] <- resolve_import_duplicates(x[[sample]],input$importDecisionAction,input$importDecisionReason)
      records(x)
      status(if (identical(x[[sample]]$status,"cancelled")) "Import cancelled. Upload corrected data and read again." else if (is.null(ready())) "Decision recorded; other samples still need review." else "Import decisions recorded. Reading data completed; all samples are ready.")
    },error=function(e) status(paste("Decision not applied:",conditionMessage(e))))
  },ignoreInit=TRUE)
  list(records=shiny::reactive(records()),ready=ready,status=shiny::reactive(status()),
    busy=shiny::reactive(busy()),decisions=decisions,directory=directory)
}
