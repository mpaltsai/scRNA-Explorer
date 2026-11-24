# Each sample owns a draft, a compact preview, one applied SCE and a decision
# history. Control IDs include the cohort generation and sample index, preventing
# another sample's numeric/slider values from being applied during a UI switch.
qc_review_server <- function(input,output,session,baselines) {
  entries <- shiny::reactiveVal(list())
  baseline_metrics <- shiny::reactiveVal(NULL)
  generation <- shiny::reactiveVal(0L)
  shiny::observeEvent(baselines(),{
    b <- baselines(); old <- shiny::isolate(entries()); fresh <- list()
    if (!is.null(b)) for (name in names(b)) {
      hash <- qc_baseline_signature(b[[name]])
      prior <- old[[name]]
      fresh[[name]] <- if (!is.null(prior) && identical(prior$baseline_hash,hash)) prior else list(
        baseline_hash=hash,metrics=qc_review_metrics(b[[name]]),
        draft=if(ncol(b[[name]])) qc_default_thresholds(baseline_qc_metrics(b[[name]])) else NULL,
        preview=NULL,filtered=NULL,applied=NULL,history=list(),message="Preview thresholds before applying filtering.")
    }
    entries(fresh);baseline_metrics(lapply(fresh,`[[`,"metrics"));generation(shiny::isolate(generation())+1L)
    choices <- names(fresh)[vapply(fresh,function(e) nrow(e$metrics)>0L,logical(1))]
    previous <- shiny::isolate(input$selectSample)
    shiny::updateSelectInput(session,"selectSample",choices=choices,
      selected=if(length(previous)==1L && previous %in% choices) previous else if(length(choices)) choices[1] else character())
  },ignoreNULL=FALSE,priority=80)
  selected <- shiny::reactive({
    name <- input$selectSample; b <- baselines()
    if (length(name)!=1L || !name %in% names(b) || !ncol(b[[name]])) return(NULL)
    name
  })
  current <- shiny::reactive({n<-selected();if(is.null(n)) NULL else entries()[[n]]})
  control_ids <- shiny::reactive({
    name <- selected();if(is.null(name)) return(NULL)
    prefix <- paste0("qc_",generation(),"_",match(name,names(baselines())),"_")
    stats::setNames(paste0(prefix,c(qc_rule_spec()$input,"nUMIs","nGenes","mitoRatio","rbcRatio","novelty")),
      c(qc_rule_spec()$input,"nUMIs","nGenes","mitoRatio","rbcRatio","novelty"))
  })
  draft <- shiny::reactive({
    ids <- control_ids();if(is.null(ids)) return(NULL)
    values <- lapply(ids[qc_rule_spec()$input],function(id) input[[id]])
    if(any(vapply(values,is.null,logical(1)))) return(NULL)
    values
  })
  shiny::observeEvent(draft(),{
    d <- draft(); n <- selected(); if(is.null(d) || is.null(n)) return()
    e <- shiny::isolate(entries())
    if (!identical(e[[n]]$draft,d)) {e[[n]]$draft<-d;entries(e)}
  },ignoreNULL=TRUE)
  output$qcReviewControls <- shiny::renderUI({
    ids <- control_ids();if(is.null(ids)) return(shiny::helpText("Create an annotated QC baseline first."))
    e <- shiny::isolate(current()); d <- e$draft;m <- e$metrics
    spec <- qc_rule_spec()
    controls <- lapply(seq_len(nrow(spec)),function(i) shiny::numericInput(ids[[spec$input[i]]],
      paste(spec$label[i],if(i>=6) "(fraction 0–1)" else "cutoff"),value=d[[spec$input[i]]],
      min=if(i==1) 100 else 0,max=if(i>=5) 1 else NA,step=if(i<=4) 1 else 0.0001))
    # Exact numeric values are authoritative. As in the original app, users
    # explicitly copy slider values into the numeric draft. This avoids rounding
    # precise numeric cutoffs to a slider's step without the user's action.
    sliders <- shiny::tags$details(shiny::tags$summary("Slider controls"),
      shiny::helpText("Move sliders, then Copy slider values to draft. Preview uses the numeric cutoffs above; dashed plot lines use the preview/applied snapshot."),
      shiny::sliderInput(ids[["nUMIs"]],"UMIs",min=100,max=max(101,m$nUMI,d$maxUMIsPerCell),value=c(d$minUMIsPerCell,d$maxUMIsPerCell),step=1),
      shiny::sliderInput(ids[["nGenes"]],"Detected genes",min=0,max=max(1,m$nGene,d$maxGenesPerCell),value=c(d$minGenesPerCell,d$maxGenesPerCell),step=1),
      shiny::sliderInput(ids[["mitoRatio"]],"Mitochondrial fraction",min=0,max=1,value=c(d$minMitoPerCell,d$maxMitoPerCell),step=0.0001),
      shiny::sliderInput(ids[["rbcRatio"]],"Hemoglobin fraction",min=0,max=1,value=c(d$minRbcPerCell,d$maxRbcPerCell),step=0.0001),
      shiny::sliderInput(ids[["novelty"]],"log10(genes) / log10(UMIs)",min=0,max=1,value=d$noveltyPerCell,step=0.0001),
      shiny::actionButton("updateThresholds","Copy slider values to draft"))
    shiny::tagList(controls,sliders)
  })
  shiny::observeEvent(input$updateThresholds,{
    ids <- control_ids();if(is.null(ids)) return()
    names <- c("nUMIs","nGenes","mitoRatio","rbcRatio","novelty")
    s <- lapply(ids[names],function(id) input[[id]])
    if (!identical(lengths(s),c(nUMIs=2L,nGenes=2L,mitoRatio=2L,rbcRatio=2L,novelty=1L))) return()
    values <- c(s$nUMIs,s$nGenes,s$novelty,s$mitoRatio,s$rbcRatio)
    for(i in seq_along(values)) shiny::updateNumericInput(session,ids[[qc_rule_spec()$input[i]]],value=values[i])
  },ignoreInit=TRUE)
  preview <- shiny::reactive({
    e <- current(); d <- draft()
    if(is.null(e) || is.null(e$preview) || is.null(d)) return(NULL)
    valid <- tryCatch(qc_validate_thresholds(d),error=function(e) NULL)
    if (!identical(valid,e$preview$thresholds) || !identical(e$baseline_hash,e$preview$baseline_hash)) return(NULL)
    e$preview
  })
  # A browser can deliver Preview, Apply and plot clicks in one input batch.
  # Order these dependent actions explicitly: Preview (20), Apply (10), plot
  # request (0). Otherwise a request can capture the previous applied revision
  # and be rejected by the stale-result guard without ever calling scater.
  shiny::observeEvent(input$qcPreview,{
    n <- selected();if(is.null(n)) return()
    e <- shiny::isolate(entries())
    tryCatch({
      if(is.null(draft())) stop("Complete the current sample's numeric controls first.")
      e[[n]]$preview <- qc_filter_proposal(e[[n]]$metrics,draft(),e[[n]]$baseline_hash)
      e[[n]]$message <- "Preview ready. No filtering has been applied by this action."
    },error=function(err) {e[[n]]$preview<<-NULL;e[[n]]$message<<-paste("Preview stopped:",conditionMessage(err))})
    entries(e)
  },ignoreInit=TRUE,priority=20)
  shiny::observeEvent(input$QCFiltering,{
    n <- selected();if(is.null(n)) return()
    e <- shiny::isolate(entries())
    tryCatch({
      p <- preview();if(is.null(p)) stop("Preview the current thresholds for this sample before applying.")
      result <- qc_apply_filter(baselines()[[n]],p,input$qcApplyReason,e[[n]]$history)
      e[[n]]$filtered <- result$object;e[[n]]$applied<-result$decision;e[[n]]$history<-result$history
      e[[n]]$message <- paste("Applied",sum(p$keep),"of",length(p$keep),"baseline cells. This sample's plots, download and downstream input use the recorded revision.")
    },error=function(err) {e[[n]]$message<<-paste("Apply stopped:",conditionMessage(err))})
    entries(e)
  },ignoreInit=TRUE,priority=10)
  selection <- shiny::reactive({
    e <- current();if(is.null(e)) return(NULL)
    if (!identical(input$qcView,"applied")) return(preview())
    if(is.null(e$applied)) return(NULL)
    p <- e$applied
    list(keep=p$per_cell$decision=="retain",per_cell=p$per_cell,thresholds=p$thresholds,
      summary=p$summary,rule_counts=p$rule_counts,overlap=p$overlap)
  })
  comparison <- shiny::reactive({
    e <- current();shiny::req(e)
    qc_comparison_metrics(e$metrics,selection(),if(identical(input$qcView,"applied")) "Applied retained" else "Preview retained")
  })
  output$QCDataFilt <- shiny::renderText({e<-current();if(is.null(e)) "Select a baseline sample." else e$message})
  output$qcComparisonStatus <- shiny::renderText({
    e <- current();if(is.null(e)) return("Create an annotated QC baseline, then select a sample.")
    p <- selection()
    paste(selected(),"—",nrow(e$metrics),"baseline cells.",
      if(is.null(p)) "No current selection for this view; only baseline is shown. Preview again after editing cutoffs." else
        paste(if(identical(input$qcView,"applied")) "Applied:" else "Preview:",sum(p$keep),"retained;",sum(!p$keep),"excluded."),
      if(is.null(e$applied)) "No applied filtering decision." else paste("Stored applied revision:",substr(e$applied$revision,1,12),".",ncol(e$filtered),"cells."))
  })
  output$qcSampleSummary <- shiny::renderTable({
    e<-entries();if(!length(e)) return(NULL)
    do.call(rbind,lapply(names(e),function(n) data.frame(sample=n,baseline_cells=nrow(e[[n]]$metrics),
      applied_cells=if(is.null(e[[n]]$filtered)) NA_integer_ else ncol(e[[n]]$filtered),
      median_UMIs=if(nrow(e[[n]]$metrics)) stats::median(e[[n]]$metrics$nUMI) else NA_real_,
      unavailable_mito_cells=sum(!is.finite(e[[n]]$metrics$mitoRatio)),
      unavailable_hemoglobin_cells=sum(!is.finite(e[[n]]$metrics$rbcRatio)))))
  },rownames=FALSE)
  # Compute graphics only for the chosen panel. Tables/masks are independent
  # of panel visibility; hidden plots do not consume time on threshold edits.
  visible <- function(tab) shiny::req(identical(input$qcReviewTab,tab))
  output$qcCohortOverview <- shiny::renderPlot({
    visible("overview")
    m<-baseline_metrics();shiny::req(length(m)>0)
    qc_plot_overview(do.call(rbind,m))
  })
  output$qcSelectedOverview <- shiny::renderPlot({visible("overview");qc_plot_overview(comparison(),compare=TRUE)})
  joint <- shiny::reactive(qc_joint_plot_data(comparison(),if(is.null(input$qcJointX)) "nGene" else input$qcJointX))
  output$qcJointPlot <- shiny::renderPlot({visible("joint");qc_plot_joint(joint(),
    if(is.null(selection())) NULL else selection()$thresholds,if(is.null(input$qcJointX)) "nGene" else input$qcJointX)})
  output$qcHoveredCells <- shiny::renderTable({
    shiny::req(input$qcJointHover)
    d<-shiny::nearPoints(joint(),input$qcJointHover,xvar="qc_x",yvar="qc_y",maxpoints=5,threshold=10)
    d[,c("cell_id","cells","sample","nUMI","nGene","mitoRatio","selection","reasons"),drop=FALSE]
  },rownames=FALSE)
  output$qcRibosomalPlot <- shiny::renderPlot({visible("ribosomal");qc_plot_ribosomal(comparison())})
  output$qcRibosomalCoverage <- shiny::renderTable({e<-current();shiny::req(e);qc_ribosomal_coverage(e$metrics)},rownames=FALSE)
  output$qcAnnotationCoverage <- shiny::renderPlot({visible("annotation");qc_plot_annotation_coverage(comparison())})
  output$qcRuleCounts <- shiny::renderTable({p<-selection();shiny::req(p);p$rule_counts},rownames=FALSE)
  output$qcRuleImpact <- shiny::renderPlot({visible("impact");p<-selection();shiny::req(p);qc_plot_rule_impact(p)})
  output$qcRuleOverlap <- shiny::renderPlot({
    visible("impact")
    p<-selection();shiny::req(p)
    qc_plot_overlap(p)
  })
  sensitivity <- shiny::reactive({
    p<-selection();e<-current();shiny::req(p,e)
    d<-qc_threshold_sensitivity(e$metrics,p$thresholds)
    d$view<-if(identical(input$qcView,"applied")) "applied" else "preview"
    d$baseline_hash<-e$baseline_hash
    d$selection_id<-if(identical(input$qcView,"applied")) e$applied$revision else p$proposal_id
    d
  })
  output$qcSensitivity <- shiny::renderPlot({visible("sensitivity");qc_plot_sensitivity(sensitivity())})
  cell_table <- shiny::reactive({
    e<-current();shiny::req(e);p<-selection();d<-e$metrics
    d$decision<-if(is.null(p)) "not previewed" else ifelse(p$keep,"retain",if(identical(input$qcView,"applied")) "excluded (applied)" else "proposed exclusion")
    d$reasons<-if(is.null(p)) "" else p$per_cell$reasons
    d
  })
  output$qcCells <- DT::renderDT(DT::datatable(cell_table(),rownames=FALSE,options=list(scrollX=TRUE,pageLength=10),selection="single"))
  output$qcHistory <- shiny::renderTable({
    e<-current();if(is.null(e) || !length(e$history)) return(NULL)
    do.call(rbind,lapply(e$history,function(d) data.frame(revision=substr(d$revision,1,12),applied_at=d$applied_at,
      reason=d$reason,retained=d$summary$retained,excluded=d$summary$excluded)))
  },rownames=FALSE)
  output$qcCellsDownload <- shiny::downloadHandler(filename=function() paste0(selected(),"-qc-review.csv"),
    content=function(path) utils::write.csv(cell_table(),path,row.names=FALSE))
  output$qcSensitivityDownload <- shiny::downloadHandler(filename=function() paste0(selected(),"-threshold-sensitivity.csv"),
    content=function(path) utils::write.csv(sensitivity(),path,row.names=FALSE))
  output$qcDecisionsDownload <- shiny::downloadHandler(filename=function() "qc-filtering-decisions.rds",
    content=function(path) saveRDS(lapply(entries(),`[[`,"history"),path))
  output$QCDownload <- shiny::downloadHandler(filename=function() paste0(selected(),"-applied-qc.rds"),content=function(path) {
    e<-current();shiny::req(e$filtered);saveRDS(e$filtered,path)
  })
  # Preserve original plotting routines in an additional review panel. Pairing
  # via a single faceted plot gives both stages identical axes/colour scales.
  legacy <- list(UMIsPlot=function(d,t) createUMIsPlot(d,t$minUMIsPerCell,t$maxUMIsPerCell),
    GenesHistPlot=function(d,t) createGenesHistPlot(d,t$minGenesPerCell,t$maxGenesPerCell),
    GenesBoxPlot=function(d,t) createGenesBoxPlot(d),
    UMIsGenesPlot=function(d,t) createUMIsGenesPlot(d,t$minGenesPerCell,t$maxGenesPerCell),
    UMIsGenesPlotRBC=function(d,t) createUMIsGenesPlotRBC(d,t$minGenesPerCell,t$maxGenesPerCell),
    mitoRatioPlot=function(d,t) createmitoRatioPlot(d,t$minMitoPerCell,t$maxMitoPerCell),
    rbcRatioPlot=function(d,t) createrbcRatioPlot(d,t$minRbcPerCell,t$maxRbcPerCell),
    noveltyPlot=function(d,t) createnoveltyPlot(d,t$noveltyPerCell,1))
  for(id in names(legacy)) local({
    name<-id;make<-legacy[[id]]
    output[[name]] <- shiny::renderPlot({
      visible("existing")
      d<-comparison();e<-current();t<-if(is.null(selection())) e$draft else selection()$thresholds
      make(d,t)+ggplot2::facet_wrap(~stage,scales="fixed")
    })
  })
  # Heavy expression plots are requested explicitly. A request belongs to one
  # sample/revision; changing that identity clears the view without computing a
  # new expression plot or constructing a proposed SCE.
  expression_plot <- function(button,id,applied,make) {
    requested <- shiny::reactiveVal(NULL)
    identity <- shiny::reactive({e<-current();if(is.null(e)) return(NULL)
      list(sample=selected(),baseline=e$baseline_hash,revision=if(applied && !is.null(e$applied)) e$applied$revision else "baseline")})
    shiny::observeEvent(input[[button]],{requested(identity())},ignoreInit=TRUE,priority=0)
    output[[id]] <- shiny::renderPlot({
      visible("existing")
      e<-current()
      shiny::validate(shiny::need(!is.null(e),"Select a sample with a QC baseline first."))
      object <- if(applied) current()$filtered else baselines()[[selected()]]
      # Validation messages are visible below the button. Silent req() checks
      # made missing/empty/stale requests indistinguishable from a broken plot.
      shiny::validate(shiny::need(!is.null(object),"Apply filtering to this sample before requesting an applied plot."))
      shiny::validate(shiny::need(ncol(object)>0L,"No cells remain in the applied result. Revise thresholds, preview and apply again."))
      shiny::validate(shiny::need(!is.null(requested()),"Click the plot button to generate this view."))
      shiny::validate(shiny::need(identical(requested(),identity()),
        "The sample or applied filtering changed. Click the plot button again for the current result."))
      make(object)
    })
  }
  expression_plot("makeLibrary","libraryPlot",FALSE,createlibraryPlot)
  expression_plot("makeLibraryFiltered","libraryPlotFiltered",TRUE,createlibraryPlot)
  expression_plot("makeHighExpr","highExprPlot",FALSE,createhighExprPlot)
  expression_plot("makeHighExprFiltered","highExprPlotFiltered",TRUE,createhighExprPlot)
  filtered <- shiny::reactive({e<-entries();if(!length(e)) NULL else lapply(e,`[[`,"filtered")})
  revisions <- shiny::reactive(vapply(entries(),function(e) if(is.null(e$applied)) "" else e$applied$revision,character(1)))
  list(selected=selected,metrics=shiny::reactive({e<-current();if(is.null(e)) NULL else e$metrics}),
    draft=draft,preview=preview,selection=selection,comparison=comparison,sensitivity=sensitivity,control_ids=control_ids,
    filtered=filtered,revisions=revisions,history=shiny::reactive(lapply(entries(),`[[`,"history")),
    require_ready=shiny::reactive(qc_require_applied(baselines(),filtered())))
}
