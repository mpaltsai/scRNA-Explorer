# Optional diagnostics own tables/configuration, not another collection of SCEs.
# Changing thresholds only subsets those tables; baseline/method-setting changes
# invalidate evidence. Each Run captures settings before deferring computation.
qc_diagnostics_server <- function(input,output,session,baselines,review,app_root=".",services=list()) {
  registry<-shiny::reactiveVal(list());drafts<-shiny::reactiveVal(list());uploads<-shiny::reactiveVal(list())
  marker_reviews<-shiny::reactiveVal(list());history<-shiny::reactiveVal(list())
  generation<-shiny::reactiveVal(0L);visit<-shiny::reactiveVal(0L)
  busy<-shiny::reactiveVal(character());notice<-shiny::reactiveVal("")
  library_path<-qc_optional_library(app_root)
  methods<-c("doublets","ambient","cell_cycle")
  runner<-if(is.null(services$run)) qc_run_optional else services$run
  marker_runner<-if(is.null(services$markers)) qc_cycle_marker_review else services$markers
  baseline_ids<-review$baseline_ids
  shiny::observeEvent(baselines(),{
    ids<-baseline_ids();old<-shiny::isolate(registry());old_drafts<-shiny::isolate(drafts());old_uploads<-shiny::isolate(uploads())
    fresh<-list();d<-list();u<-list()
    for(n in names(ids)) {
      unchanged<-!is.null(old[[n]]) && identical(old[[n]]$baseline_hash,ids[[n]])
      fresh[[n]]<-if(unchanged) old[[n]] else list(baseline_hash=ids[[n]],results=list())
      if(unchanged) {d[[n]]<-old_drafts[[n]];u[[n]]<-old_uploads[[n]]}
    }
    registry(fresh);drafts(d);uploads(u);marker_reviews(list());busy(character())
    generation(shiny::isolate(generation())+1L)
  },ignoreNULL=FALSE,priority=70)
  selected<-review$selected
  raw_fields<-function(fields,p) stats::setNames(lapply(fields,function(n) input[[paste0(p,n)]]),fields)
  previous<-NULL
  shiny::observeEvent(list(selected(),generation()),{
    # Save a departing sample's last input values before changing namespaces,
    # including edits co-delivered with the sample switch. A fresh visit token
    # prevents a queued Run click from targeting a previously visited sample
    # whose old ready fields still exist in Shiny's input cache.
    if(!is.null(previous) && identical(previous$generation,generation())) {
      d<-shiny::isolate(drafts())
      for(k in c("common",methods)) if(identical(input[[paste0(previous$prefix,"ready_",k)]],previous$prefix)) {
        fields<-if(k=="common") .qc_common_controls else .qc_method_controls[[k]]
        d[[previous$sample]][[k]]<-raw_fields(fields,previous$prefix)
      }
      drafts(d)
    }
    v<-shiny::isolate(visit())+1L;visit(v);n<-selected()
    previous<<-if(is.null(n)) NULL else list(sample=n,generation=generation(),
      prefix=paste0("qc_optional_",generation(),"_",v,"_",match(n,names(baselines())),"_"))
  },ignoreNULL=FALSE,priority=60)
  prefix<-shiny::reactive({n<-selected();if(is.null(n)) return(NULL)
    paste0("qc_optional_",generation(),"_",visit(),"_",match(n,names(baselines())),"_")})
  for(kind in c("common",methods)) local({
    k<-kind;fields<-if(k=="common") .qc_common_controls else .qc_method_controls[[k]]
    shiny::observeEvent({p<-prefix();if(is.null(p) || !identical(input[[paste0(p,"ready_",k)]],p)) NULL else raw_fields(fields,p)}, {
      n<-selected();p<-prefix();d<-shiny::isolate(drafts());value<-raw_fields(fields,p)
      if(!identical(d[[n]][[k]],value)) {d[[n]][[k]]<-value;drafts(d)}
    },ignoreNULL=TRUE)
    output[[paste0("qcOptionalControls_",k)]]<-shiny::renderUI({
      p<-prefix();n<-selected();if(is.null(p)) return(shiny::helpText("Create a baseline and select a sample first."))
      # Results/drafts do not rebuild controls on every keystroke. Upload changes
      # rebuild only the corresponding file controls, restoring saved values.
      u<-if(k %in% c("doublets","cell_cycle")) uploads()[[n]] else shiny::isolate(uploads()[[n]])
      qc_diagnostic_widgets(p,k,baselines()[[n]],shiny::isolate(drafts()[[n]][[k]]),library_path,u)
    })
  })
  for(type in c("rates","markers")) local({
    kind<-type
    shiny::observeEvent({p<-prefix();if(is.null(p)) NULL else input[[paste0(p,kind,"_upload")]]},{
      n<-selected();p<-prefix();u<-shiny::isolate(uploads())
      tryCatch({u[[n]][[kind]]<-qc_read_diagnostic_table(input[[paste0(p,kind,"_upload")]]);uploads(u);notice("Diagnostic table stored for this sample.")},
        error=function(e) {u[[n]][[kind]]<-NULL;uploads(u);notice(paste("Upload stopped:",conditionMessage(e)))})
    },ignoreNULL=TRUE)
    shiny::observeEvent(input[[if(kind=="rates") "qcClearRates" else "qcClearMarkers"]],{
      n<-selected();if(is.null(n)) return();u<-shiny::isolate(uploads());u[[n]][[kind]]<-NULL;uploads(u)
      notice("Stored diagnostic table cleared for this sample.")
    },ignoreInit=TRUE,priority=30)
  })
  capture_context<-function(sample) {
    other<-setdiff(names(baselines()),sample);d<-drafts()
    declarations<-stats::setNames(lapply(other,function(n) {
      x<-d[[n]]$common
      if(is.null(x$capture_mode) || identical(x$capture_mode,"unresolved")) return(list(mode="unresolved"))
      list(mode=x$capture_mode,id=x$capture_id,column=x$capture_column,complete=x$capture_complete,evidence=x$capture_evidence)
    }),other)
    list(declarations=declarations,baseline_ids=baseline_ids()[other])
  }
  configs<-stats::setNames(lapply(methods,function(m) {force(m);shiny::reactive({
    n<-selected();p<-prefix();if(is.null(p)) return(NULL)
    saved<-drafts()[[n]]
    common<-if(identical(input[[paste0(p,"ready_common")]],p)) raw_fields(.qc_common_controls,p) else saved$common
    specific<-if(identical(input[[paste0(p,"ready_",m)]],p)) raw_fields(.qc_method_controls[[m]],p) else saved[[m]]
    if(is.null(common) || is.null(specific)) return(NULL)
    proxy<-stats::setNames(c(common,specific),paste0(p,c(names(common),names(specific))))
    proxy[[paste0(p,"ready_common")]]<-p;proxy[[paste0(p,"ready_",m)]]<-p
    cfg<-qc_diagnostic_config(proxy,p,m,uploads()[[n]])
    if(m=="doublets" && !is.null(cfg)) cfg$cohort_context<-capture_context(n)
    cfg
  })}),methods)
  current_results<-function(sample) {
    r<-registry()[[sample]]$results
    for(m in names(r)) {
      cfg<-if(identical(sample,selected())) configs[[m]]() else {
        p<-paste0("stored_");d<-drafts()[[sample]]
        if(is.null(d$common) || is.null(d[[m]])) NULL else {
          proxy<-stats::setNames(c(d$common,d[[m]]),paste0(p,c(names(d$common),names(d[[m]]))))
          proxy[[paste0(p,"ready_common")]]<-p;proxy[[paste0(p,"ready_",m)]]<-p
          cfg<-qc_diagnostic_config(proxy,p,m,uploads()[[sample]])
          if(m=="doublets" && !is.null(cfg)) cfg$cohort_context<-capture_context(sample)
          cfg
        }
      }
      if(is.null(cfg) || !identical(r[[m]]$parameter_hash,digest::digest(cfg,algo="sha256"))) {
        r[[m]]$previous_revision<-r[[m]]$revision;r[[m]]$status<-"stale_settings"
        r[[m]]$reasons<-"Settings changed since this fit. Run the check again; previous scores are not current evidence."
        for(field in setdiff(names(r[[m]]$values),c("cell_id","cells"))) r[[m]]$values[[field]][]<-NA
        r[[m]]$revision<-digest::digest(list(r[[m]]$previous_revision,cfg,status="stale_settings"),algo="sha256")
      }
    }
    r
  }
  output$qcOptionalNotice<-shiny::renderText(notice())
  output$qcOptionalBusy<-shiny::renderUI({if(length(busy())) shiny::tags$div(role="status",shiny::icon("spinner",class="fa-spin"),paste("Running optional QC:",paste(busy(),collapse=", ")) )})
  output$qcOptionalPackages<-shiny::renderTable({
    p<-prefix();if(is.null(p)) return(NULL)
    path<-input[[paste0(p,"library_path")]]
    libs<-if(isTRUE(input[[paste0(p,"use_library")]]) && length(path)==1L && dir.exists(path)) c(path,.libPaths()) else .libPaths()
    do.call(rbind,lapply(c("scDblFinder","decontX","Seurat","SoupX"),function(package) {
      found<-find.package(package,lib.loc=libs,quiet=TRUE)
      data.frame(package=package,available=length(found)>0L,version=if(isNamespaceLoaded(package)) as.character(getNamespaceVersion(package)) else if(length(found)) read.dcf(file.path(found,"DESCRIPTION"),fields="Version")[1] else NA_character_)
    }))
  },rownames=FALSE)
  output$qcOptionalStatus<-shiny::renderTable({
    n<-selected();if(is.null(n)) return(NULL);r<-current_results(n)
    do.call(rbind,lapply(methods,function(m) {
      x<-r[[m]]
      data.frame(check=m,status=if(paste(n,m,sep="|") %in% busy()) "running" else if(is.null(x)) "not run" else x$status,
        notes=if(is.null(x)) "No scores assigned. Enable/review settings and run this optional check." else paste(c(x$reasons,x$warnings,x$diagnostics$interpretation),collapse=" | "))
    }))
  },rownames=FALSE)
  live_controls<-function(method,p) identical(input[[paste0(p,"ready_common")]],p) && identical(input[[paste0(p,"ready_",method)]],p)
  queue_run<-function(method,markers_only=FALSE) {
    n<-selected();p<-prefix();if(is.null(n) || !live_controls(method,p)) {notice("Wait for the selected sample's diagnostic controls, then run the check.");return()}
    cfg<-configs[[method]]();b<-baselines()[[n]];hash<-baseline_ids()[[n]]
    others<-baselines()[setdiff(names(baselines()),n)]
    # Known declarations in other samples affect eligibility, including when
    # their capture metadata columns have different names. Unresolved inputs
    # are retained as unresolved; completeness remains a user declaration.
    if(method=="doublets") for(other in names(cfg$cohort_context$declarations)) {
      declaration<-cfg$cohort_context$declarations[[other]];values<-character()
      if(identical(declaration$mode,"single") && length(declaration$id)==1L && nzchar(declaration$id)) values<-trimws(declaration$id)
      if(identical(declaration$mode,"metadata") && length(declaration$column)==1L && declaration$column %in% names(SummarizedExperiment::colData(baselines()[[other]])))
        values<-as.character(baselines()[[other]][[declaration$column]])
      values<-values[!is.na(values)&nzchar(values)]
      if(length(values)) others<-c(others,list(list(declared_captures=values)))
    }
    key<-paste(n,method,sep="|");if(key %in% busy()) return()
    busy(c(busy(),key));notice(paste("Preparing",method,"on the fixed baseline. No filtering or count correction will be applied."))
    session$onFlushed(function() shiny::isolate({
      on.exit(busy(setdiff(busy(),key)),add=TRUE)
      if(session$isClosed() || !identical(baseline_ids()[[n]],hash)) return(invisible(NULL))
      tryCatch({
        if(method=="cell_cycle") {
          checked<-tryCatch(marker_runner(b,cfg,app_root),error=function(e) list(error=conditionMessage(e)))
          mr<-marker_reviews();mr[[n]]<-list(parameter_hash=digest::digest(cfg,algo="sha256"),review=checked);marker_reviews(mr)
          if(markers_only) {notice(if(!is.null(checked$error)) checked$error else "Marker coverage reviewed. No scores assigned by this action.");return(invisible(NULL))}
        }
        result<-runner(b,method,cfg,hash,app_root,others)
        if(!identical(baseline_ids()[[n]],hash)) return(invisible(NULL))
        store<-registry();store[[n]]$results[[method]]<-result;registry(store)
        h<-history();h[[n]]<-c(h[[n]],list(result[setdiff(names(result),"values")]));history(h)
        notice(paste(method,":",result$status,"—",paste(c(result$reasons,result$warnings),collapse=" | ")))
      },error=function(e) notice(paste("Optional check stopped:",conditionMessage(e))))
    }),once=TRUE)
  }
  for(m in methods) local({method<-m;shiny::observeEvent(input[[paste0("qcRun_",method)]],{queue_run(method)},ignoreInit=TRUE)})
  shiny::observeEvent(input$qcReviewCycleMarkers,{queue_run("cell_cycle",TRUE)},ignoreInit=TRUE)
  fitted<-function(method) {
    n<-selected();shiny::validate(shiny::need(!is.null(n),"Create a QC baseline and select a sample."))
    r<-current_results(n)[[method]]
    shiny::validate(shiny::need(!paste(n,method,sep="|") %in% busy(),"This optional check is running."),
      shiny::need(!is.null(r),"Review settings and run this optional check first."))
    shiny::validate(shiny::need(r$status %in% c("scored","estimated"),paste(c(r$status,r$reasons),collapse=" — ")))
    r
  }
  compared<-function(method) qc_diagnostic_comparison(review$comparison(),fitted(method))
  show<-function(method) shiny::req(identical(input$qcReviewTab,"optional"),identical(input$qcOptionalTab,method))
  x_axis<-function() if(identical(input$qcOptionalX,"nGene")) "nGene" else "nUMI"
  plots<-list(qcDoubletDistribution=function() qc_plot_doublet_distribution(compared("doublets")),
    qcDoubletJoint=function() qc_plot_doublet_joint(compared("doublets"),x_axis()),qcDoubletCounts=function() qc_plot_doublet_counts(compared("doublets")),
    qcAmbientDistribution=function() qc_plot_ambient_distribution(compared("ambient")),qcAmbientJoint=function() qc_plot_ambient_joint(compared("ambient"),x_axis()),
    qcAmbientPopulations=function() qc_plot_ambient_populations(compared("ambient")),qcAmbientFit=function() {
      r<-fitted("ambient");tryCatch(qc_plot_ambient_fit(r),error=function(e) shiny::validate(shiny::need(FALSE,conditionMessage(e))))
    },qcCycleScores=function() qc_plot_cycle_scores(compared("cell_cycle")),qcCycleCounts=function() qc_plot_cycle_counts(compared("cell_cycle")),qcCycleQC=function() qc_plot_cycle_qc(compared("cell_cycle")))
  for(id in names(plots)) local({
    name<-id;make<-plots[[id]];method<-if(grepl("Doublet",id)) "doublets" else if(grepl("Ambient",id)) "ambient" else "cell_cycle"
    output[[name]]<-shiny::renderPlot({show(method);make()})
  })
  output$qcAmbientFits<-shiny::renderTable({qc_ambient_fit_table(fitted("ambient"))},rownames=FALSE)
  marker_review<-shiny::reactive({
    n<-selected();shiny::validate(shiny::need(!is.null(n),"Create a baseline and select a sample first."))
    entry<-marker_reviews()[[n]];cfg<-configs$cell_cycle()
    shiny::validate(shiny::need(!is.null(entry) && !is.null(cfg),"Review marker coverage or run cell-cycle scoring first."))
    shiny::validate(shiny::need(identical(entry$parameter_hash,digest::digest(cfg,algo="sha256")),"Marker settings changed. Review coverage again."))
    shiny::validate(shiny::need(is.null(entry$review$error),entry$review$error));entry$review
  })
  output$qcCycleCoverage<-shiny::renderTable(marker_review()$coverage,rownames=FALSE)
  output$qcCycleMapping<-DT::renderDT(DT::datatable(marker_review()$mapping,rownames=FALSE,options=list(scrollX=TRUE,pageLength=10)))
  cell_table<-shiny::reactive({
    n<-selected();shiny::req(n);d<-review$metrics();r<-current_results(n)
    chosen<-review$selection()
    d$qc_selection<-if(is.null(chosen)) "baseline only" else ifelse(chosen$keep,"retained",
      if(identical(input$qcView,"applied")) "excluded (applied)" else "proposed exclusion")
    d$qc_filtering_reasons<-if(is.null(chosen)) "" else chosen$per_cell$reasons
    for(m in names(r)) {
      i<-match(d$cell_id,r[[m]]$values$cell_id);if(anyNA(i)) stop("Diagnostic cell identifiers differ from baseline.")
      for(field in setdiff(names(r[[m]]$values),c("cell_id","cells"))) d[[field]]<-r[[m]]$values[[field]][i]
      d[[paste0("qc_diag_",m,"_status")]]<-r[[m]]$status
      d[[paste0("qc_diag_",m,"_revision")]]<-r[[m]]$revision
    }
    d
  })
  output$qcOptionalCells<-DT::renderDT(DT::datatable(cell_table(),rownames=FALSE,options=list(scrollX=TRUE,pageLength=10)))
  output$qcOptionalCSV<-shiny::downloadHandler(filename=function() paste0(selected(),"-optional-qc.csv"),content=function(path) utils::write.csv(cell_table(),path,row.names=FALSE))
  output$qcOptionalEvidence<-shiny::downloadHandler(filename=function() "optional-qc-evidence.rds",content=function(path) {
    samples<-names(registry());results<-stats::setNames(lapply(samples,current_results),samples)
    saveRDS(list(results=results,history=history()),path)
  })
  output$qcCycleMappingDownload<-shiny::downloadHandler(filename=function() "cell-cycle-marker-mapping.csv",content=function(path) utils::write.csv(marker_review()$mapping,path,row.names=FALSE))
  list(results=current_results,history=shiny::reactive(history()),prefix=prefix,configs=configs,marker_review=marker_review,
    compared=compared,cell_table=cell_table,busy=shiny::reactive(busy()),
    decorate=function(object,sample) qc_attach_optional(object,current_results(sample),baseline_ids()[[sample]]))
}
