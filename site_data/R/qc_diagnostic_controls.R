# File inputs are small metadata/marker tables, not another expression import.
qc_read_diagnostic_table <- function(upload) {
  if(!is.data.frame(upload) || nrow(upload)!=1L || !all(c("name","datapath") %in% names(upload))) stop("Upload one CSV or TSV table.")
  p<-upload$datapath[1]
  if(!file.exists(p) || file.info(p)$size>5*1024^2) stop("Diagnostic tables must exist and be at most 5 MB.")
  ext<-tolower(tools::file_ext(upload$name[1]))
  if(!ext %in% c("csv","tsv","txt")) stop("Use CSV or tab-separated TSV/TXT.")
  d<-utils::read.table(p,header=TRUE,sep=if(ext=="csv") "," else "\t",quote="\"",comment.char="",check.names=FALSE,stringsAsFactors=FALSE,colClasses="character")
  if(anyDuplicated(names(d))) stop("Table columns must have unique names.")
  # Capture IDs and symbols stay literal (including leading zeroes). Only
  # explicitly numeric rate columns are parsed, with invalid text rejected.
  for(n in intersect(c("expected_rate","uncertainty"),names(d))) {
    parsed<-suppressWarnings(as.numeric(d[[n]]))
    if(any(!is.na(d[[n]]) & nzchar(trimws(d[[n]])) & is.na(parsed))) stop("Rate/uncertainty columns must be numeric.")
    d[[n]]<-parsed
  }
  list(table=d,source=list(name=upload$name[1],md5=unname(tools::md5sum(p))))
}

.qc_common_controls<-c("seed","use_library","library_path","capture_mode","capture_id","capture_column","capture_evidence","capture_complete","technology")
.qc_method_controls<-list(doublets=c("doublet_mode","provisional_group_mode","provisional_group_column","rates_mode","expected_rate","uncertainty","rate_source"),
  ambient=c("ambient_method","pooled","clusters_column","maxIter","delta1","delta2","estimateDelta","convergence","varGenes","dbscanEps","droplet_compatible"),
  cell_cycle=c("marker_mode","marker_name","marker_species","nbin","ctrl_auto","ctrl","scale_factor","minimum_markers"))

qc_diagnostic_config <- function(input,prefix,method,uploads) {
  v<-function(n) input[[paste0(prefix,n)]]
  text<-function(n) {x<-v(n);if(is.null(x)) "" else as.character(x)}
  number<-function(n) {x<-v(n);if(is.null(x)) NA_real_ else as.numeric(x)}
  if(!identical(v("ready_common"),prefix) || !identical(v(paste0("ready_",method)),prefix)) return(NULL)
  cfg<-list(seed=number("seed"),library_path=if(isTRUE(v("use_library"))) text("library_path") else "",parameters=list())
  if(isTRUE(v("use_library")) && !nzchar(cfg$library_path)) cfg$invalid_reason<-"Provide an existing isolated library path or disable that option. No packages will be installed."
  if(method %in% c("doublets","ambient")) {
    cfg$capture<-list(mode=text("capture_mode"),id=text("capture_id"),column=text("capture_column"),
      complete=isTRUE(v("capture_complete")),evidence=text("capture_evidence"))
    cfg$technology<-text("technology")
  }
  if(method=="doublets") {
    cfg$analysis_mode<-if(nzchar(text("doublet_mode"))) text("doublet_mode") else "capture"
    if(cfg$analysis_mode=="exploratory") {
      # Ignore unused capture declarations; provisional groups never establish
      # physical capture identity, even if an earlier form contained an ID.
      cfg$capture<-NULL
      cfg$group_mode<-text("provisional_group_mode");cfg$group_column<-text("provisional_group_column")
    }
    cfg$rates_mode<-text("rates_mode")
    if(cfg$rates_mode!="unknown") {
      cfg$expected_rate<-number("expected_rate");cfg$uncertainty<-number("uncertainty");cfg$rate_source<-text("rate_source")
    }
    if(cfg$rates_mode=="csv") {
      if(is.null(uploads$rates)) cfg$invalid_reason<-"Upload a per-group rate table first (the legacy capture column is also accepted)." else {
        cfg$rate_table<-uploads$rates$table;cfg$rate_table_source<-uploads$rates$source
      }
    }
  } else if(method=="ambient") {
    cfg$ambient_method<-text("ambient_method");cfg$exploratory_pooled<-isTRUE(v("pooled"))
    cfg$clusters_column<-text("clusters_column");cfg$droplet_compatible<-isTRUE(v("droplet_compatible"))
    cfg$parameters$model_args<-list(maxIter=number("maxIter"),delta=c(number("delta1"),number("delta2")),
      estimateDelta=isTRUE(v("estimateDelta")),convergence=number("convergence"),varGenes=number("varGenes"),dbscanEps=number("dbscanEps"))
  } else {
    cfg$minimum_markers<-number("minimum_markers")
    cfg$parameters<-list(nbin=number("nbin"),ctrl=if(isTRUE(v("ctrl_auto"))) NULL else number("ctrl"),scale_factor=number("scale_factor"))
    if(text("marker_mode")=="custom") {
      u<-uploads$markers
      if(is.null(u) || !all(c("program","gene_symbol") %in% names(u$table)) ||
        !nrow(u$table) || anyNA(u$table[,c("program","gene_symbol")]) ||
        !setequal(unique(u$table$program),c("S","G2M")) || any(!nzchar(trimws(u$table$gene_symbol))) ||
        !nzchar(text("marker_name")) || !nzchar(text("marker_species"))) cfg$invalid_reason<-"Provide a named species-matched marker CSV/TSV with program (S/G2M) and gene_symbol columns." else {
        cfg$parameters$reference<-list(name=text("marker_name"),species=text("marker_species"),
          s_genes=unique(u$table$gene_symbol[u$table$program=="S"]),g2m_genes=unique(u$table$gene_symbol[u$table$program=="G2M"]),source=u$source)
      }
    }
  }
  cfg
}

qc_diagnostic_widgets <- function(prefix,kind,object,saved=list(),library_path="",uploads=list()) {
  id<-function(n) paste0(prefix,n)
  g<-function(n,default) {
    if(!n %in% names(saved)) return(default)
    value<-saved[[n]]
    if(!is.null(value)) value else if(is.numeric(default)) NA_real_ else if(is.logical(default)) FALSE else ""
  }
  number<-function(n,label,value,min=0,max=NA,step=1) shiny::numericInput(id(n),label,g(n,value),min=min,max=max,step=step)
  text<-function(n,label,value="") shiny::textInput(id(n),label,g(n,value))
  check<-function(n,label,value=FALSE) shiny::checkboxInput(id(n),label,g(n,value))
  select<-function(n,label,choices,value="") shiny::selectInput(id(n),label,choices=choices,selected=g(n,value))
  cols<-names(SummarizedExperiment::colData(object));cols<-cols[!grepl("^qc_|^cell_cycle_|^ambient_|^scDblFinder",cols)]
  ready<-shiny::tags$div(style="display:none",shiny::textInput(id(paste0("ready_",kind)),NULL,prefix))
  # Capture-informed declarations remain visible; public-data exploration
  # explicitly leaves these unresolved and selects its groups below.
  if(kind=="common") return(shiny::tagList(ready,
    shiny::tags$section(shiny::h4("Capture and runtime settings (per sample)"),
      shiny::helpText("A capture is a physical cell-capture experiment. A donor, input file or GEO accession does not establish capture identity."),
      shiny::helpText("For public data with unresolved captures, leave Capture assignment unresolved and select Public-data exploratory in the Doublets panel below."),
      select("capture_mode","Capture assignment",c("Unresolved"="unresolved","Explicitly one complete capture"="single","Existing cell metadata column"="metadata"),"unresolved"),
      text("capture_id","Capture ID for the explicit single-capture declaration"),
      select("capture_column","Capture column",c("Select a column"="",cols)),
      text("capture_evidence","Evidence supporting the capture assignment"),
      check("capture_complete","I confirm this input contains the complete baseline population for each declared capture"),
      text("technology","Protocol / chemistry (if known)"),number("seed","Diagnostic random seed",1,max=.Machine$integer.max),
      check("use_library","Use an existing isolated R library for optional checks",nzchar(library_path)),
      text("library_path","Existing R library path",library_path),
      shiny::helpText("No packages are installed or upgraded. Library search paths and RNG state are restored after each run."))))
  if(kind=="doublets") return(shiny::tagList(ready,
    select("doublet_mode","Doublet analysis mode",c("Capture-informed: use reviewed capture declarations above"="capture","Public-data exploratory: capture identity remains unknown"="exploratory"),"capture"),
    shiny::helpText("Capture-informed mode uses the declarations above. Exploratory mode does not require a capture ID, evidence text or completeness confirmation. Its groups and predictions are provisional. Neither mode removes cells."),
    shiny::conditionalPanel(condition=paste0("input['",id("doublet_mode"),"'] === 'exploratory'"),
      select("provisional_group_mode","Provisional analysis groups",c("Selected input matrix as one provisional group"="input","Groups from an existing cell-metadata column"="metadata"),"input"),
      select("provisional_group_column","Provisional grouping column (metadata option only)",c("Select a column"="",cols)),
      shiny::helpText("Prefer original sample/library groups when available. Groups may contain multiple or partial captures. A donor, condition or cell type is not a confirmed capture. Use the broad input population before cell-type subsetting; review limitations before interpreting predictions.")),
    select("rates_mode","Expected doublet rates",c("Use the explicit rate below for every analysis group"="shared","Upload one rate per group"="csv","Unknown expected rate: classification-error-based thresholding"="unknown"),"shared"),
    shiny::helpText("If the rate is unknown, select Unknown expected rate: no numeric rate, uncertainty or rate-source text is required. Thresholding uses dbr.sd=1; scDblFinder retains its package defaults for internal calculations, so this is not protocol-calibrated evidence."),
    shiny::conditionalPanel(condition=paste0("input['",id("rates_mode"),"'] === 'shared'"),
    number("expected_rate","Expected doublet rate (fraction 0–1; no automatic 10X assumption)",NA_real_,max=1,step=.001),
    number("uncertainty","Allowed deviation around the expected rate (± fraction)",NA_real_,max=1,step=.001),text("rate_source","Rate and uncertainty source / reasoning")),
    shiny::conditionalPanel(condition=paste0("input['",id("rates_mode"),"'] === 'csv'"),
    shiny::fileInput(id("rates_upload"),"Per-group rates: group, expected_rate, uncertainty, source (capture column also accepted)",accept=c(".csv",".tsv")),
    shiny::helpText(if(is.null(uploads$rates)) "No stored rate table." else paste("Stored rate table:",uploads$rates$source$name)),
    shiny::actionButton("qcClearRates","Clear stored rate table"))))
  if(kind=="ambient") return(shiny::tagList(ready,
    select("ambient_method","Ambient-RNA check",c("DecontX estimates (original counts retained)"="DecontX","SoupX input readiness only"="SoupX"),"DecontX"),
    shiny::helpText("DecontX estimates depend on population structure. Without a reviewed capture assignment, pooled exploration must be explicitly enabled. SoupX is readiness-only here: paired raw/background-droplet input is not yet connected."),
    check("pooled","Allow exploratory pooled DecontX with unresolved capture identity"),
    select("clusters_column","Preliminary population labels",c("DecontX heuristic clustering"="",cols)),
    check("droplet_compatible","For SoupX readiness: droplet/protocol compatibility has been established"),
    shiny::tags$details(shiny::tags$summary("DecontX model settings"),
      number("maxIter","Maximum iterations",500,min=1),number("varGenes","Variable genes for preliminary clustering",5000,min=1),
      number("delta1","Native-count prior",10,min=.001,step=.1),number("delta2","Contamination prior",10,min=.001,step=.1),
      check("estimateDelta","Estimate prior parameters",TRUE),number("convergence","Convergence setting",.001,min=.000001,step=.0001),number("dbscanEps","Preliminary DBSCAN epsilon",1,min=.001,step=.1))))
  species<-S4Vectors::metadata(object)$scrna_explorer_qc$config$organism
  shiny::tagList(ready,
    shiny::helpText("Markers define S/G2M programs. The baseline's pinned AnnotationHub table maps those markers to original feature IDs. Phase alone never excludes cells; no regression is performed."),
    select("marker_mode","Cell-cycle marker reference",c("Seurat human cc.genes.updated.2019"="seurat","Explicit species-matched marker table"="custom"),if(identical(species,"Homo sapiens")) "seurat" else "custom"),
    shiny::fileInput(id("markers_upload"),"Custom markers: program (S/G2M), gene_symbol",accept=c(".csv",".tsv")),
    shiny::helpText(if(is.null(uploads$markers)) "No stored custom marker table." else paste("Stored marker table:",uploads$markers$source$name)),
    shiny::actionButton("qcClearMarkers","Clear stored marker table"),
    text("marker_name","Custom reference name"),select("marker_species","Custom reference species",c("Select species"="","Homo sapiens","Mus musculus")),
    shiny::tags$details(shiny::tags$summary("Scoring settings"),
      number("minimum_markers","Minimum matched markers in each program",1,min=1),number("nbin","Expression bins",24,min=2),
      check("ctrl_auto","Control-gene count: use the smaller mapped program size",TRUE),number("ctrl","Control-gene count if automatic is disabled",25,min=1),
      number("scale_factor","Temporary LogNormalize scale factor",10000,min=1)))
}
