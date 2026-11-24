# Backward-compatible standalone entry point. The shared implementation lives
# inside site_data so a deployed Shiny app does not depend on analysis/ files.
(function(target) {
  paths<-unlist(lapply(sys.frames(),function(frame) {
    p<-frame$ofile;if(is.character(p)) p else character()
  }),use.names=FALSE)
  paths<-paths[basename(paths)=="qc_optional_checks.R"]
  relative<-if(length(paths)) file.path(dirname(utils::tail(paths,1)),"..","..","site_data","R","qc_optional_checks.R") else character()
  candidates<-c(relative,file.path("site_data","R","qc_optional_checks.R"),file.path("..","..","site_data","R","qc_optional_checks.R"))
  candidates<-candidates[file.exists(candidates)]
  if(!length(candidates)) stop("Cannot locate the shared site_data/R/qc_optional_checks.R implementation.")
  sys.source(normalizePath(candidates[1],mustWork=TRUE),envir=target)
})(environment())
