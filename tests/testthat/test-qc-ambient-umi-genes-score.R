testthat::test_that("ambient UMI-gene colours follow cell-ID fractions with common baseline/retained scales", {
  f<-qc_diag_fixture();on.exit(unlink(f$root,recursive=TRUE))
  b<-f$baseline;before<-digest::digest(b);m<-qc_review_metrics(b)
  # Include both endpoints and a missing estimate: missing must not appear as 0%.
  fractions<-c(0,.01,.1,.2,.4,.6,1,NA_real_)
  values<-data.frame(cell_id=m$cell_id,qc_diag_ambient_fraction=fractions,
    qc_diag_ambient_fit_group="pooled",qc_diag_ambient_population=rep(c("1","2"),4))
  r<-list(status="estimated",values=values[8:1,,drop=FALSE])
  t<-qc_default_thresholds(m);t$maxUMIsPerCell<-floor(stats::median(m$nUMI))
  selection<-qc_filter_proposal(m,t,qc_baseline_signature(b))
  for(label in c("Preview retained","Applied retained")) {
    d<-qc_diagnostic_comparison(qc_comparison_metrics(m,selection,after_label=label),r)
    plot<-qc_plot_ambient_umi_genes_score(d);built<-ggplot2::ggplot_build(plot)
    expected<-100*values$qc_diag_ambient_fraction[match(plot$data$cell_id,values$cell_id)]
    testthat::expect_identical(plot$data$ambient_plot_percent,expected)
    testthat::expect_equal(nrow(plot$data),ncol(b)+sum(selection$keep))
    testthat::expect_equal(as.numeric(built$data[[1]]$x),log10(plot$data$nUMI))
    testthat::expect_equal(as.numeric(built$data[[1]]$y),log10(plot$data$nGene))
    colour<-built$plot$scales$get_scales("colour")
    testthat::expect_equal(colour$get_limits(),c(0,100))
    testthat::expect_identical(built$data[[1]]$colour,colour$map(expected))
    testthat::expect_true(all(built$data[[1]]$colour[is.na(expected)]=="#777777"))
    testthat::expect_true(all(diff(stats::na.omit(plot$data$ambient_plot_percent))>=0))
    testthat::expect_length(built$layout$panel_scales_x,1L)
    testthat::expect_length(built$layout$panel_scales_y,1L)
  }
  testthat::expect_identical(digest::digest(b),before)
})
