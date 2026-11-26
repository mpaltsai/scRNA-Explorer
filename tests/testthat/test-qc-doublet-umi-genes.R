testthat::test_that("UMI-gene colour plots reuse cell-ID scores with shared count and colour scales", {
  f<-qc_diag_fixture();on.exit(unlink(f$root,recursive=TRUE))
  b<-f$baseline;before<-digest::digest(b);hash<-qc_baseline_signature(b)
  cfg<-qc_diag_config("doublets",f$app_root)
  cfg$analysis_mode<-"exploratory";cfg$group_mode<-"input";cfg$rates_mode<-"unknown"
  r<-qc_run_optional(b,"doublets",cfg,hash,f$app_root,execute=qc_diag_fake_fit)
  # Shuffling the exported table must not change which cell gets each colour.
  r$values<-r$values[nrow(r$values):1,,drop=FALSE]
  m<-qc_review_metrics(b);t<-qc_default_thresholds(m)
  t$maxUMIsPerCell<-floor(stats::median(m$nUMI))
  selection<-qc_filter_proposal(m,t,hash)
  d<-qc_diagnostic_comparison(qc_comparison_metrics(m,selection),r)
  plot<-qc_plot_doublet_umi_genes(d);built<-ggplot2::ggplot_build(plot)
  expected<-r$values$qc_diag_doublet_score[match(plot$data$cell_id,r$values$cell_id)]
  testthat::expect_identical(plot$data$qc_diag_doublet_score,expected)
  testthat::expect_equal(nrow(plot$data),ncol(b)+sum(selection$keep))
  testthat::expect_equal(as.numeric(built$data[[1]]$x),log10(plot$data$nUMI))
  testthat::expect_equal(as.numeric(built$data[[1]]$y),log10(plot$data$nGene))
  colour<-built$plot$scales$get_scales("colour")
  testthat::expect_equal(colour$get_limits(),c(0,1))
  testthat::expect_identical(built$data[[1]]$colour,colour$map(expected))
  testthat::expect_length(built$layout$panel_scales_x,1L)
  testthat::expect_length(built$layout$panel_scales_y,1L)
  testthat::expect_true(all(grepl("Provisional group:",plot$data$qc_plot_group,fixed=TRUE)))
  testthat::expect_identical(digest::digest(b),before)
})
