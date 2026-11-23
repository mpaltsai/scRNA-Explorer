# Existing QC with optional diagnostics

Keep the app's sliders, cutoff lines, recalculated metrics, and side-by-side pre/post-filtering plots. The standalone code tests additions at this same QC stage before changing the app. Technical tests run automatically; the user reviews one QC proposal and its scientific choices, rather than approving every assertion or helper invocation.

## Code reuse

`make_legacy_qc_steps()` sources the existing `site_data/R/createMetadata.R`, `filter_cells.R`, and `createQCPlots.R` into a private environment. It calls their original functions. The adapter supplies their existing R dependencies and redirects the relative annotation-file lookup without changing the working directory. Those source files are unchanged.

The metadata step adds the existing metrics to the existing SCE and preserves GEO/other metadata. If the old helper removes/remaps features or changes counts, the adapter stops for review. Existing gene mapping is not silently replaced. The filtering step calls `filter_cells()` with all nine explicit thresholds and can render the same ten plots before and after filtering. It requires an accepted parent when used through the runner.

Existing `cells` or `sample` fields with conflicting identities cause an error rather than being overwritten. Select the existing sample mapping explicitly when reusing an object that already has those fields.

The initial `initial_filt()` routine is unchanged in the app. This diagnostic example does not silently apply its `nUMI > 100` rule to an already supplied object. The pilot is before that rule; the report shows its prospective effect. Model-based diagnostics on low-quality cells need review; this example is not a claim that the pilot has passed QC.

## Run one combined QC proposal

In an R session at the repository root, reuse `sce` if it is already loaded:

```r
source("analysis/R/analysis_runner.R")
source("analysis/R/legacy_qc_adapter.R")
source("analysis/R/qc_optional_checks.R")
source("analysis/examples/preprocessing_qc_step.R")

# Packages installed with approval are isolated here. This affects this R session only.
project_library <- "analysis/local/R-library"
if (dir.exists(project_library)) .libPaths(c(normalizePath(project_library), .libPaths()))

input_path <- "analysis/local/GSE176031-r1/pilot_sce.rds"
# sce <- readRDS(input_path)  # only if not already available
directory <- "analysis/local/my-qc-review"
legacy <- make_legacy_qc_steps("site_data")

parameters <- list(
  app_root = "site_data",
  legacy = list(organism = "Homo sapiens", sample_column = "geo_accession",
                plot_dir = file.path(directory, "existing-qc-plots")),
  cell_cycle = list(organism = "Homo sapiens",
                    reference = seurat_cell_cycle_reference(), method_seed = 1L),
  doublets = list(capture_column = "geo_capture_id", technology = "Seq-Well S3",
                  method_seed = 1L),
  ambient = list(method = "DecontX", method_seed = 1L,
                 exploratory_pooled = TRUE),
  diagnostic_dir = file.path(directory, "optional-diagnostics")
)
proposal <- propose_analysis_step(
  sce, directory, "preprocessing-qc-review", preprocessing_qc_step,
  parameters = parameters, input_path = input_path,
  source_files = c(legacy$source_files, "analysis/R/analysis_runner.R",
                   "analysis/R/legacy_qc_adapter.R", "analysis/R/qc_optional_checks.R",
                   "analysis/examples/preprocessing_qc_step.R")
)
proposal$diagnostics$optional_checks
write_legacy_qc_review(sce, proposal, file.path(directory, "review.html"))
```

Use new evidence paths for each revision. Proposed objects stay in memory; no full object is automatically saved. The report is a local inspection view, not a replacement Shiny app. It references the immutable evidence recorded by the runner. Decisions still use `decide_analysis_step(proposal, "accept"/"reject", reason)`; no live pilot decision has been recorded on the user's behalf.

## Check definitions and limits

- **Existing gene sets:** mitochondrial genes are those annotated on `MT`, including mitochondrial coding genes, rRNAs and tRNAs. Hemoglobin genes have descriptions starting with `hemoglobin`, without a pseudogene exclusion. `audit_legacy_qc_metrics()` records every defined/matched/missing gene and checks seven per-cell metrics independently. Fractions are 0–1, not percentages. The bundled annotation file does not record its release/build; matching to hg19 remains unverified. The user chose to retain these definitions during this review.
- **Cell cycle:** `cell_cycle_qc_step()` calls the installed Seurat `NormalizeData()` and `CellCycleScoring()` on a temporary sparse representation, then returns scores and predicted phases on the original object. No scaling, regression, filtering, or final clustering is performed. Internal reversible feature aliases prevent Seurat's underscore replacement from modifying/collapsing original identifiers; the matched-marker map is recorded. Reference species, gene lists, overlap, normalization, control sampling and seed are recorded. Missing programs do not become zero scores. `Undecided` phases remain undecided. Scores alone do not establish how strongly cell cycle affects downstream PCs. [Seurat documentation](https://satijalab.org/seurat/articles/cell_cycle_vignette)
- **Doublets:** `doublet_qc_step()` operates on original counts by confirmed capture, with explicit expected rates, uncertainty and provenance. It appends `scDblFinder.score` and `scDblFinder.class`, with no exclusions. Missing capture/rate inputs give `not_assessed`, not invented singlet labels. The initial implementation does not automatically apply a 10x expected rate to Seq-Well. [scDblFinder documentation](https://bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html)
- **SoupX readiness:** `ambient_rna_readiness_step()` records unavailable/mismatched background inputs and missing capture information; it does not run estimation/correction or assign zero contamination. SoupX needs appropriate background counts and preliminary clusters for automatic estimation. [SoupX documentation](https://github.com/constantAmateur/SoupX)

The GSE176031 publication reports one to four arrays per sample, followed by separate amplification/library pools. A GEO accession is therefore not automatically a physical capture. The authors used DoubletFinder with a theoretical 5% rate; that historical setting is source evidence for review, not an automatically accepted rate in this workflow. [Study methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC8748675/)

For matrix-only ambient-RNA estimation, the user authorized installing/testing DecontX while retaining original counts. It estimates contamination from cross-cluster profiles without requiring empty droplets, but depends on preliminary cluster quality. Such clusters must stay separate from the final clustering/annotation workflow. [DecontX documentation](https://bioconductor.org/packages/release/bioc/vignettes/decontX/inst/doc/decontX.html)

The example uses estimate-only pooled DecontX. Select its mode before proposing the step:

```r
parameters$ambient <- list(method = "DecontX", method_seed = 1L,
                           exploratory_pooled = TRUE)
# With verified technical batches, use batch_column instead of pooled exploration.
# Optional clusters_column supplies reviewed preliminary labels; otherwise DecontX
# performs heuristic clustering solely for estimation.
```

Set an optional check's parameter entry to `NULL` to disable it. For SoupX input-readiness review instead, use `ambient = list(method = "SoupX", technology = "Seq-Well S3", droplet_compatible = FALSE, capture_column = "geo_capture_id")`; this reports eligibility only and does not run SoupX.

`decontx_qc_step()` keeps the original `counts` assay selected and unchanged, adds `ambient_decontX_fraction` and a separately named preliminary-cluster field, and discards its temporary corrected matrix. It records fit controls, warnings, iteration-limit flags and likelihood traces. The initial module estimates contamination only; applying or retaining corrected counts is not implemented by this step. `model_args` exposes `maxIter`, `delta`, `estimateDelta`, `convergence`, `varGenes`, and `dbscanEps`. The pilot explicitly uses pooled exploration because capture/batch mapping remains unresolved; its estimates must not be presented as validated capture-specific contamination or QC truth. Doublet scores remain unassigned for this pilot until capture and rate inputs are reviewed.

## Verification

```bash
Rscript --vanilla -e '.libPaths(c("analysis/local/R-library", .libPaths())); testthat::test_dir("tests/testthat", reporter="summary")'
```

Tests use small controlled objects and the installed backends. The real pilot separately checks unchanged counts/identities/GEO metadata and evidence generation. Passing these checks establishes software behavior, not measured biological accuracy of predictions or optimal thresholds. Scores, proposed exclusions, and any future count correction require review before proceeding downstream.
