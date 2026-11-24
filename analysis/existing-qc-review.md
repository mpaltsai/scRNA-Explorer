# Existing QC with optional diagnostics

Keep the app's sliders, cutoff lines, recalculated metrics, and side-by-side pre/post-filtering plots. The standalone code tests additions at this same QC stage before changing the app. Technical tests run automatically; the user reviews one QC proposal and its scientific choices, rather than approving every assertion or helper invocation.

The development app now has a [QC review checkpoint](../docs/qc-review-workflow.md)
with new plots and explicit per-sample Preview/Apply decisions. That guide
describes the current QC tab, which replaces separate metrics/redraw button
sequences. The standalone diagnostic instructions below remain available;
optional doublet, ambient-RNA and cell-cycle panels are now connected as
diagnostic evidence in the same QC tab; see [optional QC panels](../docs/optional-qc-panels.md).

## Shiny input integration: first review checkpoint

Development of the app now takes place in the separate `feat/shiny-input-integration` worktree. The original checkout's tracked `site_data` remains at its previous state; it is not the place to run the new input workflow. From the development checkout, run:

```r
shiny::runApp("site_data")
```

For this checkpoint, test **Read data** and the import review panel. The existing input choices are retained. CSV/TSV readers preserve the first feature-ID column and cell column names; 10X readers accept uncompressed/gzipped files, use primary feature IDs and retain symbols separately. Supplied barcodes remain intact; missing barcodes receive reproducible generated IDs marked as generated. ZIP extraction is session-specific. The single-sample RDS path requires a SingleCellExperiment counts assay and retains the original object and metadata in the import records. This is not arbitrary processed-Seurat or multi-sample RDS restoration.

Negative, missing, nonfinite or fractional UMI counts stop import with a message. Duplicate primary feature IDs pause import and display their rows and totals. The initial resolution controls offer **Keep all rows separately** with unique internal IDs, or **Cancel import; correct the input**. Both require a reason, which can be downloaded. Keeping separate rows preserves counts and records original IDs; this first version does not merge counts or drop rows. Every sample must be ready before records reach the annotation input. A failed new read clears prior imported and dependent analysis state.

The imported SCE, metadata and identity maps are retained in `my_react_vals$import_records`. Baseline propagation is covered by the second checkpoint below; explicit filtering decisions are covered by the QC review guide. Downstream Seurat conversion remains a later step. An import test does not validate those stages or biological accuracy of optional QC methods.

Changing the input mode or an active upload clears stale completion/results and shows **Not yet read**. **Read data** sends a **Reading data…** label and busy indicator before starting synchronous file processing. Completion and the import summary are published together afterward. The controller captures the input selection on click and ignores superseded work; inactive upload controls and annotation settings do not reset a valid import. This improves feedback, not file-reading speed or background concurrency.

Manual checks:

1. Import the existing pilot SCE RDS and verify the import summary reports 19,665 features and 766 cells, with no generated cell IDs.
2. Import a small count table with a repeated primary feature ID; verify that the annotation input remains blocked until a decision with a reason is recorded. Inspect/download the decision record.
3. Start another read with invalid UMI counts; verify that an error is displayed and the previous input is not reused.
4. For a ZIP containing multiple samples, check that each sample and its counts are reported separately.

The focused automated checks run without starting a web server:

```bash
Rscript --vanilla -e 'testthat::test_dir("tests/testthat", filter="input-readers|input-import-server|import-feedback", reporter="summary", stop_on_failure=TRUE)'
```

They exercise the readers, compressed inputs, metadata preservation, duplicate decisions, session isolation, and actual Shiny server input wiring through `shiny::testServer()`. No full expression checkpoint is saved per import decision.

## Gene annotation and baseline SCE: second review checkpoint

After **Read data**, explicitly select the organism and each sample's feature-identifier type. The 10X reader declares gene IDs from its feature file; table and RDS inputs require a choice rather than an ENSG-prefix guess. Enable Ensembl-version stripping only when wanted; this changes lookup keys, never the count row names. Click **Add gene annotations** and inspect the reference provenance, mapping/metric-availability summary and first 20 mapping issues. The complete mapping and pre-initial-filter metrics can be downloaded.

`site_data/R/gene_annotations.R` provides the shared cache loader, mapping and metric preparation, and baseline constructor. Cell-cycle reference loading now uses the same cache loader. The app reads the bundled AnnotationHub-derived caches without a network lookup or new package installation. Resource ID, snapshot, release and genome build remain unknown until established; a user-confirmed hg19 count alignment does not verify the annotation cache's build.

Features are never dropped or renamed as an annotation shortcut. Ambiguous or unmapped annotations remain visible, including candidate gene IDs. The existing MT-chromosome and hemoglobin-description gene definitions determine metric membership in the explicitly selected identifier type. In symbol mode these are the same symbol sets as the legacy app; set membership does not establish a unique gene-ID assignment for an ambiguous symbol. Both source rows retained by an explicit duplicate-ID decision remain distinct and contribute their counts. An absent gene set or no matching features produces unavailable (`NA`) fractions, whereas a supplied matching gene with zero counts yields observed zero expression. The existing nine-metric filter is blocked when required evidence is unavailable; adapting its plots/controls to reduced-evidence QC is a later checkpoint.

Click **Minimal filtering** to create each sample's baseline SCE using `nUMI >= 100`. The original SCE is reused where supplied, preserving other assays, reductions and source row/cell/experiment metadata. Matrix inputs become an SCE at this stage. Preparation retains only annotation and metric tables; no extra full expression checkpoint is written. Empty baselines are reported and omitted from the QC sample selector. Replacing the import or changing annotation settings invalidates baseline and active downstream inputs.

Current metrics and annotations use `qc_*` fields, including `qc_nUMI`, `qc_sample`, `qc_gene_id`, `qc_source_feature_id` and sample-qualified `qc_cell_id`. Imported fields such as `sample`, `nUMI` or `gene_name` remain intact, even when they conflict with the new current values; conflicts are reported. Legacy aliases are added only where absent. QC views use `baseline_qc_metrics()` to read current values. Unexpected collisions with the reserved `qc_*` namespace stop for review; recognised app-generated metadata can be updated while preserving the earlier QC recipe. Reference/code checksums, gene sets, mappings summary, import decisions and initial exclusions are recorded in `metadata(sce)$scrna_explorer_qc`.

Ribosomal feature annotations are now added for later review, without gene exclusion, new count denominators or a variable-feature/PCA mask:

- `qc_ribosomal_protein`: protein-coding Entrez genes supported by GOALL structural-constituent-of-ribosome annotations (`GO:0003735`) in the installed species-specific OrgDb.
- `qc_cytosolic_ribosomal_protein` and `qc_mitochondrial_ribosomal_protein`: subsets also annotated to `GO:0022626` and `GO:0005761`, respectively. They overlap the parent protein group, so their counts must not be added together.
- `qc_nuclear_rRNA`, `qc_mitochondrial_rRNA` and `qc_rRNA_pseudogene`: separate cached biotype/chromosome-derived RNA categories. Functional rRNA and pseudogenes are not combined. A mitochondrial ribosomal-protein annotation does not imply a gene lies on chromosome MT; the existing mitochondrial fraction definition remains unchanged.

These are reference-supported classifications, not a complete inventory of all ribosomal functions. False means not annotated under this definition, not a confirmed biological negative. Missing Entrez links/resources and conflicting candidate classifications remain `NA`/unresolved; no first candidate is selected. Broad symbol prefixes are deliberately insufficient: RPS6 kinases must not become structural ribosomal proteins, while fusion genes such as FAU/UBA52 can be recognised through their reference IDs. Installed `org.Hs.eg.db`/`org.Mm.eg.db` resources provide species-specific memberships without ortholog conversion, package installation or data download. Resource versions, source dates, GO terms, explicit gene sets and their hash are saved in the QC recipe; the small public reference sets are cached in memory only. If resources are unavailable, protein classifications remain unavailable and RNA biotypes can still be reviewed.

The review panel summarises each group and unresolved counts; the full mapping download and baseline SCE contain the feature flags. Selecting any group for exclusion remains a separate future decision, including whether exclusion should affect variable-feature selection or an explicitly derived analysis matrix. Original counts remain the reference.

Definition references: [GO structural constituent of ribosome](https://amigo.geneontology.org/amigo/term/GO:0003735), [cytosolic ribosome](https://amigo.geneontology.org/amigo/term/GO:0022626), [mitochondrial ribosome](https://amigo.geneontology.org/amigo/term/GO:0005761), and [Bioconductor annotation interfaces](https://bioconductor.org/packages/release/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf).

For the pilot: choose **Homo sapiens** and **Gene symbols**, add annotations, then run minimal filtering. Expected baseline: **19,665 features and 595 cells**. All five earlier basic QC metrics agree for all 766 input cells. The cached reference yields 13,115 uniquely matched, 2,048 ambiguous and 4,502 unmapped features, with 28 MT and five hemoglobin feature matches. These are mapping diagnostics, not biological quality labels. Review the downloaded baseline SCE before proceeding to filtering orchestration or optional diagnostic controls.

Focused tests:

```bash
Rscript --vanilla -e 'testthat::test_dir("tests/testthat", filter="gene-annotations|annotation-qc-server|cell-cycle-annotation", reporter="summary", stop_on_failure=TRUE)'
```

## Code reuse

`make_legacy_qc_steps()` sources the existing `site_data/R/createMetadata.R`, `filter_cells.R`, and `createQCPlots.R` into a private environment. It also records/loads the shared `gene_annotations.R` source. The standalone legacy adapter calls the original metric/filter/plot functions, supplies their existing R dependencies and redirects the relative annotation-file lookup without changing the working directory. The three legacy helper files are unchanged. The development Shiny app uses the new identifier-preserving preparation/baseline path described above.

The metadata step adds the existing metrics to the existing SCE and preserves GEO/other metadata. If the old helper removes/remaps features or changes counts, the adapter stops for review. Existing gene mapping is not silently replaced. The filtering step calls `filter_cells()` with all nine explicit thresholds and can render the same ten plots before and after filtering. It requires an accepted parent when used through the runner.

Existing `cells` or `sample` fields with conflicting identities cause an error rather than being overwritten. Select the existing sample mapping explicitly when reusing an object that already has those fields.

The initial `initial_filt()` routine keeps cells with `nUMI >= 100`, as explicitly confirmed by the user: cells with exactly 100 UMIs remain. This diagnostic example does not silently apply initialization to an already supplied object. The pilot is before that rule; newly generated reports show its prospective effect. Earlier saved pilot reports used the previous strict `nUMI > 100` boundary and remain historical evidence. Model-based diagnostics on low-quality cells need review; this example is not a claim that the pilot has passed QC.

The baseline UMI filter remains in the workflow. Doublet detection, ambient-RNA assessment, cell-cycle scoring and experimental threshold suggestions are optional additions at the QC stage. Running a diagnostic does not automatically exclude cells; applying additional exclusions requires an explicit user decision. Existing sliders and pre/post-filtering plots remain the manual review workflow.

## Run one combined QC proposal

In an R session at the repository root, reuse `sce` if it is already loaded:

```r
source("analysis/R/analysis_runner.R")
source("analysis/R/legacy_qc_adapter.R")
source("analysis/R/qc_optional_checks.R")
source("analysis/examples/preprocessing_qc_step.R")

# Reuse the original checkout's data/packages from the isolated development
# checkout. Set data_root to "." when running in the original checkout instead.
data_root <- "../scRNA-Explorer"
project_library <- file.path(data_root, "analysis/local/R-library")
if (dir.exists(project_library)) .libPaths(c(normalizePath(project_library), .libPaths()))

input_path <- file.path(data_root, "analysis/local/GSE176031-r1/pilot_sce.rds")
# sce <- readRDS(input_path)  # only if not already available
directory <- "analysis/local/my-qc-review"
legacy <- make_legacy_qc_steps("site_data")

parameters <- list(
  app_root = "site_data",
  legacy = list(organism = "Homo sapiens", sample_column = "geo_accession",
                plot_dir = file.path(directory, "existing-qc-plots")),
  cell_cycle = list(organism = "Homo sapiens",
                    annotation = load_cell_cycle_annotation("Homo sapiens", "gene_name", "site_data"),
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
                   "analysis/R/legacy_qc_adapter.R", "analysis/R/qc_optional_checks.R", "site_data/R/qc_optional_checks.R",
                   "analysis/examples/preprocessing_qc_step.R")
)
proposal$diagnostics$optional_checks
write_legacy_qc_review(sce, proposal, file.path(directory, "review.html"))
```

Use new evidence paths for each revision. Proposed objects stay in memory; no full object is automatically saved. The report is a local inspection view, not a replacement Shiny app. It references the immutable evidence recorded by the runner. Decisions still use `decide_analysis_step(proposal, "accept"/"reject", reason)`; no live pilot decision has been recorded on the user's behalf.

## AnnotationHub mapping for cell-cycle scoring

Cell-cycle marker membership and gene annotation are separate inputs. The Seurat reference supplies S/G2M marker symbols; the AnnotationHub-derived table maps them to the object's original identifiers. Scoring now requires an explicit annotation input and identifier type rather than assuming row names are symbols. The example above reuses the app's existing cached annotation table without contacting AnnotationHub. Its file checksum is recorded; its missing resource ID, snapshot date, release and genome build remain unknown. A confirmed hg19 count alignment alone does not establish the cache's build. The loader also accepts an explicit `reference_id` and `cache_dir` for separately downloaded references; see [local cache maintenance](../docs/annotation-reference-cache.md).

For an already loaded AnnotationHub-derived table with `gene_id` and `gene_name`, use:

```r
parameters$cell_cycle$annotation <- list(
  table = annotations,
  species = "Homo sapiens",
  feature_id_type = "gene_id",  # use gene_name for symbol row names
  provenance = list(source = "User-supplied AnnotationHub-derived table")
)
# Include the actual resource ID, snapshot date, release and build in provenance
# when available; do not infer them from this example.
parameters$cell_cycle$strip_ensembl_version <- TRUE
```

Version stripping is opt-in, applies only to Ensembl gene-ID lookup keys and never renames matrix rows. Duplicate annotation pairs are collapsed. Multiple IDs for a marker, conflicting symbols for an ID, multiple matching matrix rows, and features shared by both programs are excluded from scoring as ambiguous. Unannotated symbols do not fall back to direct row-name matching. Missing annotations and missing matrix features have separate statuses. A program with no usable markers stops rather than assigning zero scores. Partial coverage remains diagnostic and is not evidence of accurate phase assignment. Mouse annotations require an explicit mouse-compatible marker reference; automatic ortholog conversion is still disabled.

The combined example exports `cell-cycle-marker-coverage.csv` with marker symbols, candidate IDs/features, mapping statuses and matched original identifiers, plus `cell-cycle-annotation-provenance.txt`. The diagnostics retain the annotation-table SHA-256 fingerprint, matched/missing lists and scoring settings. Existing saved pilot outputs are historical; this update does not overwrite them or accept a new runner decision.

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
