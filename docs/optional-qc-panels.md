# Optional QC panels in Shiny

Restart the development app and open **Quality Control Plots → Optional checks**.
Doublets, ambient RNA and cell cycle are nested panels inside that same QC tab.
The existing overview, threshold controls, Preview/Apply and plots remain.

These methods run only on request, on the selected sample's fixed
**UMIs ≥100 baseline**, before integration. Preview/applied comparisons select
the matching rows of the baseline scores by stable cell ID. They do not refit a
model on each proposed subset. If only baseline is shown, select/prepare a
current preview or switch the comparison to the last applied result.

No optional score automatically excludes a cell, regresses out a signal or
selects corrected counts. The existing nine threshold rules still control
filtering. Model parameters, seed, reference identity, software versions,
warnings and status are recorded separately from a filtering acceptance.

## Doublets

1. Find the always-visible **Capture and runtime settings (per sample)** section
   above the package/status tables and method tabs. Keep capture assignment unresolved
   unless its identity is known. A file, donor or GEO accession does not establish
   a physical capture. For a complete single capture, explicitly name it and
   record its supporting evidence; alternatively select a cell-metadata column.
2. Confirm that each declared capture's complete baseline population is present.
   Use consistent capture IDs across inputs and distinct IDs for independent
   captures. A capture known to span input samples is blocked here; its complete
   population must be assembled/reviewed before scoring. Different files are
   not silently treated as independent captures.
3. Supply expected doublet rate, uncertainty and source. The shared-rate option
   applies the user-entered values to every declared capture; alternatively
   upload a CSV/TSV with `capture,expected_rate,uncertainty,source`. Fractions
   use 0–1. Literal capture IDs, including leading zeroes, are preserved.
   There is no automatically assumed 10X rate.
4. Click **Run doublet check**. scDblFinder scores each confirmed capture
   separately. Inspect score histograms, score versus UMIs/detected genes, and
   predicted counts per capture, alongside the retained subset.

Scores/classes are predictions for review, not confirmed cell identities or a
claim of calibrated classification accuracy. High counts alone do not establish
a doublet. Unresolved capture/rate evidence produces **not_assessed** and `NA`
scores; it never produces zero doublets or a singlet label by default. Capture
completeness remains a user declaration rather than independently verified truth.

For GSE176031, the pilot capture identity is still unresolved. Its doublet check
must remain unassessed until supported capture information and rate assumptions
are provided.

## Ambient RNA

Choose **DecontX estimates** and a reviewed capture assignment, or explicitly
enable **exploratory pooled DecontX** when capture identity is unresolved. A
pooled fit is labelled as pooled exploration, not assigned a fabricated capture.
Optional preliminary labels can come from metadata; otherwise DecontX performs
heuristic clustering. These are not final biological cluster annotations.

Advanced settings expose iterations, variable genes, native/contamination priors,
prior estimation, convergence setting and DBSCAN epsilon. The initial values
reuse the existing helper defaults, not empirically validated settings for every
dataset. **Run ambient-RNA check** produces contamination distributions,
relationships with UMIs/genes, preliminary-population comparisons, iteration
information and the reported likelihood history when available. Fit groups are
kept separate. Warnings and iteration-limit indicators remain visible; stopping
does not establish model validity or independent contamination truth.

Without background droplets, a fit with fewer than two preliminary populations
is unavailable evidence rather than a zero-contamination result. The method is
population-dependent even when it returns numeric estimates. Original counts
stay selected; corrected DecontX counts are neither applied nor checkpointed.

**SoupX input readiness only** reports missing prerequisites. Matching raw
background-droplet upload/pairing and SoupX estimation are not implemented by
this panel. An ordinary baseline SCE is not raw droplet background.

## Cell cycle

The default human program reference is Seurat's `cc.genes.updated.2019`.
**Marker membership and gene annotation are separate inputs**: S/G2M markers
define programs; the baseline's pinned AnnotationHub-derived table maps their
symbols to original feature identifiers. Its file checksum and table hash must
still match. There is no silent reference refresh or automatic ortholog mapping.

For another species/reference, upload a CSV/TSV with `program,gene_symbol` rows
for both `S` and `G2M`, give it a name and explicitly select its matching species.
The mouse workflow requires such an explicit species-matched reference. Retained
duplicate source features are treated as ambiguous marker matches rather than
silently choosing the first unchanged internal feature name.

Use **Review marker coverage** to inspect matched, missing and ambiguous
markers before scoring. Advanced controls expose minimum matched markers per
program, expression bins, automatic/explicit control-gene count and the temporary
LogNormalize scale factor. A minimum of one matched marker preserves the helper's
eligibility boundary; it is not a validation claim for sparse marker coverage.

**Run cell-cycle scoring** shows S versus G2M scores, predicted phase counts and
QC distributions by phase, with marker tables/download. These are program scores
and predicted phases, not independently confirmed phases. Phase alone does not
exclude a cell. Checking whether cell cycle drives downstream PCs/clustering and
deciding on regression remain later analysis decisions.

## State, dependencies and exports

- Settings, uploaded rate/marker tables and latest results are kept per sample.
  Clear stored tables explicitly or upload replacements. Tables are limited to
  5 MB and are metadata/marker inputs, not expression uploads.
- New sample visits use fresh control identities. A queued Run cannot consume
  another sample's cached controls. Changing relevant settings marks previous
  evidence **stale_settings**; plots explain the required rerun and exports use
  `NA` for stale scores. Baseline changes invalidate dependent evidence. Shared
  capture declarations in other inputs are part of doublet eligibility context.
- Busy feedback is sent before synchronous model work. A closed/superseded
  baseline cannot receive a late result. Progress feedback is not a promise of
  background execution or interruptible model fitting.
- The optional library selector reuses an existing library. In the development
  worktree it detects the already-approved common repository's
  `analysis/local/R-library`. An explicit `scrnaExplorer.qcLibrary` option or
  `SCRNA_EXPLORER_QC_LIBRARY` environment variable can supply another existing
  path. No package is installed/upgraded; library search paths and RNG state are
  restored after each run. Actual loaded method versions are recorded. SoupX
  may remain unavailable. Missing dependencies give a readable status.
- The per-cell table/download includes stable IDs, ordinary QC metrics,
  filtering selection/reasons, diagnostic status, score fields and run IDs.
  Evidence RDS downloads retain parameters, diagnostics and compact histories.
  The latest per-cell results are retained once per sample/method; histories
  retain metadata rather than another expression object per run.
- Baseline/applied SCE downloads receive `qc_diag_*` fields and reference/method
  evidence under `metadata(object)$scrna_explorer_qc$optional_diagnostics`.
  Applying filtering records the optional evidence IDs/status available at that
  decision and attaches those scores without using them as extra exclusion
  rules. Later diagnostic runs can refresh exported metadata while preserving
  the applied filtering revision. Existing foreign metadata/assays and counts
  are preserved; an unrecognised reserved-field collision stops decoration.
- Latest compatible diagnostics are attached to normalization inputs. The
  existing Seurat conversion/integration helpers still require their planned
  metadata-transfer/performance refactor; this step does not claim to fix those
  helpers or their known two-sample path.

Baseline models use the ≥100-UMI population. Earlier exploratory fits on all
766 pilot cells used a different population, so their numerical scores are not
an equivalence target for fits on the 595-cell baseline. Counts, IDs and the
unchanged ordinary QC calculations remain the preservation checks.

## Shared code and checks

`site_data/R/qc_optional_checks.R` is the shared method implementation, packaged
inside Shiny for deployment. `analysis/R/qc_optional_checks.R` remains the
standalone entry point and loads that implementation without changing the
working directory. Runner examples list both source files for provenance.
`qc_diagnostics.R` validates inputs and creates compact evidence;
`qc_diagnostic_controls.R`, `qc_diagnostic_plots.R` and
`qc_diagnostic_server.R` provide controls, plots and orchestration.

```r
testthat::test_dir("tests/testthat", filter="qc-diagnostic|qc-optional-checks",
                   reporter="summary", stop_on_failure=TRUE)
```

Tests include real cell-cycle scoring and a known synthetic-capture doublet fit,
as well as explicit model test seams for orchestration. They check unknown
capture/rate abstention, pooled consent, marker/species/reference validity,
duplicate mapping uncertainty, stable-ID joins, RNG/library restoration,
stale-score invalidation, foreign metadata protection, shared captures, export
round trips and baseline/retained plotting without subset refits.

### Local validation, 24 November 2026

The complete `tests/testthat` suite passed **823 expectations**, with no
failures, errors or skips, using the existing approved isolated library.
The full UI also rendered with all three optional panels inside the QC tab.

The real app controller was exercised with the existing GSE176031 pilot SCE:
766 imported cells produced the expected 595-cell baseline. Cell-cycle scoring
mapped 31/43 S and 47/54 G2M markers. The unresolved capture declaration left
doublet scores `NA`. Exploratory pooled DecontX was exercised with `maxIter=50`
and `varGenes=1000`; these are test settings, not a filtering recommendation.
The fit reached its iteration limit, which the diagnostics reported. All four
ambient and three cell-cycle graphics rendered. Doublet graphics and an actual
scDblFinder fit were checked separately with a known synthetic capture.

A temporary threshold preview checked baseline/retained score comparisons.
No filtering decision or expression checkpoint was saved. Baseline counts,
source-file checksums and the five ordinary QC metrics matched their preserved
inputs/historical CSV by cell identifier (floating-point tolerance `1e-12`).
These checks establish software behavior on the tested inputs, not independent
biological validation of the model estimates.

### Startup validation fix, 26 November 2026

The plot renderer now checks readiness before entering plotting functions.
With `SingleCellExperiment` attached, Bioconductor's S4 dispatch previously
wrapped the normal "Review settings and run this optional check first" Shiny
validation into an ordinary error when forcing lazy plotting arguments.
The regression test opens all ten optional plots with that package attached,
before any check is requested. It requires ordinary Shiny validation and checks
that no model runs, no filtering is applied and the baseline remains unchanged.
Existing tests still exercise scored plots, stale settings and subset reuse.

### Deferred-check cleanup fix, 26 November 2026

Missing capture/rate information leaves the doublet check `not_assessed`, with
unavailable scores and an explanatory message. The deferred callback now resets
its busy state inside a separate `isolate()`: its exit handler runs after the
main isolated body has returned. Reading reactive state there without its own
context previously caused a fatal error, even after a message was prepared.

`test-qc-diagnostic-callbacks.R` replays the controller's actual deferred callback
with a session domain but no reactive consumer, matching the real app. The
standard Shiny mock normally supplies isolation and masked this failure.
The test checks missing captures, missing rates, a simulated runner failure,
successful retry after correction and marker-only review's early return. It
requires cleared busy state and an open session, preserves baseline counts and
never applies filtering. Its model-fitting seam is synthetic; it does not
assign real capture identities or validate biological doublet predictions.

Validation: all 174 optional-QC expectations passed, with no failures, errors
or skips, including 36 callback-lifecycle checks. The callback test reproduced
the reported reactive-context/`setdiff` crash before the cleanup fix.
