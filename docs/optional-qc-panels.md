# Optional QC panels in Shiny

Restart the development app and open **Quality Control Plots → Optional checks**.
Doublets, ambient RNA and cell cycle are nested panels inside that same QC tab.
The existing overview, threshold controls, Preview/Apply and plots remain.

These methods run only on request, on the selected sample's fixed
**UMIs ≥100 baseline**, before integration. Preview/applied comparisons select
the matching rows of the baseline scores by stable cell ID. They do not refit a
model on each proposed subset. If only baseline is shown, select/prepare a
current preview or switch the comparison to the last applied result.

Running a check does not itself exclude cells, regress out a signal or select
corrected counts. Filtering uses the existing nine rules plus any explicitly
enabled optional exclusions described below. Model settings and evidence remain
separate from the user's acceptance of a filtering decision.

## Optional exclusions and QC layout

**Overview** stays first, with **QC plots** immediately beside it. The **Cells
and decisions** tab is removed; per-cell metrics/reasons and decision histories
remain downloadable from the sidebar, and Optional checks still has its
per-cell evidence table.

In the sidebar, **Optional diagnostic exclusions** defaults to off:

- Choose **Exclude cells predicted as doublets**, or **Exclude scores above my
  maximum** and enter a maximum retained doublet score.
- Enable ambient exclusion and enter a maximum retained DecontX fraction (0–1).
  There is no preselected scientific cutoff. A cell exactly at a maximum stays.
- Click **Preview filtering**, review combined retained/excluded counts and rule
  overlaps, then **Apply preview to this sample** with a reason if accepted.

These rules use the current baseline diagnostic results, including explicitly
exploratory results if the user selects them. Missing, stale, incomplete,
mismatched or nonfinite evidence blocks only an enabled optional rule; it is
never converted into zero contamination or a singlet classification. Disabling
optional rules restores the original nine-rule workflow without requiring a
model fit. Cell-cycle phase and preliminary cluster membership are not exclusion
rules. Overlapping failures exclude a cell once and preserve every reason.

Changing a selected cutoff or rerunning an enabled check invalidates the preview.
The last applied SCE remains intact until a new preview is applied. Applying
rechecks rules, cell-ID alignment, baseline identity, run IDs and used-value
hashes against current evidence. Each decision records selected rules, model
settings/provenance and the actual per-cell scores or labels used (`*_used`
fields), so later reruns or refreshed diagnostic exports do not rewrite the
accepted evidence. No expression checkpoint is added. Sensitivity varies the
nine standard cutoffs while holding the snapshot's optional-exclusion mask fixed.

## Doublets

The **Doublet analysis mode** selector provides two explicit paths. Switching
mode, grouping or rate handling invalidates earlier scores until the next run.

**Capture-informed** (the existing default) uses the always-visible capture
declarations above the package/status tables. Name a confirmed single capture
or select a capture metadata column, record supporting evidence and confirm
the complete baseline population. A capture known to span multiple inputs is
blocked on this path until assembled/reviewed together. Files, donors and GEO
accessions are not silently declared physical captures.

**Public-data exploratory** allows scoring when capture information is missing:

1. Leave the shared capture assignment unresolved. Select **Public-data
   exploratory: capture identity remains unknown** in the Doublets panel.
2. Select **Selected input matrix as one provisional group**, or choose an
   existing cell-metadata column for provisional groups. Every cell must have a
   nonempty label in the chosen column. GEO metadata already attached to an SCE
   can supply sample/library groups; this panel does not add a series-matrix
   importer or establish that a GEO sample is one physical capture.
3. Under **Expected doublet rates**, choose **Unknown expected rate** if there is
   no supported estimate. No rate, uncertainty, rate-source text, capture ID,
   capture-evidence text or completeness confirmation is required on this path.
4. Click **Run doublet check**. Inspect the distributions, relationships with
   UMIs/genes and predicted counts, with facets labelled **Provisional group**.

The Doublets panel also includes **UMIs versus detected genes: doublet score**
alongside the original three plots. Like the existing mitochondrial-coloured
QC scatterplot, it uses log-scaled UMIs on x and detected genes on y. Point
colour is the continuous baseline scDblFinder score, with a shared 0–1 colour
scale and fixed count axes across baseline/retained facets and analysis groups.
Higher scores are drawn last for visibility. This plot's two count axes stay
fixed when the separate diagnostic-relationship selector changes. Plotting or
changing a QC subset reuses the existing scores and does not rerun detection.
The plot extension passed 145 targeted expectations, including score-to-cell
alignment, log-count coordinates, shared facet/colour scales, pre-run validation
and rendering in capture-informed and exploratory modes. A full UI render
confirmed that all four plots remain inside the Doublets tab.

Groups are scored separately within the selected sample's fixed baseline.
Exploratory groups may mix or split physical captures. Use the broad available
population before cell-type subsetting; small or unrepresentative populations
can make predictions unreliable or cause the method to fail with a message.
No cell is dropped to make a grouping complete. Missing metadata labels are
unavailable grouping information, not a new biological group.

Both modes also support explicit shared rates or uploaded per-group rates.
The CSV/TSV columns are `group,expected_rate,uncertainty,source`; the earlier
`capture` column is accepted as an alias for `group`. Literal identifiers,
including leading zeroes, are retained. Conflicting group/capture columns are
rejected. Explicit rates need reviewed numeric fractions and source/reasoning;
uncertainty is an absolute +/- margin, not a confidence interval.

Unknown-rate handling passes `dbr=NULL, dbr.sd=1`: thresholding does not enforce
an expected doublet count. This is **not a claim that all rate assumptions have
been removed**. scDblFinder retains its internal package defaults, including
`dbr.per1k`, which are recorded with software versions and the threshold policy;
they are not protocol calibration for Seq-Well or another unverified platform.
See [the authors' parameter guidance](https://bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html).

Scores/classes are predictions for review, not confirmed cell identities or a
claim of calibrated classification accuracy. High counts alone do not establish
a doublet. Incomplete capture declarations still produce **not_assessed** on
the capture-informed path. Exploration is an explicit opt-in, not an automatic
conversion of missing information into singlet labels. Capture completeness
remains a user declaration rather than independently verified truth.

For GSE176031, the pilot capture identity remains unresolved. It can now receive
exploratory predictions by its input matrix (or reviewed metadata groups), with
unknown-rate handling. Exports retain `qc_diag_capture=NA`, capture status
`unknown`, the actual provisional group, grouping basis, analysis mode and rate
mode. Confirmed-capture predictions require their separate documented inputs.
Neither path automatically accepts a filtering decision or changes counts.

Validation of the exploratory extension (26 November 2026): all 955 recorded
expectations in the full suite passed, with no failures, errors or skips.
The real app controller scored the existing 595-cell GSE176031 baseline using
the input-matrix grouping and unknown-rate handling, with no capture declaration
or invented numeric rate. It produced 96 provisional doublet predictions;
these are not independent quality truth or an exclusion recommendation. All
three doublet graphics rendered. A temporary 535-cell QC preview reused the
baseline scores. Capture IDs stayed unavailable in exports; original counts,
baseline and source-file checksum stayed unchanged. No filtering decision or
expression checkpoint was saved. The synthetic-backend tests cover grouping
changes, missing metadata labels, old rate-table aliases and capture-informed
abstention; a separate real synthetic-data fit exercises the shared method.
A follow-up compatibility run passed 122 expectations, including preservation
of the earlier capture-keyed rate-table output for capture-informed runs.

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

**UMIs versus detected genes: preliminary ambient-RNA clusters** uses log-scaled
count axes and categorical colours for DecontX's preliminary population labels.
Baseline and retained panels share count scales and cluster colours; fits are
faceted separately because a cluster label is local to its fit. These are not
final clustering/cell-type annotations. Cluster membership itself is not an
exclusion rule: the optional ambient rule uses each cell's estimated fraction.

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
  decision and attaches current scores. Selected optional rules additionally
  retain their decision-time values and references in the filtering history.
  Later diagnostic runs can refresh current exported diagnostics while preserving
  the applied filtering revision and its immutable `*_used` values. Foreign metadata/assays and counts
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
The regression test opens all optional plots with that package attached,
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
