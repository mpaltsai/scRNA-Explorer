# QC metric review and filtering

This development step preserves the existing data lifecycle: one fixed baseline
SCE per sample with `nUMI >= 100`, and one current applied QC-filtered SCE per
sample. Multiple samples are filtered separately. The filtered SCEs remain the
inputs to the existing Seurat workflow. Neither previews nor decision histories
contain additional expression objects.

## Try the development app

Restart the development app after editing; use the isolated worktree, not the
original checkout:

```r
shiny::runApp("/home/ismini/Documents/RNA_cor/iRNA/scRNA_Explorer_gitHub/scRNA-Explorer-shiny-input-integration/site_data")
```

1. Read data, choose annotation settings, **Add gene annotations**, inspect the
   mapping/reference and perform **Minimal filtering**. The original baseline
   boundary remains inclusive at 100 UMIs.
2. Open **Quality Control Plots** and select a sample. Metrics are available
   automatically for the selected baseline; separate create/recompute/redraw
   button sequences are no longer required.
3. Review **Overview**, **Joint metrics**, **Ribosomal fractions** and
   **Annotation coverage**. The overview includes other baseline samples for
   comparison, but filtering settings remain per sample.
4. Edit numeric draft cutoffs. The initial settings retain all baseline cells
   when the required metrics are available; they are not recommendations.
   Optional sliders remain available under **Slider controls**. As in the
   original app, explicitly copy their values into the numeric draft. Exact
   numeric cutoffs are not silently rounded to slider precision.
5. Click **Preview filtering**. Inspect retained/excluded counts, **Filtering
   impact**, rule overlaps, **Sensitivity** and the per-cell reasons table.
   Hover over joint-metric points to inspect stable cell IDs. A preview changes
   tables/masks/plots only; the applied SCE and its download remain unchanged.
6. Give a reason and click **Apply preview to this sample**. Apply verifies the
   preview against the current baseline and numeric cutoffs. Editing a cutoff,
   switching to an unreviewed sample or changing the baseline cannot apply an
   old preview. All nine existing rules have inclusive bounds.
7. Use **Compare baseline with** to distinguish the current preview from the
   last applied result. An edited/unpreviewed draft shows only the baseline in
   preview mode. Applied mode continues to show the recorded applied result.
   The comparison status also identifies the stored applied revision.
8. Switch samples and repeat. Other applied samples remain intact. Previously
   entered drafts are retained per sample. Newly rendered controls use distinct
   sample/generation IDs so values from another sample/import cannot leak in.
9. Download the applied SCE, per-cell review CSV, sensitivity CSV or compact
   decision histories as needed. Decision histories also travel in the applied
   SCE's `metadata(object)$scrna_explorer_qc$filtering`. They include thresholds,
   reasons, membership, previous/current revisions, code and reference identity.

Apply replaces that sample's current filtered SCE and invalidates combined
normalization/clustering/correlation results. Normalization requires an applied,
nonempty result for **every** sample; a preview is insufficient. The resulting
Seurat object records the applied sample revision IDs in
`object@misc$scrna_explorer_qc$filtered_sample_revisions`.

An empty retained selection can be inspected and applied explicitly, but blocks
downstream analysis until revised. Apply and Preview failures preserve the
previous applied data and report the reason. Baseline/reference invalidation
clears affected QC state through the existing input/annotation lifecycle.

## Interpretation of the plots

- Overview: boxplots for UMIs, detected genes, mitochondrial/hemoglobin fractions
  and the original complexity metric, `log10(nGene)/log10(nUMI)`. Counts are
  shown on log10 axes. This complexity statistic is not `nGene/nUMI` or sequencing
  saturation.
- Joint metrics: mitochondrial UMI percentage against log10 genes or UMIs,
  with proposed/applied exclusions highlighted on the baseline and retained
  cells alongside. Stages share axes and colours. Dashed cutoffs belong to the
  selected preview/applied snapshot, rather than unreviewed edits.
- Ribosomal fractions: the six existing annotation groups remain separate.
  Every denominator is all original UMIs. Positive annotations contribute to
  the numerator; unresolved-feature UMI shares are reported alongside. Partial
  annotated shares are not complete biological classification. An unavailable
  protein reference yields unavailable protein fractions, not zero. Parent and
  subgroup protein fractions overlap and must not be summed. These plots do not
  remove genes, introduce ribosomal filtering or change variable-feature choices.
- Annotation coverage: matched/ambiguous/partial/unmapped gene annotations are
  weighted by the original UMI counts they contribute, not just feature numbers.
- Filtering impact: per-rule failures can overlap. Exact rule combinations count
  each excluded cell once; the graph shows the top 15 combinations and the cell
  export retains all reasons. Zero exclusions remains a valid preview.
- Sensitivity: one cutoff varies while the other eight stay fixed at the
  selected preview/applied snapshot. The export records the mode, baseline hash
  and selection/revision ID. Fraction
  cutoff axes use 0–1, while retained-cell proportions use percent. This measures
  sensitivity to settings, not classification accuracy or optimality.
- Existing QC plots: the original plotting functions remain available in their
  own panel. Baseline/retained facets use common scales. Library composition and
  highly expressed gene plots run on request; retained expression plots use the
  applied SCE. A sample/revision change clears an incompatible expression view.

Missing/nonfinite evidence in any of the nine filtering metrics blocks preview
instead of silently removing those cells. Count-only plots remain inspectable;
the sample summary reports unavailable mitochondrial/hemoglobin metrics. A
reduced-evidence filter requires a separate explicit design and validation.

## Code reused and new code

- Reuse `baseline_qc_metrics()` and existing SCE assays/annotations. Existing
  imported cell/feature/experiment metadata and additional assays survive Apply.
- Keep `createQCPlots.R` and `filter_cells.R` unchanged. The original plot
  functions are called directly. The new mask evaluator expresses the same nine
  inclusive rules and is checked against `filter_cells()`; masks also provide
  per-rule explanations and preview tables without creating an SCE.
- `site_data/R/qc_review.R`: metric tables, validation, masks, application,
  sensitivity and all-sample eligibility.
- `site_data/R/qc_review_plots.R`: table-only scientific plotting helpers.
- `site_data/R/qc_review_server.R`: per-sample state, sample-specific controls,
  plots, hover information, downloads and decision histories.
- The old QC event block in `server.R` is replaced by the controller and small
  compatibility readers. Downstream normalization methods are unchanged; this
  step adds eligibility checks and input revision tracking, not a redesign of
  normalization/integration. Known issues in those helpers, including the
  two-sample path and dense conversion, still need the later Seurat stage.

Doublet, ambient-RNA and cell-cycle controls/plots are now available under
**Optional checks** inside this QC tab; see [their workflow](optional-qc-panels.md).
They are not used to exclude cells by this controller. Model-based threshold
suggestions remain a later integration step.
Top-gene concentration extensions and a full analysis-bundle restore workflow
are also separate work. No reference download, dependency install, source-data
filtering, commit or push is performed by these tests.

## Automated and manual checks

```r
testthat::test_dir("tests/testthat", filter="qc-review",
                   reporter="summary", stop_on_failure=TRUE)
```

Tests cover legacy-mask parity at boundaries, overlapping reasons, unavailable
evidence, stale/tampered previews, retained metadata/extra assays, compact
history round-trip, empty selections, common plot scales, sample isolation and
combined-output invalidation through the actual app server.

Validation on 2026-11-24: the full regression run recorded 673 passing
expectations, with no failures, errors or skips. The final QC-focused run after
aligning sensitivity with the selected preview/applied view recorded 114
passing expectations, also without failures, errors or skips. These runs
overlap; their counts must not be added.

The existing GSE176031 pilot was reused without saving another expression
object. Its baseline remained 19,665 features by 595 cells, and the five
historical QC metrics matched by stable cell ID within `1e-12`. Default cutoffs
retained all 595 cells. A temporary upper-UMI cutoff exercised a 535-retained,
60-excluded preview and sensitivity export; it was a technical check, not a
biological recommendation or an applied/accepted decision. Joint and ribosomal
plots were rendered and visually inspected. Pilot/source metrics and both
bundled annotation checksums were unchanged; the original checkout's tracked
`site_data` still matched `main`. Browser hover/layout remains a manual check.

For manual review, move a cutoff, preview, inspect an excluded cell and switch
between preview/applied modes. Confirm that an unapplied preview leaves the
download and applied cell count unchanged. Switch samples, apply one, and verify
the other sample's count/history remains unchanged. Review the plot labels and
point hover behaviour in the browser; automated tests do not replace that check.

## Applied expression plot timing fix

The reported blank applied library-composition and highly-expressed-gene plots
were reproduced with nonempty applied data when Apply and plot clicks arrived
in one input batch. A plot request could record the old revision before Apply
created the new one. The stale-result guard then stopped rendering with an empty
validation message. Baseline requests were unaffected because their identity
did not change on Apply.

The controller now orders dependent observers explicitly: Preview at priority
20, Apply at 10, plot requests at 0. Baseline invalidation remains higher
priority. Missing, empty, unrequested and stale expression views have visible
messages. A changed applied result still requires a fresh plot click; this fix
does not reuse an image from an older sample or filtering decision.

`tests/testthat/test-qc-expression-plots.R` checks both real scater helpers,
first Apply/re-Apply, co-delivered Preview/Apply/plot actions, baseline requests,
sample changes and explanatory validation messages. Its synthetic data are
separate from user datasets; the installed scater's known internal
`perCellQCMetrics` deprecation warning is handled specifically in the tests.

Restart the development app before retesting. Apply filtering, then click both
applied expression buttons. Reapply a changed filter and request the plots
again; their data must match the current applied revision.

Fix validation: all 709 recorded expectations in the full regression suite
passed, with zero failures, errors or skips. Both applied expression images
also rendered through the actual app server for a temporary 424-cell pilot
subset using the originally failing input order. The source file checksum and
595-cell baseline were unchanged; no filtered source file or runner decision
was saved. The new tests reproduced missing images before the fix and passed
after it. The original checkout and protected user edits remained unchanged.

The work is versioned locally on `feat/shiny-input-integration`. Commit
`3043510` checkpoints the input/annotation/QC implementation before this fix.
The timing fix, regression test and these notes form a separate following
commit, so that commit can be undone with `git revert <fix-commit>` while
retaining the checkpoint. No push is required. User-owned decision-runner edits
and diagnostic scripts are kept outside these commits; no dataset or expression
checkpoint is included in them.
