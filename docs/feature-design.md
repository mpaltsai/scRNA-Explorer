# Optional analysis extensions

Working design from the feature-planning interview. Accepted scope is recorded separately from unresolved behavior; implementation has not been authorized by this document.

## Accepted direction

- Preserve all existing routines and add new capabilities alongside them.
- Initially target the new preprocessing workflow at droplet-based, UMI scRNA-seq. This does not remove any existing input support.
- Provide guided defaults and expandable expert controls.
- Extend the existing before/after filtering plots to support review of new QC and doublet decisions before exclusions are applied.
- Preserve original counts and record settings and exclusion reasons for reproducibility.
- Prioritize annotation and subclustering, followed by condition comparisons.
- Reuse compatible existing analysis objects, helpers, and fixtures across tickets. Do not create fresh objects or duplicate full matrices solely for ticket boundaries; create derived data only when required by the agreed feature and preserve original data.

## Confirmed subclustering behavior (Q7)

- Users select one or more parent clusters for a separate subanalysis; the parent analysis remains available.
- Recompute variable-feature selection, PCA, neighbors, clustering, and UMAP within the selected population.
- Expose dimension and resolution controls; the notebook's settings are not universal defaults.
- Repeating integration is an explicit user choice, not an automatic step for every subset.
- Show sample composition alongside the subclusters and retain the relationship to parent clusters.

Normalization/scaling and integration choices were subsequently settled in Q8–Q9. Label propagation was settled in Q5; downstream rerun behavior in Q11.

## Confirmed workflow choices (Q1–Q6)

- Q1: Existing defaults remain unchanged; new preprocessing steps are opt-in.
- Q2: Add an editable/uploadable sample-information table and optional barcode-level metadata for multiplexed captures. Require fields only when the selected analysis needs them. Do not infer biological replication from filenames.
- Q3: Doublet detection defaults to flagging. Exclusion is a separate user action after review, with scores, labels, exclusion reasons, and original data retained.
- Q4: Preserve enrichment-based annotation and add manual annotation and SingleR suggestions. Store automatic suggestions separately from user-approved labels and support unassigned populations.
- Q5: Applying subcluster labels to the parent dataset requires an explicit action; preserve original cluster IDs and annotations.
- Q6: Treat Merged as reference material for later correlation/enrichment extensions. Annotation and subclustering remain first; classification is deferred.

## Confirmed analysis behavior (Q8–Q12)

- Q8: Each subanalysis recomputes normalization, variable genes, scaling, and PCA on its selected cells using the chosen count matrix. Initially use the existing log-normalization workflow, with regression off by default. Preserve original counts separately when ambient correction is applied.
- Q9: Offer no integration, CCA, and RPCA for subclustering. Preselect no integration; require an explicit grouping choice when integration is enabled. Do not automatically use condition as the grouping variable. Show sample/batch and biological labels for review. Existing parent-workflow defaults remain unchanged under Q1.
- Q10: Require explicit selection of a species-compatible curated annotation reference and show tissue coverage and limitations. Keep marker evidence visible. Custom-reference uploads are deferred.
- Q11: Changes to an upstream analysis step mark affected downstream results out of date and require an explicit rerun. Keep previous results inspectable without presenting them as current. Annotation renaming updates labels without unnecessary numerical recomputation.
- Q12: Add an analysis-bundle download alongside existing outputs: analysis object, metadata, settings, seeds, package versions, filtering decisions, annotation history, and parent/subcluster relationships. Q17 below defines the initial restoration scope.

## Confirmed QC and deployment choices (Q13, Q18)

- Q13 (revised by the user): Develop a method to recommend QC thresholds using all available relevant information. A fixed three-MAD rule is no longer the accepted recommendation engine. Use and preserve scRNA-Explorer's existing QC plots, including before/after filtering displays. Applying exclusions still requires review. Q19–Q22 refine this requirement below.
- Q18: Target the existing Shiny deployment first. Provide progress reporting and configurable per-session resource limits, establish capacity through representative benchmarks, and preserve sparse matrices in new workflows. Do not promise a cell-count limit before benchmarking.

## Confirmed doublet, ambient, annotation, and restoration choices (Q14–Q17)

- Q14: Require capture identity and review of expected doublet rate. Offer automatic rate estimation only for supported experimental setups; expose expected rate and uncertainty in expert controls and record estimated versus user-supplied rates. Do not silently transfer chemistry-specific assumptions to all datasets.
- Q15: SoupX requires matching raw-droplet and called-cell matrices. Mark clusters used for estimation as preliminary. Retain original and corrected counts separately; accepting corrected counts invalidates dependent results. Missing raw droplets disables SoupX without blocking existing workflows. Exact ordering relative to doublet review remains to be specified.
- Q16: Offer per-cell SingleR annotation, preselected, and a separate cluster-level expert option. For per-cell annotation, summarize predicted-label composition and ambiguous assignments by cluster alongside marker evidence. Do not silently name heterogeneous clusters by majority prediction. Keep original cluster IDs; accepting or editing biological labels remains explicit. Diagnostic scores must not be presented as calibrated probabilities of correctness. Tumor/normal and ERG classifiers remain deferred.
- Q17: Initially restore supported, versioned app-generated bundles. Preserve existing count-import routines. Validate compatibility and report unsupported versions; general arbitrary processed-Seurat import is deferred. Restore saved results, not interrupted computations; recorded software versions alone do not guarantee identical reruns across environments.

## Confirmed QC-method principles (Q19–Q22)

- Q19: Develop and validate a recommendation workflow with established methods as foundations and benchmarks. Pursue a new statistical estimator where a demonstrated gap justifies it; do not assume a new method is better.
- Q20: Initial evidence includes counts, detected genes, complexity, mitochondrial fraction, sample/capture structure, species, chemistry, and available technical metadata, with additional metrics when available. Account for correlated metrics. Keep doublet and ambient-RNA evidence distinguishable from damaged-cell evidence. Initially use cell-type, stress, and cell-cycle information for diagnostic review rather than direct exclusion criteria.
- Q21: Editable per-metric thresholds remain the primary output. Explicitly show disagreement with joint-model decisions in the existing QC plots, distinguishing retained/excluded/review statuses. Do not imply independent cutoffs exactly reproduce a multivariate decision boundary.
- Q22: Prioritize detecting genuinely compromised cells while preserving valid and rare populations; allow insufficient-evidence outcomes. Validate using independent quality evidence where available, controlled simulations, and held-out biological samples, with threshold stability and population-specific losses assessed. Cleaner UMAPs or stronger differential expression alone are not validation of better QC.

Model agreement is evidence, not ground truth, and a generic outlier score is not a calibrated probability of poor quality.

## Confirmed QC recommendation and evaluation choices (Q23–Q25)

- Q23: Produce capture-specific recommendations. Share information across captures only when relevant tissue, chemistry, and other characteristics support a shared model. Display globally poor captures separately; do not respond merely by relaxing thresholds.
- Q24: Present a small set of candidate threshold combinations, one conservative recommendation, and alternatives showing filtering trade-offs. For each candidate show proposed excluded cells, reasons, joint-model disagreements, and sensitivity to threshold changes. Use population annotations diagnostically to highlight disproportionate losses. Do not target a predetermined retention percentage. Select the ranking rule through benchmarking rather than arbitrary weights.
- Q25: Require evidence of fewer false exclusions at comparable compromised-cell detection, or better detection at comparable false exclusion, against established baselines on held-out datasets. Assess rare populations separately and include heterogeneous tumor data. When independent quality labels are unavailable, report sensitivity/stability without claiming measured accuracy. Keep the recommender experimental until it passes the agreed evaluation.

## Confirmed prerequisite and sequencing choices (Q26–Q28)

- Q26: Correct barcode handling and metadata preservation as necessary for the new features, while preserving existing analytical routines and defaults. Verify equivalent legacy results where these corrections should not alter analysis.
- Q27: Use supported, valid evidence only. Missing mitochondrial annotation is not zero mitochondrial expression. For small captures, failed fits, or missing metrics, use a reduced-evidence mode only after validating that mode; otherwise abstain from recommendations and preserve manual filtering with an explanation.
- Q28: Develop the recommender as an experimental research component while established additions progress independently. Start with identity/metadata preservation, then doublet review, annotation, and subclustering. Benchmark the recommender alongside those additions and integrate with existing QC plots without premature validation claims. Ambient correction and bundle restoration follow prerequisites; condition comparisons remain a subsequent phase.

## Accepted recommendations to refine

- Adaptive QC alongside manual thresholds, with sample/technical-batch grouping and review of proposed thresholds.
- scDblFinder as the initial doublet method, operating by capture before integration, retaining scores and labels.
- Optional SoupX when suitable raw-droplet and called-cell inputs are available; account for its preliminary-clustering requirement.
- Cell-cycle scoring with regression as an explicit choice; preserve the signal by default.
- Improve existing enrichment-based annotation with marker review, editable labels, explicit uncertainty, and reference-based suggestions.
- Add subclustering while retaining the relationship to the original clusters.
- Add replicate-aware pseudobulk condition comparisons, distinct from descriptive cluster marker discovery.
- Extend existing per-sample cell-type proportion summaries before introducing formal differential-abundance methods.
- Defer trajectory inference and cell-cell communication until a specific biological question justifies them.

## Research outputs and later-phase specifications

- Detailed metadata schema and validation for capture, biological sample, donor, condition, and technical batch; table import/editing is confirmed.
- QC recommendation method: empirically select candidate generators and ranking rules, and specify benchmark datasets/numerical acceptance criteria before held-out evaluation. Missing-data/fit policy is confirmed in Q27; evidence and candidate presentation in Q19–Q25.
- Initial curated annotation reference list and detailed diagnostic presentation; annotation granularity and separate suggestions/manual labels are confirmed.
- Optional ambient-RNA workflow ordering relative to doublet detection and filtering; input requirements and corrected-count lineage are confirmed.
- Supported bundle-version combinations and benchmark-derived compute limits; restoration scope, deployment target/progress reporting, export contents, and downstream invalidation are confirmed.
- Exact implementation contracts and version support require fixture-based verification; accepted input-correction scope is recorded in Q26.
- Replicate eligibility, paired designs, and covariates for condition comparisons.

The current planning deliverables are [the implementation plan](implementation-plan.md) and [the QC benchmark protocol](qc-benchmark-protocol.md). Empirical choices above are explicit research deliverables, not silently selected methods. Later-phase condition-comparison and classification designs are not complete. The user confirmed shared understanding after review of these deliverables and requested a new GitHub branch for the work.

## Additional work supplied by the user

Reference directory: `/home/ismini/Documents/PCa/Merged`. Classification is deferred to a later discussion.

Source inspection found downstream analyses that load prepared Seurat objects; it did not find the upstream preprocessing, annotation, or clustering implementation in the inspected R/Rmd sources. Candidate additions, not yet accepted scope:

- Repeated-subsampling correlation with recurrence and stability summaries (`correlation_results/split/testing/cor_sampling_script_delayed_wilcoxon.R`).
- Correlation-set overlap and sensitivity to sample size (`correlation_results/split/testing/sampling/cor_multiple_samples.Rmd`).
- Comparison of correlated-gene sets with differential-expression results, module scoring, secreted-gene annotation, and gene-pathway networks (`enrichment_results/enrichment_results_Merged_Split2.Rmd`).

Before reuse, distinguish study-specific target/cohort choices from general controls, define the interpretation of zero expression, and review the statistical meaning of subsampling recurrence and summarized adjusted p-values. A subsampling stability measure must not be presented as replication across independent biological samples. These scripts are reference material, not validated app modules or authorization to migrate them.

## Upstream notebook supplied by the user

Located at `/home/ismini/Documents/PCa/PCa_analysis_merged_Split_new_TN_classifier_new_ERG_classifier.Rmd` (the supplied path omitted a directory separator). Classification sections remain deferred and were not reviewed.

- The merged notebook starts from prepared epithelial Seurat objects, preserves dataset prefixes, and integrates layers split by dataset-prefixed patient identifiers.
- Active workflow: RNA layer preparation, normalization, variable-feature selection, scaling, PCA, CCA integration, neighbors using dimensions 1:30, clustering at resolution 2.5, and UMAP. Most underlying method parameters are implicit package defaults, not pinned settings.
- Marker-signature module scores provide additional annotation evidence.
- Linked Song and Chen notebooks demonstrate epithelial cluster selection followed by repeated normalization, feature selection, scaling, PCA, CCA, clustering, and UMAP. Their fixed cluster IDs and resolution are study-specific, not proposed app defaults.
- Linked QC scripts use different study-specific thresholds and contain commented alternatives. No doublet or ambient correction was found in the inspected QC source.
- Chen source has a provenance ambiguity: an analysis loads an object named raw despite prose describing filtered data. Resolve the actual filtering lineage before treating it as a validated reference workflow.

These observations inform pending decisions about configurable integration, signature scoring, and subclustering; they do not settle those decisions or authorize code migration.

## Code observations requiring attention during implementation planning

- Normalization recreates Seurat objects from counts without carrying over all input cell metadata (`site_data/R/normalize_data.R`). Stable cell identifiers and metadata preservation are prerequisites for these extensions.
- The single-MTX input reader replaces supplied barcodes with numeric identifiers (`site_data/R/my_read10X.R`). Raw/called-cell matrix pairing and cell-metadata import require reliable barcode matching.
- An initial UMI filter precedes the existing dataset called "unfiltered" (`site_data/R/init_filter.R`). It must not be treated as a raw droplet matrix.
- Existing annotation selects the first enrichment result as a label (`site_data/R/annotateMyClusters.R`). New annotation choices must preserve access to this routine.

These observations are read-only findings, not authorization to change existing routines.

## Method references

- [miQC: joint modeling of mitochondrial fraction and detected genes](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009290)
- [SampleQC: robust multivariate, multi-cell-type, multi-sample QC](https://doi.org/10.1186/s13059-023-02859-3)
- [Adaptive QC and its assumptions](https://bioconductor.org/books/release/OSCA.advanced/quality-control-redux.html)
- [scDblFinder](https://bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html)
- [SoupX](https://github.com/constantAmateur/SoupX)
- [Cell-cycle scoring and regression](https://satijalab.org/seurat/articles/cell_cycle_vignette)
- [SingleR annotation diagnostics](https://bioconductor.org/books/release/SingleRBook/annotation-diagnostics.html)
- [Seurat differential expression and pseudobulk](https://satijalab.org/seurat/articles/de_vignette.html)
