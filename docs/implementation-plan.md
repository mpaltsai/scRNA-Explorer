# scRNA-Explorer extension plan

Status: accepted first-phase scope and sequence, based on decisions Q1–Q28 in [feature-design.md](feature-design.md) and the user's confirmation of shared understanding. Manage changes on the dedicated `feat/qc-and-downstream-extensions` branch. Dependency installation/upgrades and deployment require separate approval; the user requested publication of the planning documents on a new GitHub branch.

The approved breakdown is published as GitHub issues [#1](https://github.com/mpaltsai/scRNA-Explorer/issues/1) through [#27](https://github.com/mpaltsai/scRNA-Explorer/issues/27). Native GitHub blocking relationships determine the available work. The initial independent starting tickets are cell-identity preservation ([#1](https://github.com/mpaltsai/scRNA-Explorer/issues/1)) and the QC dataset/evidence manifest ([#20](https://github.com/mpaltsai/scRNA-Explorer/issues/20)).

## Outcome and compatibility

Development sequence clarified by the user: implement and validate features in standalone R workflows against several datasets first; only then integrate validated features into Shiny. Extending versus rebuilding the Shiny app remains undecided. Ticket #28 uses the supplied DGE files and GEO metadata, retains hg19 as user-confirmed while preserving original source declarations, and provides metadata attachment to existing SingleCellExperiment/Seurat objects. AnnData is not an input to this implementation; current Shiny scripts remain unchanged. See [standalone R usage](../analysis/README.md).

Implement features, analysis helpers, validation scripts, tests, and benchmarks in R, compatible with the Seurat/SingleCellExperiment workflow. Reuse compatible existing R objects and shared readers. Use the supplied DGE inputs for GSE176031; other datasets may enter through existing supported input formats. Do not add a Python runtime or require Python-based conversion. Record optional R artifact incompatibility instead of rebuilding a full dataset for every ticket; dependency changes still require approval.

Extend the app with optional preprocessing, annotation, and subclustering workflows. Preserve current analytical routines, defaults, plots, and downloads. Users who do not enable new options retain the existing workflow. Necessary cell-identity and metadata corrections are accepted, with legacy-result equivalence checked where no analytical difference is expected.

The threshold recommender remains an opt-in experimental research component until it passes a frozen benchmark protocol. Its research outcomes do not block established additions. Preserve original counts and make exclusion, corrected-count selection, and annotation acceptance explicit actions.

Reuse compatible existing analysis objects, helpers, and fixtures across tickets. Completing a ticket must not require constructing a fresh analysis object or duplicating full expression matrices merely to cross a ticket boundary. Store new metadata/results on the existing analysis where compatible. Materialize derived data only when the feature requires it, such as an explicit subanalysis or accepted corrected counts, while preserving original data. Reuse suitable existing study objects only with established provenance and authorized access; their availability does not establish independent QC truth.

## Verified repository constraints

- Dependencies are loaded in `site_data/R/load_libraries.R` and individual source files, without a tracked package manifest/version lock. Availability and compatibility of scDblFinder, SingleR, SoupX, and research comparators must be checked before proposing installations. The current deployed R/Bioconductor/Seurat/SeuratObject versions are unknown.
- No tracked automated test harness or CI configuration was found; `site_data/R/test.R` is commented material. Establish focused regression fixtures and checks rather than assuming an existing suite. Do not execute private-study diagnosis scripts as a substitute for fixtures.
- The README names a hosted Shiny app and Docker image, but no tracked deployment configuration establishes hardware, timeouts, or concurrency limits. The upload option in `site_data/server.R` is not a capacity benchmark.
- ZIP uploads currently reuse `./userFile` in `site_data/server.R`; session-isolated storage is a prerequisite for safely testing new multi-session workflows.
- Existing multisample RDS downloads are not reimportable through the current input flow, as noted in `site_data/tab_1.R`. Bundle restoration is new functionality.

These are source-inspection findings. No application was executed, private expression objects loaded, or software installed during planning.

## Milestone 1: identity, metadata, and reproducibility foundations

1. Preserve uploaded cell IDs/barcodes and create globally unique internal IDs with reversible mappings to their source capture and input identifiers. Detect collisions instead of silently matching unrelated cells. Non-barcode inputs retain their original identifiers with source context.
2. Support editable sample-information tables and optional barcode-level metadata. Distinguish capture, biological sample, donor, condition, and technical batch. Retain existing sample fields for compatibility. Require only the fields needed by each operation and explain unresolved assignments; multiplexed captures must support multiple donors.
3. Preserve metadata through normalization, merging, integration, subsetting, and export. Maintain sparse representations in new workflows and address shared dense conversions where necessary for the accepted extensions.
4. Distinguish original called-cell counts, optional raw-droplet matrices, and optional corrected counts. Do not relabel the current initial-filtered input as raw droplets. Record input identity, selected cells, parameters, seeds, and software/reference versions for each new analysis.
5. Keep parent and subanalysis results distinct. Changes invalidate affected downstream outputs; previous outputs remain inspectable as out of date. Label-only edits update views without unrelated numerical reruns.

Validation: use small representative single- and multiple-sample fixtures, repeated barcodes across captures, mixed donor/capture mappings, and input/export identity round trips. Compare unchanged legacy numerical outputs with an established baseline; document intentional differences from corrections. Confirm session-specific uploads and results cannot collide before shared hosting tests.

## Milestone 2: established additions

### Doublet review

Add optional scDblFinder before integration, grouped by the appropriate capture structure. Record score, predicted class, method settings, chemistry assumptions, and expected-rate provenance. Default to flagging; applying exclusions requires review. Offer automatic rate estimation only for experimentally supported setups, with expert rate/uncertainty controls. Display QC relationships and results per capture; do not interpret every high-count cell as a doublet.

### Annotation

Preserve existing enrichment-based annotation. Add editable user labels and explicit unassigned status. Add SingleR using an explicitly selected species-compatible curated reference; present tissue coverage, identifier compatibility, and reference version. The initial reference catalog is a documented implementation deliverable; custom-reference uploads are deferred.

Per-cell SingleR is preselected; cluster-level prediction is an explicit expert alternative. Summarize cell predictions and ambiguity within each cluster alongside marker evidence. Keep original cluster ID, automatic suggestions, and accepted labels separately. Do not silently replace heterogeneous clusters with a majority label or interpret scores as correctness probabilities. Annotation-service failure must not block basic cluster inspection in the new workflow.

### Subclustering

Select one or more parent clusters; create a separate subanalysis with parent-cell membership retained. Recompute log normalization, variable features, scaling, PCA, neighbors, clustering, and UMAP using selected cells and the chosen count matrix. Regression starts off. Offer no integration (preselected), CCA, and RPCA; integration grouping is explicit and must not default to condition. Expose dimensions and resolution with validity checks for the selected data. Show sample/batch and biological composition. Applying subcluster labels to the parent is explicit and preserves prior IDs/labels.

Validation: flags do not remove cells, explicit exclusions remove exactly the reviewed cells, annotation modes retain separate outputs, parent/subcluster membership is exact, no-integration paths work, and each integrated path is exercised on compatible multi-sample fixtures. Verify rerun/invalidation behavior against concrete upstream changes.

## Parallel research milestone: threshold recommendations

Follow [qc-benchmark-protocol.md](qc-benchmark-protocol.md). Joint-model evidence informs candidate editable cutoffs; the existing scRNA-Explorer QC plots show threshold effects, reasons, disagreement, and sensitivity. Capture-specific recommendations can share information only across comparable captures. Biological interpretations remain diagnostic rather than direct exclusion criteria, as recorded in [ADR 0001](adr/0001-separate-qc-evidence-from-biological-interpretation.md).

Offer a small candidate set with a conservative recommendation and explicit review/insufficient-evidence outcomes. Missing mitochondrial annotation is unavailable evidence, not a biological zero. Reduced-evidence modes need their own validation. Do not optimize for a fixed retention percentage or claim model agreement establishes true cell quality. Existing manual filtering remains available.

## Milestone 3: dependent optional workflows

- SoupX: validate matching raw-droplet and called-cell inputs; distinguish preliminary estimation clusters from final biological clusters. Preserve corrected and original counts, show effects, and invalidate downstream outputs when users accept corrected counts. The precise sequence with doublet detection/final filtering is a method-validation deliverable; record the executed order and require dependent reruns rather than silently reusing incompatible results.
- Cell-cycle scoring: optional scoring with regression only by explicit choice; scores alone do not constitute a QC failure.
- Analysis bundle: export objects, metadata, settings, seeds, software/reference versions, filtering decisions, annotation history, and parent/subcluster lineage. Restore supported versioned app bundles after compatibility validation. Existing count imports/downloads remain. Arbitrary processed-Seurat import and resuming interrupted computations are outside initial restoration scope.
- Hosted capacity: progress reporting and configurable per-session limits, established by representative runtime/memory/concurrency benchmarks. No promised cell-count capacity before measurement. Remote host capacity is not inferred from the development workstation.

## Subsequent analysis phase

Design replicate-aware pseudobulk condition comparisons after annotation/subclustering. Before implementation, specify independent biological replication, repeated/paired designs, covariates, count source, contrasts, cell-number eligibility, and multiple-testing policy. Preserve descriptive marker analysis separately. Existing proportion summaries precede formal differential-abundance methods.

Later correlation/enrichment extensions can draw on the Merged scripts after reviewing their statistical assumptions. Classification, trajectory inference, and cell-cell communication remain deferred. This plan does not claim their detailed designs are complete.

## Review checkpoint

The user confirmed the scope and sequence after reviewing the planning deliverables. Implementation can resolve routine details within accepted decisions; empirical method choices require documented evidence. Installation/upgrades/removal of dependencies and external publication outside the requested planning-branch push require separate user approval with scope and rollback information.
