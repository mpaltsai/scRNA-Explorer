# QC threshold recommender: development and evaluation protocol

Status: proposed protocol implementing accepted Q19–Q28 principles. No benchmark has been run and no improvement is established. Dataset accessions, numerical operating points, and final ranking rules are research deliverables, to be frozen before held-out evaluation.

## Question

Can an explainable recommender propose editable QC cutoffs that reduce false exclusions at comparable compromised-cell detection, or improve compromised-cell detection at comparable false exclusion, relative to established baselines? Evaluate rare populations separately. Retaining all cells is not a successful solution to this joint objective.

The output is a set of proposed thresholds and diagnostics, not an automatic truth label. Keep damaged-cell, doublet, and ambient-RNA evidence distinguishable; independently evaluate each endpoint. Preserve existing manual QC and before/after plots.

## Evidence and candidate methods

- Input: counts, detected genes, complexity, mitochondrial fraction when correctly annotated, capture/sample structure, species, chemistry, and available technical metadata. Document transformations and redundant metrics; avoid repeatedly weighting the same information.
- Additional technical metrics are optional inputs with a declared interpretation. Do not substitute zero for missing measurements. Validate reduced-evidence modes separately.
- Use cell-type annotations and stress/cell-cycle signatures to inspect disproportionate losses and uncertainty; initially exclude them from direct exclusion criteria. Do not select cutoffs to maximize downstream separation or differential-expression significance.
- Produce capture-specific suggestions. A shared fit needs a documented comparability rule across captures; globally poor captures require explicit reporting rather than automatic threshold relaxation.
- Benchmark manual/fixed cutoffs and robust univariate QC as simple comparators, miQC as a joint mitochondrial/complexity comparator where applicable, and SampleQC as a multivariate/multi-sample comparator where applicable. Methods must only be assessed on inputs satisfying their requirements.

miQC jointly models mitochondrial fraction and detected genes; SampleQC models heterogeneous QC distributions across samples. Neither a generic outlier score nor agreement among methods is an independently validated poor-quality label. Method selection and any novel estimator must be justified by observed gaps, not assumed superiority. See the primary [miQC paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009290) and [SampleQC paper](https://doi.org/10.1186/s13059-023-02859-3).

## Dataset inventory before experiments

For every candidate dataset record accession/source, license/access conditions, tissue, species, protocol, captures, independent biological units, number of cells, original filtering, raw-count availability, feature identifiers, and available reference evidence. Define eligibility and split assignments before fitting. Do not download or execute private study objects merely because their paths appear in reference notebooks.

Include heterogeneous tumor and non-tumor settings, heterogeneous cell populations, and captures of different sizes. PCa notebooks can inform scenarios but do not establish independent QC truth. Resolve their raw-versus-filtered provenance before using them as test material.

Distinguish reference evidence:

- Independent cell-quality measurements, when available: document the measurement, reliability, missingness, and relation to the targeted failure mechanism.
- Expert review: blind reviewers to method suggestions where feasible, retain disagreements, and describe labels as expert judgments rather than unquestionable truth.
- Cell-type labels: useful for population-specific loss assessment; not evidence that all labeled cells are viable or all unassigned cells are poor quality.
- Controlled count perturbations: useful stress tests with known interventions; do not claim these fully reproduce real damaged-cell biology. Use more than one perturbation mechanism and vary severity.

## Separation of development and evaluation

Split by independent donor/study as appropriate to prevent related captures or repeated cells leaking between development and held-out evaluation. Tune candidate generation, grouping rules, ranking, missing-data modes, and abstention diagnostics on development data only. Evaluate transfer across studies/protocols when the intended supported scope requires it.

Freeze before held-out evaluation: dataset eligibility, splits, software/reference versions, seeds, candidate-method settings, comparator operating points, ranking rule, metrics, uncertainty procedure, failure/abstention criteria, and numerical success criteria. Research may use development cross-validation; held-out results must not become another tuning loop.

## Candidate cutoffs and recommendation behavior

1. Validate input evidence and report unavailable or ambiguous metrics.
2. Fit applicable candidate models; diagnose failed/unstable fits and capture comparability.
3. Generate a bounded set of per-metric cutoff combinations. Record which joint-model patterns independent cutoffs cannot reproduce.
4. Compare candidate cell sets, reasons, model disagreements, and sensitivity to threshold/fit perturbation. Protect uncertainty from being converted into a confident exclusion claim.
5. Present a conservative candidate with alternatives and retained/proposed-excluded/review statuses. These are proposed decisions until the user applies filtering.

The exact generation/ranking rule is an empirical output. Never infer candidate accuracy from agreement with the same model that generated it. Do not rank solely by retention percentage, cleaner UMAP, marker separation, or more significant downstream tests. If evidence does not support a recommendation, abstain and preserve manual filtering.

## Evaluation and reporting

Where independent quality evidence exists, report false exclusion among reference acceptable cells and detection of reference compromised cells at matched operating points. State denominators and label uncertainty. Show capture/donor and population-specific results rather than relying only on pooled-cell averages. Assess uncertainty across independent biological units/studies, not by treating every cell as an independent replicate.

Across all eligible datasets report candidate retention, population-specific losses, threshold and cell-set sensitivity, fit failures, abstention frequency, missing-evidence coverage, runtime, peak memory, and sensitivity to seeds. Stability is a useful property, not proof of correctness. Report conditional performance alongside coverage so abstaining on difficult cases cannot hide failure.

Investigate failure cases: globally poor captures, rare plausible populations, high-mitochondrial but complex cells, low-RNA populations, missing/mismatched feature annotation, small captures, different chemistries, correlated QC metrics, and competing damage/doublet/ambient evidence. Do not expand supported scope on the strength of one successful dataset.

## Decision gates

1. **Specification gate:** dataset inventory, independent-evidence definitions, frozen split/tuning plan, numerical acceptance criteria, and runtime budget exist before held-out testing.
2. **Research gate:** held-out results support an accepted trade-off against baselines with uncertainty reported and no concealed population-specific failures. Otherwise revise using development evidence or retain experimental/manual-only behavior.
3. **Integration gate:** existing manual filtering and plots remain usable; suggested cutoffs reproduce the previewed cell set when applied; lineage, missing-data explanations, abstention, invalidation, and exports are verified.
4. **Hosted gate:** representative single-session and concurrency benchmarks support configured limits. Failures do not corrupt another session or silently produce incomplete recommendations.

Produce an executable benchmark configuration, dataset manifest, frozen development/held-out assignments, method/reference versions, tables/figures, a failure-case report, and a justified recommendation about supported modes. None of these empirical results exists yet.
