# Review one analysis step at a time

The standalone runner keeps the analysis code separate from Shiny. Use an existing in-memory SingleCellExperiment or Seurat object. A method produces a **step proposal**; inspect its object, diagnostics, and files before explicitly accepting or rejecting it with a reason. The runner records choices and evidence; it does not decide which scientific setting is optimal.

The first version is tested at four public boundaries: propose and inspect; accept/reject with a reason; reload history; and block downstream execution until its parent result is accepted. Revisions, failures, and changed evidence are checked through those same boundaries.

## Start with your existing object

From an R session at the repository root:

```r
source("analysis/R/analysis_runner.R")
source("analysis/examples/qc_metrics.R")

# Load once, only if the object is not already in your R session.
input_path <- "analysis/local/GSE176031-r1/pilot_sce.rds"
sce <- readRDS(input_path)
directory <- "analysis/local/GSE176031-review"  # new, dedicated log directory

proposal <- propose_analysis_step(
  sce, directory, step = "inspect-qc-metrics", execute = qc_metric_step,
  parameters = list(plot_path = tempfile("qc-", tmpdir = directory, fileext = ".pdf")),
  seed = 1L, input_path = input_path,
  source_files = c("analysis/examples/qc_metrics.R", "analysis/R/analysis_runner.R")
)
proposal$diagnostics
proposal$artifacts    # open the PDF at the recorded path
head(as.data.frame(SummarizedExperiment::colData(proposal$object)))
history <- read_analysis_history(directory)
history[, c("step", "status", "reason", "current")]
```

This example adds UMI/gene-count metrics and writes three diagnostic plots. It changes no counts, filters no cells, and proposes no cutoffs. It is an example of the runner contract, not a replacement for the existing app's pre/post-filtering QC plots. Those routines remain available in the app and are to be retained during future integration.

The baseline path and its checksum identify your declared existing checkpoint. The first object passed must be the representation you intend to use; the runner fingerprints that in-memory baseline without loading a second copy of the checkpoint to establish equality. Reuse the same checkpoint and baseline for revisions. File checksums detect later changes to that checkpoint.

## Record your decision

After inspecting the evidence, choose **one**:

```r
# Acceptance returns the reviewed object with its accepted-run marker.
accepted <- decide_analysis_step(proposal, "accept",
  reason = "Reviewed cell identities, unchanged counts, and QC metric distributions.")

# Or reject it; this returns NULL. Keep the previous accepted object/baseline.
decide_analysis_step(proposal, "reject",
  reason = "The method needs correction before continuing.")
```

These are example reasons, not recorded scientific approvals. There is no default acceptance. Decisions are immutable. To correct a choice, execute a fresh proposal and explain the revision instead of editing old records.

A later method uses the accepted object and its explicit parent:

```r
# next_method(object, parameters) is the next R analysis function you develop.
next_proposal <- propose_analysis_step(
  accepted, directory, "next-step", execute = next_method,
  parameters = list(), parent_run = proposal$run_id,
  source_files = "path/to/next_method.R"
)
```

The runner refuses pending/rejected parents, changed parent objects, and parents outside the current accepted history. It also refuses acceptance if the object, diagnostics, referenced files, or baseline checkpoint have changed. If another proposal was accepted in the meantime, regenerate your proposal against the current state before deciding. A rejected proposal may still be recorded with its reason after the history advances.

## Correct a method or parameter

Keep or reload the appropriate earlier accepted object (or the original baseline). Call `propose_analysis_step()` with corrected code/settings and that accepted parent. For a revision from the baseline, pass the original baseline object with `parent_run = NULL`. Merely generating a revision changes nothing in the accepted history.

Accepting the revision makes it the current result. Earlier choices on the replaced branch and their downstream results appear as `superseded`; their decisions and reasons remain intact. Previously accepted ancestors shared by the revision remain accepted. Downstream results must be rerun from the revised accepted object. Rejected and failed attempts remain visible.

## Method contract

`execute(object, parameters)` returns an SCE/Seurat object, or:

```r
list(object = updated_object,
     diagnostics = list(summary = "compact review evidence"),
     artifacts = c("path/to/plot.pdf", "path/to/table.csv"))
```

Only `object`, `diagnostics`, and `artifacts` are accepted return fields. Diagnostics and named parameters must contain plain data, not expression matrices, S4 objects, functions, or environments. Keep them small. Artifacts are referenced by absolute path and checksum; the runner neither copies them nor silently regenerates them. Methods must use new paths for evidence they write. A failed method is logged and cannot be accepted; its external side effects are not rolled back.

Output cell identifiers must be unique and retain input identities; subsets/reordering are allowed. Original `geo_` cell fields and the full `geo` provenance record are preserved across SCE/Seurat conversion and checked for conflicts. Methods remain responsible for retaining other metadata, choosing assays/layers explicitly, and invalidating/recomputing stale derived results. This runner does not yet enforce a complete dependency graph of assays, reductions, and annotations within an object.

## What is saved, and how to resume

The dedicated log contains metadata-only RDS records for the baseline reference, proposals, successful diagnostics/evidence, failures, and decisions. `read_analysis_history(directory)` reloads those records without loading expression data. It returns status, parent, reason, current-result flag, and list columns for parameters, diagnostics, artifacts, code, and input provenance. `current` marks the latest accepted result; `accepted` also includes its accepted ancestors.

For each proposal the log captures the R function text, declared source-file contents/checksums, git revision/dirty flag when available, R/RNG settings, seed, and loaded package versions after successful execution. Pass settings explicitly in `parameters` and list every relevant helper in `source_files`: closure variables, undeclared helpers, external reference files/services, and the full git diff are not automatically captured. A local file snapshot aids review; this is not an environment lockfile or a promise of bitwise reproducibility. The caller's random state is restored on success and failure.

Full expression objects are **not** automatically checkpointed. Keep your existing objects in the R session. If you deliberately save an accepted object to an existing analysis checkpoint location, reload it later and use its `metadata(sce)$analysis_runner$run_id` or `seu@misc$analysis_runner$run_id` as the parent. Without a saved result, a lost in-memory proposal must be rerun; the log alone cannot reconstruct expression data or the candidate object. Choose checkpoints deliberately, rather than creating one per ticket or decision. Never overwrite the referenced baseline checkpoint.

Fingerprints traverse the in-memory object and can cost time/memory on large datasets. Initial support is for ordinary in-memory SCE/Seurat objects; file-backed assays/external pointers and their backing data are not qualified by this runner. Logs and referenced evidence must remain at their recorded paths; log relocation is not supported yet. Keep logs private when metadata are sensitive. Concurrent decisions fail safely on the directory lock. A leftover `.analysis-lock` requires inspecting the writer before retrying; the runner does not automatically remove another session's lock.
