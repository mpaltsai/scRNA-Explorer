# Local annotation references

Reading counts and preparing QC use local annotations. They never contact
AnnotationHub. The bundled human and mouse `.Rda` files remain unchanged.
Their AnnotationHub origin is confirmed by the user; their exact resource IDs,
releases, snapshots and builds could not be recovered. Read-only comparison
against the locally cached human Ensembl 99/113 and mouse Ensembl 99 databases
found no exact match in the exported gene fields. These are not assigned as
the exports' sources. `site_data/annotation-provenance.csv` records this
distinction and binds the known information to each export's MD5.

In the development app's input tab:

1. Choose species and, when known, the counts' genome build. The build is a user
   declaration, not an inferred or independently verified property of the counts.
2. The bundled export remains the default. Its unknown build remains visible.
3. Expand **Annotation reference maintenance**, then **Check for reference
   updates**. Only this action checks online. It may download AnnotationHub's
   catalogue database; it does not download an EnsDb annotation resource.
4. Review the candidates. Only Ensembl EnsDb or GTF-derived gene annotations
   for the selected species and known genome build are listed.
   `hg19`/`GRCh37`, `hg38`/`GRCh38`,
   `mm10`/`GRCm38` and `mm39`/`GRCm39` are treated as build-name aliases.
   This does not perform coordinate conversion or guarantee matching annotation
   releases, gene models or identifier coverage. No compatible candidate is a
   valid outcome. Do not substitute GRCh38 for hg19 because it is newer.
5. Select a candidate, give a reason and explicitly confirm its download.
   A resource may require hundreds of MB. It is saved separately, never over the
   bundled export. Its own metadata must agree with the candidate's species and
   build. For EnsDb resources the database metadata are checked independently.
   GTF-derived GRanges must contain gene IDs, gene names and gene biotypes;
   any sequence genome recorded in the object must match the catalogue build.
   When the sequence genome is absent, the provenance explicitly records that
   the build is supported by catalogue metadata only. Release, resource ID,
   snapshot, source metadata, retrieval time, reason, table hash and file
   checksum are recorded.
6. Downloading does not select a reference or change prepared QC. To use it,
   choose it under **Annotation reference**, then **Add gene annotations** and
   review the mapping/QC again. A downloaded reference requires a declared,
   matching counts build. Changing the selected reference or declared build
   clears preparation, the baseline and affected downstream results through the
   existing invalidation path. Recreate the per-sample baseline before filtering.

The QC configuration in each SCE records the selected reference and hashes;
saved analyses retain that information. Existing imported metadata and counts
are preserved. There is no automatic reannotation of a saved analysis. Checksums
are verified when loading a versioned reference. Missing, corrupted or changed
references stop loading instead of silently falling back to another version.

## Monthly checks and storage

A due notice appears one calendar month after the last successful check for
that species/build, clamping to the last day of a shorter month. Failed checks
record the attempt/error but do not advance the last-success date. Previously
successful candidate results remain labelled stale after a failed refresh.
The user can force a check sooner with the button.

This is an app notice and an explicit check, not an installed cron job or a
background scheduler. A closed app does not perform checks. The standalone R
function below checks only when due unless `force=TRUE` is supplied.

Catalogues are scoped to the installed AnnotationHub/Bioconductor version and
the reported snapshot. No package upgrades are performed. When the bundled
release is unknown, candidates cannot be described as proven newer versions.
Ribosomal protein GO membership still uses the separately recorded installed
OrgDb snapshot; updating an EnsDb does not update that package.

Default storage is `analysis/local/annotation-cache` in this development
worktree, already ignored by Git. It contains `checks/` (small check records),
`hub/` (AnnotationHub's catalogue/resource cache) and `references/` (versioned
export tables and manifests). No expression checkpoints are created.
There is no automatic deletion or eviction. Removing cached versions is a
separate maintenance decision because analyses may refer to them. A deployment
can supply `cache_dir` or the process-local `scrnaExplorer.annotationCache`
option to use an appropriate writable location.

```r
app_root <- normalizePath("site_data")
source(file.path(app_root, "R", "annotation_reference_cache.R"))

# Explicit network-capable maintenance; skips if a successful check is not due.
state <- check_annotation_updates("Homo sapiens", "GRCh37", app_root)
state$candidates

# Loading remains offline; use a reviewed reference ID from the app/cache.
source(file.path(app_root, "R", "gene_annotations.R"))
reference <- load_gene_annotation_cache("Homo sapiens", app_root,
    reference_id = "bundled", include_ribosomal = FALSE)
```

The standalone `load_cell_cycle_annotation()` accepts the same `reference_id`
and `cache_dir`, so a selected gene annotation can be reused for marker mapping.
This does not change the independently chosen S/G2M marker reference.

## Human genome builds and suffix matching

The locally cached AnnotationHub catalogue inspected on 2026-11-27 identifies
`AH119325` as the human GRCh38 Ensembl 113 EnsDb and `AH10684` as human GRCh37
Ensembl 75 GTF-derived GRanges. The latter is not an EnsDb; restricting the
catalogue to EnsDb alone would omit it. These are build-specific references,
not a claim to reproduce the exact annotation used to generate uploaded counts.
Match the original Ensembl release too when it is known.

The GRCh38 resource already present in the user's AnnotationHub cache was
exported offline into this worktree's versioned reference cache (about 2 MB).
Its snapshot is `2026-11-06`. The large SQLite database was reused, not copied
or downloaded again. With the user's explicit approval, `AH10684` was then
downloaded into the existing AnnotationHub cache and exported into the project
cache (about 0.6 MB; checked snapshot `2026-04-23`). Its GRanges sequence genome
explicitly records GRCh37 and agrees with the catalogue. These ignored local
exports are not included in Git; other deployments must populate their own
writable reference cache through the
explicit maintenance controls and download consent.

GTF export removes repeated exon/transcript annotation combinations, retaining
conflicting gene names for the existing ambiguity review. It does not aggregate
expression rows. GTF descriptions and Entrez IDs are unavailable and remain NA.
Consequently GO/Entrez-based ribosomal protein flags may be unresolved; RNA
biotypes remain available. The existing hemoglobin definition uses descriptions
starting with “hemoglobin”, so GTF-based hemoglobin metrics are unavailable
(NA, not zero). The annotation status warns about this, and existing QC filtering
continues to block on missing evidence. A richer matched reference or an explicit
future decision about an alternative hemoglobin definition is needed before
filtering with this GTF. A different source must not silently supply missing
annotations under the GTF release's provenance.

The suffix checkbox explains that a versioned ID such as
`ENSG00000141510.18` can be matched using `ENSG00000141510`. It changes only
lookup keys: source identifiers and counts remain unchanged, and no genome
build conversion occurs.

After reading counts, **Sample names** provides one text field per imported
sample and an **Apply sample names** button. Names must be nonempty and unique;
path separators and control characters are rejected. Draft text does not change
analysis state. Applying a changed name clears dependent preparation/QC through
the existing import invalidation path. Cell IDs, source filenames, metadata and
counts remain unchanged. The original sample label and rename history are
retained in the baseline SCE's import provenance.

## Verification

Validation on 2026-11-27: the full regression suite recorded 1,082 passing
expectations with zero failures, errors or skips. After adding the explicit
missing-hemoglobin message and its tests, the final focused annotation suite
passed 127 expectations, also with zero failures, errors or skips. Tests cover
both reference formats, build separation, conflict-preserving GTF export,
unavailable hemoglobin evidence, atomic sample renaming, stable cell IDs and
clearing stale annotation/QC/Seurat results through the actual app server.
Read-only preparation of the existing GSE176031 SCE loaded both real references,
preserved UMI totals and the 595-cell initial baseline membership, and confirmed
both reference choices in the Shiny control. The source SCE and bundled exports
were unchanged; no QC filtering decision or expression checkpoint was saved.

The reference tests inject catalogue and resource retrieval functions. They
exercise build/species eligibility, monthly boundaries, failures, download
consent, immutable versions, checksum failure, preservation of counts and QC
state, and explicit reference selection through the real Shiny QC controller.
They perform no online downloads and install no dependencies.

Validation on 2026-11-24: the full regression run recorded 562 passing
expectations, with zero failures, errors or skips. The focused reference tests
also passed after the final pre-check UI guard was added. Reusing the existing
GSE176031 pilot SCE, all five historical metrics matched by cell ID within
`1e-12`; the baseline remained 19,665 features by 595 cells. Sparse retained
counts were identical, and checksums of the pilot, historical metrics and both
bundled annotation files were unchanged. The original checkout's `site_data`
still matched `main`. Live online catalogue/resource retrieval was not exercised;
no annotation database was downloaded during development validation.

```r
testthat::test_dir("tests/testthat", filter = "annotation-reference-cache",
                   reporter = "summary", stop_on_failure = TRUE)
```

See the official [AnnotationHub documentation](https://bioconductor.org/packages/release/bioc/vignettes/AnnotationHub/inst/doc/AnnotationHub.html)
for catalogue metadata, snapshots and local caching.
