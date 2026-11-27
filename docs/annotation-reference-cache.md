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
   updates**. This action checks the AnnotationHub catalogue online. It may download AnnotationHub's
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
continues to block on missing evidence for the unsupplemented reference. The
separate BioMart supplement below can restore the available descriptions and
Entrez links without changing the existing gene definitions. Its provenance
explicitly distinguishes the two sources.

## Supplementing the GRCh37 reference

With human counts and **GRCh37 / hg19** selected, expand **Annotation reference
maintenance**. **Complete missing GRCh37 annotation fields** accepts an
unsupplemented, cached human Ensembl 75 GTF reference. Give a reason, confirm
retrieval and click **Add missing fields from GRCh37 BioMart**. This is a second
explicit network action; ordinary reference loading and QC preparation still
stay offline. Errors or unavailable services leave existing references and
prepared QC active. No packages are installed and no global settings change.

The implementation verifies the dedicated service host and the returned human
dataset's GRCh37 assembly. Ensembl documents its GRCh37 gene models as frozen at
release 75, but the BioMart service may report a newer release. On 2026-11-27 it
reported `Ensembl Genes 116` with dataset `GRCh37.p13`. These are recorded
separately: this does not describe current descriptions/cross-references as a
historical 2014 snapshot. See [Ensembl's GRCh37 documentation](https://jun2026.archive.ensembl.org/info/website/tutorials/grch37.html).

Queries use BioMart's XML interface from R with the existing `httr2` and `xml2`
packages. Redirects are disabled. Registry/dataset/attribute metadata are
validated before querying, and a `[success]` completion stamp is required before
accepting the data. The source metadata, actual mart and dataset versions,
retrieval time, requested fields, query XML, metadata-response hashes and the
normalized gene-row hash are retained.
This avoids relying on a default BioMart assembly or an interrupted response.
The protocol is documented by [Ensembl](https://useast.ensembl.org/info/data/biomart/biomart_restful.html).

Join policy:

- Match by exact Ensembl gene ID. Preserve the parent annotation rows, order,
  gene names, biotypes and chromosomes; never join expression matrices.
- Fill only empty descriptions/Entrez links. Preserve existing values and
  record disagreements for review.
- Fill a missing description only when the returned nonempty descriptions
  agree. Conflicting descriptions remain NA.
- Retain all Entrez candidates as a list per gene. Ribosomal membership uses
  the existing agreement rule across candidates; it never takes the first ID.
- Keep unavailable fields unresolved. Extra BioMart genes are not added to
  the parent reference or counts.

A separate versioned reference contains `reference.rds`, a compact returned
`supplement.rds`, `supplement-review.csv` and its manifest. Parent provenance,
checksums, coverage and merge policy are recorded. Repeated identical payloads
reuse the same version; retrieval time alone does not create a new reference.
The review records original values, candidates and each per-row outcome; the
UI shows coverage and provides **Download BioMart supplement review** when a
supplemented reference is selected. Checksums are verified before export.

Creating a supplement does not select it or change prepared results. Select
**GRCh37 — Ensembl 75 (AH10684, GTF + BioMart)** under **Annotation reference**,
then **Add gene annotations**, review the mapping and rebuild the **UMIs >= 100**
baseline. This explicit selection clears dependent results through the existing
invalidation path. No count matrix or expression checkpoint is saved by the
maintenance operation. The primary resource, supplementary provider and
separate installed OrgDb GO source remain visible in provenance.

## Automatic identifier suggestions and lookup-only suffix matching

Each sample defaults to **Automatic detection** after reading data. Ensembl
gene-ID format suggests gene-ID lookup, including IDs with a numeric version
suffix. Otherwise exact matches in the selected local reference support a
gene-symbol or reference-ID suggestion. Unknown rows remain unknown; there is
no rule that every identifier without an ENS prefix must be a gene symbol.
Mixed human/mouse prefixes, mixed gene IDs/symbols, numeric IDs and
transcript/protein/exon IDs require an explicit lookup choice or corrected input.
Species and counts genome build remain explicit user declarations.

The identifier review table shows the suggested and selected type, automatic
versus manual mode, feature/Ensembl/versioned counts and exact symbol matches.
Exact symbol-match counts include names with multiple candidate genes; they
indicate presence in the reference, not a unique gene mapping.
A per-sample manual choice remains available. Browser choices are scoped to
the current source-feature list, so a stale selection from a different upload
cannot silently override the new suggestion. App-produced SCEs reuse their
recorded source IDs rather than the internal names used to retain duplicate rows.
The effective choice and its supporting evidence are recorded in the SCE QC
configuration as `identifier_decision`.

Version-suffix matching is enabled by default when versioned Ensembl gene IDs
are detected in the imported samples. The checkbox remains visible and can be
turned off. It only affects gene-ID lookup; symbol lookup remains exact. The
suffix checkbox explains that a versioned ID such as
`ENSG00000141510.18` can be matched using `ENSG00000141510`. It changes only
lookup keys: source identifiers and counts remain unchanged, and no genome
build conversion occurs. Initial browser bindings that acknowledge the same
automatic defaults do not invalidate prepared QC. Actual import/reference/
identifier/suffix changes still clear dependent results for preparation/review.
Unrelated app inputs do not change the identifier configuration.

## Unmapped, ambiguous and ribosomal feature review

The complete **Download gene mapping** table retains every input feature row,
source ID, mapped fields, candidate gene IDs and each group's TRUE/FALSE/NA
classification. `unmapped` means no reference row matched the lookup key;
`partial` means matching rows lack a usable gene ID or name; `ambiguous` means
more than one gene ID or name matched. An ambiguous gene mapping can still
support a particular group if all candidates agree on that group's membership.

Users can check identifier type, species/build, suffix matching and reference,
select the GRCh37 supplement where needed, download the full mapping and verify
IDs against the original provider's feature table. A verified source correction
can be re-uploaded. Gene symbols alone cannot always resolve a unique Ensembl
gene; the app does not offer a per-feature override editor in this version.
Keeping unresolved rows is supported: they remain in counts while annotation-
dependent information may be incomplete. No first-candidate selection, gene
merging or feature removal occurs automatically.

The ribosomal summary counts **matrix feature rows**, including rows with zero
counts, not cells or UMIs. For each group, `annotated_features` counts TRUE flags
(supported group membership) and `unresolved_features` counts NA flags (unknown
membership). FALSE flags are not annotated as members in the selected reference;
they are not proof of biological nonmembership. Missing mappings, unavailable
Entrez links/GO resources and disagreement among candidates can produce NA.

| Group | Definition used by the app |
| --- | --- |
| `ribosomal_protein` | Protein-coding structural constituents of ribosome, supported by GO:0003735. |
| `cytosolic_ribosomal_protein` | The supported structural proteins also assigned to cytosolic ribosome (GO:0022626). |
| `mitochondrial_ribosomal_protein` | The supported structural proteins also assigned to mitochondrial ribosome (GO:0005761). |
| `nuclear_rRNA` | rRNA-biotype features outside chromosome MT. |
| `mitochondrial_rRNA` | Mt_rRNA biotype, or rRNA biotype on chromosome MT. |
| `rRNA_pseudogene` | A biotype containing rRNA and pseudogene. |

Parent/subgroup protein flags overlap and their counts must not be added.
Protein and RNA groups have different evidence requirements, so their unresolved
totals can differ. Zero supported features does not prove absence while some
features remain unresolved. Users can filter the downloaded mapping by TRUE or
NA in the relevant group column to inspect the exact features counted.

After reading counts, **Sample names** provides one text field per imported
sample and an **Apply sample names** button. Names must be nonempty and unique;
path separators and control characters are rejected. Draft text does not change
analysis state. Applying a changed name clears dependent preparation/QC through
the existing import invalidation path. Cell IDs, source filenames, metadata and
counts remain unchanged. The original sample label and rename history are
retained in the baseline SCE's import provenance.

## Verification

Automatic identifier validation on 2026-11-27: the full regression suite passed
1,173 expectations with zero failures, errors or skips. Final focused identifier
and annotation-controller tests passed 69 expectations, including an explicit
mixed-identifier review/override case and the final review table. Tests verify
automatic symbol/Ensembl recognition, suffix defaults and opt-out, preserved
versioned count row names, stale-choice protection, unchanged preparation on
initial browser bindings/unrelated inputs and correct feature-count semantics
even for zero-count rows.

The existing GSE176031 pilot was correctly detected as gene symbols without a
manual identifier selection. Counts and the 595-cell baseline were unchanged,
and hemoglobin metrics remained available. Its ribosomal summary matched the
previous result. All 73 unresolved rows in each RNA group were unmapped; the
359 unresolved ribosomal-protein rows comprised 73 unmapped, 14 ambiguous and
272 mapped rows with unresolved protein classification. No expression checkpoint
or QC filtering decision was saved; the source SCE checksum was unchanged.

Validation of BioMart supplementation on 2026-11-27: the live XML query returned
66,978 annotation rows with a completion stamp. The 63,677 parent reference rows
were preserved; 40,611 missing descriptions and 27,811 sets of Entrez links were
filled, including 1,751 genes with multiple Entrez candidates. Missing values
remain missing and are reported in the downloadable review. The cached
supplemented annotation table is about 1.3 MB; all cache artifacts are ignored by
Git, with no expression checkpoints.

The existing GSE176031 SCE was reused through the real Shiny annotation and QC
controllers. Hemoglobin evidence became available (five matching features),
ribosomal classifications used the recorded OrgDb source and restored Entrez
links, and the initial baseline remained 595 cells with identical retained
counts. Default QC Preview and overview, joint-metric and ribosomal plot outputs
passed. No QC Apply decision was made. Checksums of the source SCE, bundled
exports and original GRCh37 reference were unchanged. The full regression suite
passed 1,134 expectations with zero failures, errors or skips. Focused supplement
tests also cover service failures without invalidating active QC, partial query
rejection, build checks, multiple mappings, preserved conflicts and audit-file
checksum validation.

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
