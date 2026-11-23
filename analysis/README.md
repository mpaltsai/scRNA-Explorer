# Standalone R feature development

Develop and validate R features against local datasets before integrating them into Shiny. Existing app scripts remain unchanged by this input/metadata work. Whether the app is extended or rebuilt is a later decision.

## Register GSE176031

The source is the existing compressed gene-by-cell DGE collection plus its GEO series matrix. This workflow does not use AnnData. The registry contains paths, barcode headers, sample metadata, provenance, and diagnostics; it is not a replacement for SingleCellExperiment or Seurat and contains no expression array.

From the repository root, configure the existing input directory and run:

```bash
GSE176031_DATA_ROOT=/path/to/GSE176031
Rscript --vanilla scripts/register_gse176031.R \
  --data-root "$GSE176031_DATA_ROOT" \
  --output-dir analysis/local/GSE176031 \
  --confirmed-genome-build hg19
```

The default metadata source is the official [GEO series matrix](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176031/matrix/GSE176031_series_matrix.txt.gz). For offline use or if NCBI denies a request, add `--series-matrix /path/to/GSE176031_series_matrix.txt.gz`. Paths are runtime inputs, not committed personal paths.

Outputs are local and ignored by git:

- `registry.rds`: reusable metadata/header records for subsequent tickets.
- `samples.csv`: accession-based join, original patient/tissue/specimen fields, declared and confirmed genome builds, and available DGE inputs.
- `source_metadata.csv`: every original sample field/value, including repeated fields, with source row and source field.
- `inputs.csv`: file inventory and header-derived cell counts. These are not validated expression-count values.
- `issues.csv`: missing/unmatched/duplicate inputs, invalid headers, source conflicts, and unresolved relationships.
- `artifacts.csv`: optional explicitly configured R artifact references; empty when none are supplied. R serialization headers alone do not establish class, dataset association, or count compatibility.

Nonempty output directories are never overwritten. Reuse a saved registry with `readRDS()`; choose a new output directory for a new snapshot. A failed input check produces diagnostics and a nonzero command status rather than silently selecting an ambiguous input.

The user confirmed hg19. `confirmed_genome_build` records that confirmation; the contradictory original GEO `mm8` text remains in source provenance and is marked resolved by explicit confirmation. This does not alter gene IDs or expression values.

GEO accession, donor, biological sample, capture, and sequencing pool are different concepts. Patient metadata is retained as donor evidence, but capture/pool and biological-sample mappings remain unresolved until supported. Repeated barcode strings across files require source context, not automatic merging. The archive name `RAW` is not evidence of background droplets or independent cell-quality labels.

## Reuse inputs and existing analysis objects

```r
source("analysis/R/gse176031_registry.R")
registry <- readRDS("analysis/local/GSE176031/registry.rds")
input <- resolve_gse176031_input(registry, "GSM5353214")
input$dge_path
```

Resolution checks that the registered file exists and that its size and barcode header still match. These are lightweight checks, not a full content checksum or count-value validation. After relocating the same files, pass `data_root` explicitly to the resolver. It does not load counts or create an analysis object.

Attach the metadata to an **existing** SingleCellExperiment or Seurat object:

```r
sce <- attach_geo_sample_metadata(sce, registry, "GSM5353214")
seu <- attach_geo_sample_metadata(seu, registry, "GSM5353214")
```

For multiple samples, supply one accession per cell in object order, or a character vector named by the object's cell identifiers. Named mappings must match the cells exactly. This supports existing objects from different input formats; it does not assume that a barcode itself identifies a sample.

All source sample fields are projected into namespaced cell metadata; standardized fields include `geo_accession`, `geo_donor_id`, and `geo_confirmed_genome_build`. Original source tables, source-field/column mappings, series metadata, and provenance are retained in `S4Vectors::metadata(sce)$geo` or `seu@misc$geo`. Reserved standard columns remain distinct from similarly named source characteristics. Existing fields/provenance with different values cause an error instead of being overwritten.

Counts, other assays, cell order, existing sample/QC fields, and other object metadata stay unchanged. Filtering/subsetting retains the attached fields through normal object behavior. When converting representations, explicitly retain or reattach the GEO record using the current cell mapping; the existing Shiny normalization code has not yet been changed to do this automatically. The helper does not normalize, filter, cluster, or create a new expression-data container.

Optionally pass `--r-artifact /path/to/existing.rds` to register an existing R file without deserializing it. No AnnData files are inspected or converted, and no Python runtime is required.

## Checks

```bash
Rscript --vanilla tests/run_registry_tests.R
```

The registry command uses base R. Tests use an existing `testthat` installation; object-attachment tests additionally use installed Matrix, SingleCellExperiment, and SeuratObject packages. Nothing installs dependencies automatically.

The checks cover accession joins, missing/duplicate/unmatched inputs, reference provenance, changed headers/relocation, source-field collisions, SCE/Seurat attachment with unchanged counts, and command output reuse. Real-data checks can resolve every registered accession without loading expression arrays. This establishes an input/metadata foundation, not biological validation of a new QC or annotation method.
