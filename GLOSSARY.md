# scRNA-Explorer

Language for interactive single-cell RNA-seq quality control and analysis.

## Language

**Capture**:
A cell-capture experiment within which doublets can form. A capture can contain material from more than one donor when samples are multiplexed.
_Avoid_: Using sample, donor, and capture interchangeably

**Donor**:
The biological individual from which material originates. Multiple captures or specimens from one donor do not automatically constitute independent biological replicates.

**Biological sample**:
A specimen collected from a biological source under a defined experimental context. Samples from the same donor can represent different tissues, conditions, or time points and are not automatically independent replicates.

**Condition**:
The biological or experimental group of interest in a comparison, such as control or treatment.

**Raw droplet matrix**:
An unnormalized count matrix containing both called-cell barcodes and background droplets. This differs from called-cell counts that have not yet undergone user-selected QC.
_Avoid_: Unfiltered data without specifying whether background droplets are present

**Doublet prediction**:
A computational score or label indicating that a barcode may represent more than one captured cell. A prediction is evidence for review rather than a confirmed identity.

**QC threshold recommendation**:
A data-informed proposal for quality-control cutoffs that the user can inspect and edit before applying filtering. A recommendation does not itself exclude cells or establish their true quality.

**QC baseline**:
Called-cell data remaining after the initial requirement of at least 100 UMIs per cell. It provides the reference for a sample's subsequent QC filtering and is distinct from a raw droplet matrix.

**QC-filtered data**:
The cells retained from a sample's QC baseline by the user's applied QC filtering rules. Each sample is filtered separately before the retained data enter downstream analysis.

**Cluster annotation**:
A biological interpretation assigned to a cluster using marker evidence or reference data, with uncertainty retained when the evidence is insufficient.

**Annotation suggestion**:
A proposed biological label supported by an automated method or reference. It remains distinct from the label accepted or edited by the user.

**User-approved annotation**:
A biological label explicitly accepted or assigned by the user. It does not replace the original cluster identifier or erase automatic suggestions.

**Subclustering**:
An analysis of a selected cell population to investigate finer structure while retaining its relationship to the original population.

**Cluster marker analysis**:
An analysis identifying genes that distinguish a cluster from other cells in the dataset. It is distinct from testing a condition effect across biological replicates.

**Pseudobulk profile**:
Gene counts summed across cells belonging to a specified biological sample and cell population. Independent biological replicates, rather than individual cells, support condition comparisons.

**Analysis bundle**:
A downloadable record of an analysis containing its data object, experimental metadata, settings, software versions, filtering decisions, and annotation history. Parent and subcluster relationships remain part of that record.

**Step proposal**:
A candidate analysis result with its input reference, method settings, diagnostic evidence, and provenance, awaiting a user's decision. Producing a proposal does not approve it for downstream analysis.

**Analysis decision**:
An explicit acceptance or rejection of a step proposal, accompanied by a reason. Acceptance records the user's choice and supporting evidence, not a claim of statistical optimality.
