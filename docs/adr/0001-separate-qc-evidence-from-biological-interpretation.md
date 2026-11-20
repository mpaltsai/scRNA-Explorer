# Separate QC evidence from biological interpretation

The new QC recommender uses technical QC evidence and experimental context to propose editable thresholds. Initially, cell-type annotations and stress/cell-cycle signatures inform diagnostic review rather than directly determining exclusions; doublet and ambient-RNA evidence remain distinguishable from damaged-cell evidence.

Using biological interpretations directly could make recommendations appear more specific, but could also remove genuine populations and create circular support for the downstream results. We accept a more conservative recommender that can report insufficient evidence, preserves the existing manual workflow, and requires independent evaluation rather than judging success by cleaner clustering or stronger differential expression. This records the accepted Q19–Q22 design, not a restriction on existing routines.
