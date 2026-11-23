# Validate the R analysis core before replacing Shiny logic

Develop and test analysis features as standalone R functions operating on the existing SingleCellExperiment/Seurat objects. Record method settings, provenance, diagnostic evidence, and explicit user decisions through a shared runner. Shiny integration follows validation across datasets; whether to replace the current app or integrate incrementally remains a later choice.

This separates scientific decisions from reactive UI behavior so the tested implementation can become the future app's analysis core. It costs a small runner/log contract now, but avoids implementing the same algorithm separately for experiments and Shiny. Existing expression objects are reused; decision records do not contain automatic full-object checkpoints. Accepted settings are recorded as accepted, not presumed statistically optimal.
