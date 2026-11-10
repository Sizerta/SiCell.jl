# Changelog

All notable changes to SiCell.jl are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `CODE_OF_CONDUCT.md` (Contributor Covenant 2.1), previously only linked
  from `CONTRIBUTING.md`.
- Regression tests for the differential-expression statistics (Scanpy-convention
  log2FC and per-cluster Benjamini-Hochberg over all genes) and for the
  orientation conventions of the `.h5ad` reader (obsm embeddings and dense X).
- Optional Muon.jl interoperability through `from_anndata` and `to_anndata`,
  including AnnData metadata, layers, embeddings, pairwise annotations and
  unstructured metadata preservation.
- Tests for the Trajectory Uncertainty Framework (`trajectory_uncertainty!`):
  TES/TDS/LPS columns are created, TES and TDS lie in [0, 1], and pseudotime
  outside [0, 1] raises an error.
- Synthetic multi-scale scaling benchmark (`benchmarks/synthetic_scaling_benchmark.jl`)
  covering 500–40,000 cells with replicate datasets, plus result tables.
- `CONTRIBUTING.md` with guidelines for bug reports, feature requests and pull requests.
- `CHANGELOG.md` (this file).

### Changed
- **`run_harmony!` is renamed to `run_centroid_alignment!`.** The
  implementation is a centroid-based alignment heuristic (global clustering
  plus batch-centroid shift), not the Harmony algorithm; the new name and
  docstring describe what the code actually does, and the docstring now states
  explicitly that results will differ from `scanpy.external.pp.harmony_integrate`.
  `run_harmony!` remains available as a deprecated alias that emits a warning
  and forwards to `run_centroid_alignment!`.
- `run_bbknn!` default `reduction` changed from `"harmony"` to `"centroid"`,
  matching the default storage key of `run_centroid_alignment!`. Pass
  `reduction` explicitly to use a custom embedding key.
- `find_all_markers` now computes `avg_log2FC` following the Scanpy
  convention, log2 of the ratio of linear-scale means
  (`log2((expm1(mean1) + 1e-9) / (expm1(mean2) + 1e-9))` on log1p-normalized
  means), instead of taking log2 of the log1p means directly.
- `find_all_markers` now applies Benjamini-Hochberg correction within each
  cluster over ALL genes before thresholding (matching Scanpy's
  `rank_genes_groups`), so `p_val_adj` no longer depends on `logfc_threshold`
  or `min_pct`; previously the adjustment pooled the already-filtered markers
  across clusters.
- `Project.toml`: declare Muon.jl as an optional weak dependency and register
  `SiCellMuonExt` as a Julia package extension.
- README: document Muon.jl-based AnnData import/export and count-layer selection.
- `CITATION.cff` updated to CFF 1.2.0 with full software metadata and a
  reference to the TUF preprint (doi:10.64898/2026.07.31.742012); removed a
  placeholder `preferred-citation` for a SiCell article that is not yet published.
- README: citation section now points to `CITATION.cff` and the TUF preprint;
  fixed the PanglaoDB spelling and the nesting of the feature lists.
- `Project.toml`: author field now lists the maintainer's full name.

### Fixed
- `run_harmony!` now retains its historical default reduction key (`"harmony"`) while emitting a deprecation warning; the new `run_centroid_alignment!` API defaults to `"centroid"`.
- `annotate_clusters!` validates the requested cluster metadata column before taking the empty-marker fallback, so a misspelled key is reported clearly.
- Muon interoperability now preserves nullable `DataFrame` annotations in `obsm`/`varm`; nullable numeric arrays in the SiCell passthrough payload are represented as DataFrames on export, matching Muon.jl AnnData’s documented constructor types.
- `read_h5ad` now raises an informative error for square dense matrices whose orientation cannot be inferred from dimensions, and accepts `dense_orientation=:obs_by_var` or `:var_by_obs` to resolve the ambiguity explicitly.
- CI now generates a small `.h5ad` fixture with Python `anndata` and checks that the native reader can import it.
- Static tutorials and README now describe centroid-based alignment accurately and distinguish cell-level marker discovery from replicate-aware condition-level inference.
- Removed an unused copied `reference/Project.toml.github-reference` file containing placeholder author metadata.
- `find_all_markers`: log2 fold change is now on the correct scale (see
  Changed above); the previous values were log2 of log1p means.
- `read_h5ad`: obsm embeddings are stored in the cells × dims orientation
  used by SiCell reductions. Previously the branch logic transposed embeddings
  whose observation axis is first (the AnnData convention), silently breaking
  `run_umap!` and `filter_cells!` for imported files.
- `read_h5ad`: dense (non-sparse) X matrices are now oriented against the
  obs/var metadata. Python's anndata stores dense X as obs × var while
  Muon.jl writes the transpose (genes × cells); both layouts are now read
  correctly. Previously the dense path applied no transpose and failed the
  dimension check for any non-square matrix.
- `read_h5ad`: obsm entries are converted to Float64 (h5ad commonly stores
  Float32) and 1-D entries are accepted as single-column embeddings; the
  `X_` prefix is now stripped only from the start of the key.
- `annotate_clusters!` now handles marker frames without usable genes (for
  example an empty `find_all_markers` result): it assigns 'Unknown' to all
  cells and warns instead of crashing with an ArgumentError.
- `find_all_markers` returns an empty DataFrame with the documented column
  schema (including `p_val_adj`) when no gene passes the thresholds.
- `run_bbknn!` builds its neighbor graph with `add_edge!` on a `SimpleGraph`
  instead of a dense cells × cells Bool adjacency matrix, removing the
  O(n²) memory allocation (~10 GB at 100k cells). The resulting graph is
  identical to the previous construction.
- `filter_cells!` now subsets cell embeddings along their observation axis
  (rows), matching the cells × dimensions layout used by AnnData and SiCell.

## [0.1.1] - 2026-06-18

First release registered in the Julia General registry.

### Added
- `SingleCellObject` data structure for counts, normalized data, metadata and reductions.
- Data loading for 10x Genomics matrices, HDF5 (`.h5`) and AnnData (`.h5ad`) files.
- Quality control: QC metrics, mitochondrial content, cell and gene filtering.
- Preprocessing: library-size normalization, log transformation, highly variable
  gene selection and feature scaling.
- Batch integration: centroid-based alignment heuristic, BBKNN, batch merging and a mixing score.
- Dimensionality reduction: randomized PCA, UMAP and diffusion maps.
- Clustering: KNN graph construction, graph-based clustering and k-means.
- Differential expression with a sparse-aware Wilcoxon rank-sum test and
  multiple-testing correction, and pseudobulk aggregation.
- Marker-based cell-type annotation using CellMarker and PanglaoDB references.
- Trajectory analysis: diffusion pseudotime and the Trajectory Uncertainty
  Framework (TES, TDS).
- Plotting: embeddings, feature plots, violin plots, heatmaps, volcano plots,
  PAGA-style connectivity graphs and trajectory uncertainty maps.
- Example pipelines (PBMC 3k, breast cancer, glioblastoma), documentation site,
  test suite, and CI/CompatHelper/TagBot workflows.

[Unreleased]: https://github.com/Sizerta/SiCell.jl/compare/v0.1.1...HEAD
[0.1.1]: https://github.com/Sizerta/SiCell.jl/releases/tag/v0.1.1
