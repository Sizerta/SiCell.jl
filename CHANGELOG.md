# Changelog

All notable changes to SiCell.jl are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Tests for the Trajectory Uncertainty Framework (`trajectory_uncertainty!`):
  TES/TDS/LPS columns are created, TES and TDS lie in [0, 1], and pseudotime
  outside [0, 1] raises an error.
- Synthetic multi-scale scaling benchmark (`benchmarks/synthetic_scaling_benchmark.jl`)
  covering 500–40,000 cells with replicate datasets, plus result tables.
- `CONTRIBUTING.md` with guidelines for bug reports, feature requests and pull requests.
- `CHANGELOG.md` (this file).

### Changed
- `CITATION.cff` updated to CFF 1.2.0 with full software metadata and a
  reference to the TUF preprint (doi:10.64898/2026.07.31.742012); removed a
  placeholder `preferred-citation` for a SiCell article that is not yet published.
- README: citation section now points to `CITATION.cff` and the TUF preprint;
  fixed the PanglaoDB spelling and the nesting of the feature lists.
- `Project.toml`: author field now lists the maintainer's full name.

## [0.1.1] - 2026-06-18

First release registered in the Julia General registry.

### Added
- `SingleCellObject` data structure for counts, normalized data, metadata and reductions.
- Data loading for 10x Genomics matrices, HDF5 (`.h5`) and AnnData (`.h5ad`) files.
- Quality control: QC metrics, mitochondrial content, cell and gene filtering.
- Preprocessing: library-size normalization, log transformation, highly variable
  gene selection and feature scaling.
- Batch integration: Harmony, BBKNN, batch merging and a mixing score.
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
