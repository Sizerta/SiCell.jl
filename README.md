# SiCell.jl — High-Performance Single-Cell RNA-seq Analysis in Native Julia
![SiCell Poster](docs/Poster_SiCell.png)

[![CI](https://github.com/Sizerta/SiCell.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/Sizerta/SiCell.jl/actions/workflows/CI.yml)
[![Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://sizerta.github.io/SiCell.jl)
[![SiCell](https://img.shields.io/badge/registry-General-blue)](https://juliahub.com/ui/Packages/General/SiCell)
[![Version](https://juliahub.com/docs/SiCell/version.svg)](https://juliahub.com/ui/Packages/General/SiCell)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)



SiCell.jl is a modern, Julia-native toolkit for single-cell RNA sequencing (scRNA-seq) analysis. It provides a complete end-to-end workflow from raw count matrices to biological interpretation, combining high performance, scalable sparse computations, and an intuitive API.

---



# Installation

```julia
using Pkg
Pkg.add("SiCell")

```

---

# Quick Start

```julia
using SiCell

# Load data
obj = read_10x("path/to/filtered_feature_bc_matrix")

# Quality control
calculate_qc_metrics!(obj)
filter_cells!(
    obj;
    min_genes=200,
    max_mito=5.0
)

# Preprocessing
normalize_data!(obj)
find_variable_features!(obj)

# Dimensionality reduction
run_pca!(obj)

# Graph construction and clustering
find_neighbors!(obj, k=20)
run_graph_clustering!(obj)

# Visualization
run_umap!(obj)

dim_plot(
    obj,
    reduction="umap",
    group="graph_cluster"
)
```

---


#  Documentation
Please look at the site below:

https://sizerta.github.io/SiCell.jl/

Documentation includes:

* Getting started tutorials
* Complete analysis workflows
* Example Applications
* Case studies
* Visualization examples
* API reference
---
#  Example Applications

SiCell has been tested on multiple biological systems, including:

* Human PBMC datasets
* Large-scale tumor microenvironment datasets
  * Breast cancer (>34,000 cells)
  * Glioblastoma datasets

Example analyses include:

* Automated cell-type annotation
* Marker discovery
* Differential expression
* Trajectory reconstruction
* Identification of uncertain transition states using TUF


---
## Related Software

Interested in quantifying uncertainty in trajectory inference?

**TUF (Trajectory Uncertainty Framework)** is a standalone Python package that introduces two complementary uncertainty metrics for single-cell trajectory analysis:

- **Temporal Entropy Score (TES)** — measures local temporal mixing.
- **Trajectory Divergence Score (TDS)** — quantifies directional ambiguity near branching regions.

TUF is designed for the Scanpy ecosystem and can be used alongside any pseudotime inference method that stores results in an `AnnData` object.

**Repository:** https://github.com/Sizerta/tuf_python

---
# 🚀 Key Features

##  Complete Single-Cell Analysis Pipeline

- Data loading:
  - 10x Genomics matrices
  - Hierarchical Data Format(`.h5`)
  - AnnData (`.h5ad`) files

- Quality control:
  - Cell filtering
  - Mitochondrial content analysis
  - Gene and count statistics

- Preprocessing:
  - Library-size normalization
  - Log transformation
  - Highly variable gene selection
  - Feature scaling

- Dimensionality reduction:
  - Randomized PCA
  - UMAP
  - Diffusion Maps

- Clustering and batch integration:
  - Graph-based clustering
  - K-Means clustering
  - Efficient KNN graph construction
  - BBKNN and centroid-based batch alignment (a lightweight heuristic, not the Harmony algorithm)

- Differential expression:
  - Sparse-aware Wilcoxon rank-sum marker testing
  - Multiple-testing correction
  - Pseudobulk count aggregation by metadata group

`find_all_markers` is a cell-level, cluster-versus-rest marker-discovery method and is best treated as exploratory. Cells from the same biological sample are not independent biological replicates. For condition-level differential inference with replicated samples, aggregate raw counts by a sample/donor identifier using `get_pseudobulk!` and fit an appropriate replicate-aware model to the sample-level counts. `get_pseudobulk!` performs aggregation only; it does not fit that model or calculate condition-level differential-expression statistics.

- Cell type annotation:
  - Marker-based annotation using CellMarker and PanglaoDB as references

- Trajectory analysis:
  - Diffusion pseudotime
  - **Trajectory Uncertainty Framework (TUF)**:
    - TES (Temporal Entropy Score) — local temporal heterogeneity
    - TDS (Trajectory Divergence Score) — directional uncertainty and branching

- Publication-quality visualization:
  - UMAP embeddings
  - Feature plots
  - Violin plots
  - Heatmaps
  - Volcano plots
  - PAGA-style connectivity graphs
  - Trajectory uncertainty maps

---

---

# AnnData interoperability (optional)

SiCell provides an optional interoperability layer for [Muon.jl](https://github.com/scverse/Muon.jl), which supplies Julia-side AnnData and MuData data structures. This makes it possible to move an in-memory AnnData object into SiCell for analysis and export the results back to AnnData. Muon.jl is an optional dependency; the core SiCell workflow does not require it.

Install Muon.jl in the environment where you use SiCell:

```julia
using Pkg
Pkg.add("Muon")
```

Load both packages to activate the interoperability extension:

```julia
using SiCell
using Muon
```

## Import an AnnData object

```julia
# Read an AnnData file with Muon.jl
adata = Muon.readh5ad("pbmc.h5ad")

# Import adata.X into SiCell
obj = SiCell.from_anndata(adata)

# Alternatively, use a raw-count layer when adata.X is normalized
# obj = SiCell.from_anndata(adata; layer="counts")
```

## Analyze and export

```julia
# Run a SiCell workflow
calculate_qc_metrics!(obj)
normalize_data!(obj)
run_pca!(obj)
run_umap!(obj)

# Export the normalized matrix and analysis results as AnnData
adata_out = SiCell.to_anndata(obj; X=:normalized)
Muon.writeh5ad("pbmc_sicell.h5ad", adata_out)
```

## Data mapping and current limitations

AnnData stores expression matrices as **cells × genes**, whereas SiCell stores them as **genes × cells**. `SiCell.from_anndata` transposes the selected matrix during import; `SiCell.to_anndata` restores the AnnData orientation during export.

The adapter carries the observation and variable tables and their names, numeric matrix layers, embeddings, pairwise annotations (`obsp` and `varp`), and unstructured metadata (`uns`). DataFrame-valued annotations, including columns with `missing` values, are preserved. If a nullable numeric array is present in SiCell's passthrough payload, it is represented as a DataFrame on export because Muon.jl's documented AnnData constructor accepts numeric arrays or DataFrames. SiCell reductions are written to `obsm`; SiCell neighbor connectivity and distance information are exported to `obsp["connectivities"]` and `obsp["distances"]`, respectively.

The selected input matrix (`adata.X` by default, or the layer specified with `layer`) becomes SiCell's count matrix. If `adata.X` already contains normalized expression, explicitly select the raw-count layer when importing, for example `layer="counts"`.

The native `read_h5ad` reader infers dense matrix orientation from the observation and variable dimensions. For square dense matrices, where shape alone is ambiguous, pass `dense_orientation=:obs_by_var` or `dense_orientation=:var_by_obs` explicitly.

`SiCell.to_anndata` supports `X=:counts` (the default), `X=:normalized`, and `X=:scaled`. When `X` is `:normalized` or `:scaled`, the original counts are also written to `layers["counts"]`. Exporting `X=:scaled` requires scaled data containing all genes; scaled data containing only a subset of highly variable genes cannot be exported as `X`. Backed AnnData objects are not currently supported; read the AnnData object into memory before calling `SiCell.from_anndata`.

---

# ⚡ Performance

SiCell is designed around efficient sparse matrix operations, multithreading, and memory-aware algorithms.

### Performance highlights

* Native Julia implementation
* Sparse-aware algorithms
* Multi-threaded computations
* Randomized SVD-based PCA for large datasets
* Efficient graph construction and downstream analysis

Benchmarks are currently being expanded on datasets ranging from PBMCs to large-scale tumor atlases.

---
#  Citation

If you use SiCell.jl in your research, please cite the software. Citation metadata is provided in [`CITATION.cff`](CITATION.cff) (GitHub's "Cite this repository" button generates BibTeX/APA from it):

> Mahdavifar M. *SiCell.jl: A High-Performance, Native Julia Framework for End-to-End Single-Cell RNA-Sequencing Analysis* (Version 0.1.1) [Computer software]. https://github.com/Sizerta/SiCell.jl

If you use the Trajectory Uncertainty Framework (TES/TDS) in SiCell, please also cite the TUF preprint:

> Mahdavifar M, Mohammadifar Z, Iranpourtari T. Trajectory Uncertainty Framework (TUF): A Modular Framework for Identifying Transitional and Branch-Point Cell States in Single-Cell Trajectory Analysis. *bioRxiv* (2026). https://doi.org/10.64898/2026.07.31.742012

---

#  Contributing

Bug reports, feature requests, and pull requests are welcome. See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

---

# License

SiCell.jl is released under the MIT License.
