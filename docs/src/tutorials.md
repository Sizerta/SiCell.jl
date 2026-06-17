# Tutorials

This section provides step-by-step walkthroughs for common single-cell RNA sequencing (scRNA-seq) analysis workflows using SiCell.jl.

---

# 1. Basic Analysis Workflow

This tutorial demonstrates the standard analysis pipeline from raw count matrices to clustered cell populations and biological interpretation.

## 📂 Loading Data

SiCell.jl supports commonly used single-cell data formats.

```julia
using SiCell

# Load 10x Genomics output
obj = read_10x("data/filtered_feature_bc_matrix/")

# Load an AnnData (.h5ad) file
# obj = read_h5ad("data/dataset.h5ad")
```

## 🔬 Quality Control

Calculate quality-control metrics and remove low-quality cells.

```julia
calculate_qc_metrics!(obj)

filter_cells!(
    obj;
    min_genes=200,
    max_genes=5000,
    max_mito=10.0
)
```

## ⚖️ Normalization and Feature Selection

Normalize expression values and identify highly variable genes.

```julia
# Library-size normalization followed by log transformation
normalize_data!(obj)

# Identify highly variable genes
find_variable_features!(
    obj;
    n_features=2000
)

# Optional: scale gene expression values
scale_data!(obj)
```

## 📊 Dimensionality Reduction

Reduce the dimensionality of the dataset and construct the neighborhood graph.

```julia
# Principal component analysis
run_pca!(obj)

# Construct K-nearest-neighbor graph
find_neighbors!(
    obj;
    k=20,
    dims=30
)

# Compute UMAP embedding
run_umap!(obj)
```

## 🧩  Clustering

Identify cellular populations using graph-based community detection.

```julia
run_clustering!(obj, k=12)
run_graph_clustering!(obj, method="label_propagation", key="graph_cluster")
run_graph_clustering!(obj, method="louvain", key="graph_cluster")

```

## 🎨 Visualization

Visualize cell clusters and gene expression patterns.

```julia
# Visualize clusters
dim_plot(
    obj,
    reduction="umap",
    group="graph_cluster"
)

# Visualize expression of a marker gene
feature_plot(
    obj,
    "CD3D"
)
```

---

# 2. 🔗 Batch Integration with Harmony

For datasets containing multiple samples, patients, or experimental batches, SiCell.jl provides Harmony-based batch correction and integration.

```julia
# Merge datasets
merged_obj = merge_batches!(
    obj1,
    obj2;
    batch1_name="Patient1",
    batch2_name="Patient2"
)

# Standard preprocessing
calculate_qc_metrics!(merged_obj)
normalize_data!(merged_obj)
find_variable_features!(merged_obj)
scale_data!(merged_obj)
run_pca!(merged_obj)

# Harmony correction
run_harmony!(
    merged_obj,
    batch_key="batch"
)

# Construct graph using Harmony embedding
find_neighbors!(
    merged_obj,
    reduction="harmony"
)

# Visualization and  etc
run_umap!(
    merged_obj,
    reduction="harmony"
)


```

---

# 3. 🌊 Trajectory Inference

SiCell.jl provides diffusion-based pseudotime analysis for studying continuous biological processes such as differentiation, cellular activation, and disease progression.

## Diffusion Maps

Compute a low-dimensional manifold representation suitable for trajectory analysis.

```julia
run_diffusion_map!(
    obj;
    n_components=10
)
```

## ⏱️ Pseudotime

Select a biologically meaningful root cell or population and compute pseudotime.

```julia
# Example using a root cell index
root_cell = 1

run_pseudotime!(
    obj,
    root_cell;
    method="graph"
)
```

Pseudotime values are stored in:

```julia
obj.meta_data.pseudotime
```

## 🎨 Visualize Pseudotime

```julia
dim_plot(
    obj,
    reduction="umap",
    group="pseudotime"
)
```

---

# 4. 🌱 Trajectory Uncertainty Framework (TUF)

Traditional pseudotime methods assign cells a position along a trajectory but do not quantify the local uncertainty associated with that trajectory.

The Trajectory Uncertainty Framework (TUF) introduces two complementary metrics:

* **Temporal Entropy Score (TES):** quantifies temporal inconsistency among neighboring cells.
* **Trajectory Divergence Score (TDS):** quantifies disagreement among forward developmental directions.

Together, TES and TDS provide complementary views of local trajectory uncertainty.

## 🧮 Compute TES and TDS

After constructing a neighborhood graph and computing pseudotime:

```julia
trajectory_uncertainty!(
    obj;
    pseudotime_key="pseudotime",
    reduction="diffusion"
)
```

The computed scores are stored in:

```julia
obj.meta_data.traj_unc_tes
obj.meta_data.traj_unc_tds
```

## 🎨 Visualize TES and TDS

```julia
    trajectory_uncertainty_plot(obj; reduction="umap", key_prefix="traj_unc")
```

## Interpreting TUF

| TES  | TDS  | Interpretation                                                        |
| ---- | ---- | --------------------------------------------------------------------- |
| Low  | Low  | Stable and well-defined cellular trajectory                           |
| High | Low  | Temporal mixing or heterogeneous developmental states                 |
| Low  | High | Directional divergence and potential branching                        |
| High | High | Complex transitions involving both temporal and directional ambiguity |

---

# 5. 🧬 Differential Expression Analysis

Identify genes associated with specific cellular populations or states.

```julia
find_markers!(
    obj;
    group="graph_cluster"
)
```

The output includes:

* Differentially expressed genes
* Effect sizes
* Statistical significance values
* Multiple-testing corrected p-values

---

# 6. 🏷️ Automated Cell Type Annotation

SiCell.jl supports marker-based cell type annotation using Cell marker 2.

```julia
db = load_cellmarker(species="Hs")

annotate_clusters!(obj, markers, db,
    species="Hs",
    min_score=0.03)
```

---

# 🚀 What's Next?

After completing the basic workflow, you can:

* Perform differential expression analysis to identify marker genes.
* Annotate cellular populations using PanglaoDB.
* Integrate multi-sample datasets with Harmony.
* Study developmental and disease trajectories using pseudotime.
* Identify transitional and branching states using TES and TDS.
* Generate high quality visualizations.

For more advanced examples, consult the full API documentation and case-study notebooks.
