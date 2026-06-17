module SiCell

# ------------------------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------------------------
using SparseArrays
using DataFrames
using MatrixMarket
using CSV
using Clustering
using GZip
using Statistics
using HDF5
using Plots
using StatsPlots
using MultipleTesting
using Distributions
using StatsBase
using LinearAlgebra
using Graphs
using Arpack
using NearestNeighbors
using UMAP
using Base.Threads
using Random


# ------------------------------------------------------------------------------
# Module Includes
# ------------------------------------------------------------------------------
include("types.jl")
include("io.jl")                    # 10x, h5ad, and h5 readers
include("qc.jl")                    # Metric calculation and filtering
include("batch_normalization.jl")    # Harmony, BBKNN, and Batch alignment
include("analysis.jl")              # PCA, KNN, Clustering, UMAP
include("differential_expression.jl")
include("annotation.jl")            
include("trajectory_inference.jl")   # Trajectory + TUF (Core Feature)
include("plotting.jl")              # High-quality scRNA-seq visualizations

# ------------------------------------------------------------------------------
# Exported API
# ------------------------------------------------------------------------------

# Core & IO
export SingleCellObject, subset_object, read_10x, read_h5ad, read_h5

# Quality Control
export calculate_qc_metrics!, filter_cells!

# Preprocessing & Integration
export normalize_data!, find_variable_features!, scale_data!
export run_harmony!, run_bbknn!, merge_batches!, calculate_mixing_score

# Dimensionality Reduction & Clustering
export run_pca!, find_neighbors!, run_clustering!, run_graph_clustering!, run_umap!

# Differential Expression & Pseudobulk
export find_all_markers, get_pseudobulk!

# Trajectory Inference
export run_diffusion_map!, run_pseudotime!, trajectory_uncertainty!, trajectory_uncertainty_plot, uncertainty_plot, pseudotime_uncertainty!, highlight_top_uncertainty

# Annotation 
export load_cellmarker, load_panglaodb, annotate_clusters!, plot_annotation_heatmap

# Visualization
export dim_plot, volcano_plot, feature_plot, violin_plot, plot_pseudobulk_heatmap, paga_plot

end # module SiCell