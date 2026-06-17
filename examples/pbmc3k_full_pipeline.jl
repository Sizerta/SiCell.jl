# ==============================================================================
# SiCell — Full Pipeline Test
# Covers: QC, Clustering, Annotation, DE, Trajectory Inference,
#         Uncertainty-Aware Pseudotime (UAP), Branch Entropy
# Dataset: PBMC 3k (10x Genomics)
# ==============================================================================

using SiCell
using DataFrames
using Plots
using CSV
using Statistics
using Base.Threads
println("Julia threads available: ", Threads.nthreads())
# ── Helper ─────────────────────────────────────────────────────────────────────
function section(title::String)
    bar = "="^60
    println("\n$bar")
    println("  $title")
    println(bar)
end

function step(msg::String)
    println("\n▸ $msg")
end

# ==============================================================================
# 0. LOAD & PREPROCESS
# ==============================================================================
section("0. LOADING DATA & PREPROCESSING")

DATA_PATH = ".\\test_data\\pbmc3k_filtered_gene_bc_matrices\\filtered_gene_bc_matrices\\hg19"

obj = read_10x(DATA_PATH)
step("Loaded: $(size(obj.counts, 2)) cells × $(size(obj.counts, 1)) genes")

calculate_qc_metrics!(obj)
filter_cells!(obj, min_genes=200, max_mito=5.0)
step("After QC: $(size(obj.counts, 2)) cells retained")

normalize_data!(obj)
find_variable_features!(obj)

run_pca!(obj)
step("Preprocessing complete (normalize → HVG → PCA)")

# ==============================================================================
# 1. GRAPH, CLUSTERING & EMBEDDING
# ==============================================================================
section("1. CLUSTERING & EMBEDDING")

find_neighbors!(obj, k=30)
step("KNN graph built (k=30)")

step("Running K-Means clustering (k=8)...")
run_clustering!(obj, k=8)

step("Running Louvain clustering (try different resolutions)...")
run_graph_clustering!(obj, method="louvain", resolution=1.2)

step("Running UMAP...")
run_umap!(obj)

# ==============================================================================
# 2. ANNOTATION
# ==============================================================================
section("2. MARKER DETECTION & CELL TYPE ANNOTATION")

step("Finding markers for all K-Means clusters...")
markers = find_all_markers(obj, group="cluster")
println("  → Markers found for $(length(unique(markers.cluster))) clusters")

step("Annotating clusters using cell marker...")
db = load_cellmarker(species="Hs")

annotate_clusters!(obj, markers, db,
    species="Hs",
    min_score=0.03)

println("\n  Annotation summary:")
summary = combine(groupby(obj.meta_data, [:cluster, :cell_type]), nrow => :n_cells)
show(summary, allrows=true)

# ==============================================================================
# 3. STANDARD VISUALIZATION
# ==============================================================================
section("3. STANDARD VISUALIZATION")

step("UMAP: K-Means clusters")
savefig(dim_plot(obj, reduction="umap", group="cluster",
        title="UMAP — K-Means Clusters"), "umap_kmeans.png")

step("UMAP: Louvain clusters")
savefig(dim_plot(obj, reduction="umap", group="graph_cluster",
        title="UMAP — Louvain Clusters"), "umap_louvain.png")

step("UMAP: Cell type annotations")
savefig(dim_plot(obj, reduction="umap", group="cell_type",
        title="UMAP — Cell Type Annotations"), "umap_annotated.png")

step("Violin plot: CD79A & CD3D by cell type")
savefig(violin_plot(obj, features=["CD79A", "CD3D"], group="cell_type"),
    "violin_plot_annotated.png")

step("Pseudobulk heatmap")
pb_data = get_pseudobulk!(obj, "cell_type")
CSV.write("pseudobulk_counts.csv", pb_data)
savefig(plot_pseudobulk_heatmap(pb_data, ranking=:range, n_genes=40,
        title="Pseudobulk Signature Genes"), "pseudobulk_heatmap.png")

step("Volcano plot (cluster 1 markers)")
cluster_ids = unique(markers.cluster)
target_cluster = string(cluster_ids[1])
cluster_mrks = filter(row -> string(row.cluster) == target_cluster, markers)
savefig(volcano_plot(cluster_mrks, logfc_thresh=0.25, pval_thresh=0.05, top_n=10),
    "volcano_cluster_$(target_cluster).png")

step("PAGA connectivity graph (cell types)")
savefig(paga_plot(obj, reduction="umap", cluster_col="cell_type",
        min_conn_weight=0.02, title="Lineage Graph — Cell Types"),
    "paga_celltype_graph.png")

println("\n  ✓ Standard plots saved.")

# ==============================================================================
# 4. STANDARD TRAJECTORY INFERENCE
# ==============================================================================
section("4. STANDARD TRAJECTORY INFERENCE")

step("Building Diffusion Map (10 components)...")
run_diffusion_map!(obj, n_components=10)

# Auto-select a biologically meaningful root cell type
ROOT_CANDIDATES = [
    "Naive T cells", "CD4 Naive", "Naive CD4 T",
    "Naive CD8 T", "HSC", "Progenitor"
]
available_types = unique(obj.meta_data.cell_type)
root_type = something(
    findfirst(t -> t in available_types, ROOT_CANDIDATES),
    1   # fallback: first available type
) |> i -> (i isa Int ? available_types[i] : ROOT_CANDIDATES[i])

# Cleaner fallback: just pick first match or first available
root_type = let found = filter(t -> t in available_types, ROOT_CANDIDATES)
    isempty(found) ? available_types[1] : first(found)
end
step("Root cell type selected: '$root_type'")

step("Computing pseudotime from root...")
run_pseudotime!(obj, root_type, cluster_col="cell_type")

step("Plotting pseudotime on UMAP...")
savefig(feature_plot(obj, :pseudotime,
        title="Pseudotime from '$root_type'"), "umap_pseudotime.png")

# Gene dynamics over pseudotime (binned)
step("Generating pseudotime gene dynamics heatmap...")
obj.meta_data[!, :pseudotime_bin] = map(obj.meta_data.pseudotime) do t
    t < 0.33 ? "Early" : t < 0.66 ? "Mid" : "Late"
end
pb_time = get_pseudobulk!(obj, "pseudotime_bin")
savefig(plot_pseudobulk_heatmap(pb_time, n_genes=40,
        title="Gene Dynamics Over Pseudotime"), "pseudotime_heatmap.png")

println("\n  ✓ Standard trajectory plots saved.")

# ==============================================================================
# 5. TRAJECTORY UNCERTAINTY FRAMEWORK (TES + TDS + LPS)
# ==============================================================================
section("5. TRAJECTORY UNCERTAINTY FRAMEWORK")

step("Computing new trajectory uncertainty metrics (TES + TDS + LPS)...")
trajectory_uncertainty!(obj, key_prefix="traj_unc")

# Summary statistics
println("\n  Uncertainty summary:")
for metric in ["traj_unc_tes", "traj_unc_tds", "traj_unc_lps"]
    if metric in names(obj.meta_data)
        vals = obj.meta_data[!, metric]
        println("    $metric → μ=$(round(mean(vals), digits=3))   max=$(round(maximum(vals), digits=3))")
    end
end

println("\n  Mean metrics by cell type:")
summary_df = combine(
    groupby(obj.meta_data, :cell_type),
    :traj_unc_tes => mean => :mean_TES,
    :traj_unc_tds => mean => :mean_TDS,
    :traj_unc_lps => mean => :mean_LPS
)
sort!(summary_df, :mean_TES, rev=true)
show(summary_df, allrows=true)

step("Saving trajectory uncertainty plots...")
savefig(
    trajectory_uncertainty_plot(obj, title="Trajectory Uncertainty Framework"),
    "trajectory_uncertainty.png"
)

# Export full table
step("Exporting uncertainty metrics table...")
CSV.write("trajectory_uncertainty_metrics.csv",
    select(obj.meta_data, [:cell_type, :pseudotime, :traj_unc_tes, :traj_unc_tds, :traj_unc_lps]))
println("    ✓ trajectory_uncertainty_metrics.csv")
# ==============================================================================
# 6. SUMMARY
# ==============================================================================
section("TEST COMPLETE — OUTPUT FILES")

outputs = [
    ("umap_kmeans.png", "UMAP colored by K-Means clusters"),
    ("umap_annotated.png", "UMAP colored by cell type"),
    ("violin_plot_annotated.png", "Violin: CD79A & CD3D by cell type"),
    ("pseudobulk_heatmap.png", "Pseudobulk signature genes"),
    ("pseudobulk_counts.csv", "Raw pseudobulk counts table"),
    ("volcano_cluster_$(target_cluster).png", "Volcano plot (cluster $target_cluster)"),
    ("paga_kmeans_graph.png", "PAGA connectivity graph"),
    ("umap_pseudotime.png", "Standard pseudotime on UMAP"),
    ("pseudotime_heatmap.png", "Gene dynamics heatmap over pseudotime"),
    ("uap_uncertainty.png", "Trajectory confidence: color=pseudotime, opacity=confidence"),
    ("uap_cell_uncertainty.csv", "Per-cell pseudotime and trajectory confidence scores"),
]

for (fname, desc) in outputs
    marker = isfile(fname) ? "✓" : "·"
    println("  $marker  $fname  —  $desc")
end