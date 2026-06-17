# ==============================================================================
# SiCell — Breast Cancer Pipeline Test
# Covers: QC, Clustering, Annotation, DE, Trajectory Inference,
#         Trajectory Uncertainty Framework (TES + TDS + LPS)
# Dataset: 320k scFFPE Breast Cancer (10x GEM-X FLEX, 16-plex)
# ==============================================================================

using SiCell
using DataFrames
using Plots
using CSV
using Statistics

# ── Helpers ────────────────────────────────────────────────────────────────────
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

DATA_PATH = ".\\test_data\\320k_scFFPE_16-plex_GEM-X_FLEX_BreastCancer1_BC7-8_count_sample_raw_feature_bc_matrix.h5"

obj = read_h5(DATA_PATH)
step("Loaded: $(size(obj.counts, 2)) cells × $(size(obj.counts, 1)) genes")

calculate_qc_metrics!(obj)
# FFPE data typically has lower RNA quality — relax mito threshold to 10%
filter_cells!(obj, min_genes=200, max_mito=10.0)
step("After QC: $(size(obj.counts, 2)) cells retained")

normalize_data!(obj)
find_variable_features!(obj)
run_pca!(obj)
step("Preprocessing complete (normalize → HVG → PCA)")

# ==============================================================================
# 1. GRAPH, CLUSTERING & EMBEDDING
# ==============================================================================
section("1. CLUSTERING & EMBEDDING")

find_neighbors!(obj, k=20, dims=30)
step("KNN graph built (k=20, dims=30)")

step("Running K-Means clustering (k=8)...")
run_clustering!(obj, k=8)

step("Running Label Propagation clustering...")
run_graph_clustering!(obj, method="label_propagation", key="graph_cluster")

step("Running UMAP...")
run_umap!(obj)

# ==============================================================================
# 2. ANNOTATION
# ==============================================================================
section("2. MARKER DETECTION & CELL TYPE ANNOTATION")

step("Finding markers for all K-Means clusters...")
markers = find_all_markers(obj, group="cluster")
println("  → Markers found for $(length(unique(markers.cluster))) clusters")

step("Annotating clusters using PanglaoDB...")
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
        title="UMAP — K-Means Clusters (k=8)"), "umap_kmeans.png")

step("UMAP: Graph clusters (Label Propagation)")
savefig(dim_plot(obj, reduction="umap", group="graph_cluster",
        title="UMAP — Graph Clusters"), "umap_graph_cluster.png")

step("UMAP: Cell type annotations")
savefig(dim_plot(obj, reduction="umap", group="cell_type",
        title="UMAP — Cell Type Annotations"), "umap_annotated.png")

step("Pseudobulk heatmap: cell type signatures")
pb_celltype = get_pseudobulk!(obj, "cell_type")
CSV.write("pseudobulk_counts.csv", pb_celltype)
savefig(plot_pseudobulk_heatmap(pb_celltype, ranking=:range, n_genes=40,
        title="Cell Type Signature Genes"), "celltype_heatmap.png")

step("Volcano plot (cluster 1 markers)")
target_cluster = string(unique(markers.cluster)[1])
cluster_mrks = filter(row -> string(row.cluster) == target_cluster, markers)
savefig(volcano_plot(cluster_mrks, logfc_thresh=0.25, pval_thresh=0.05, top_n=15),
    "volcano_cluster_$(target_cluster).png")

step("PAGA connectivity graph (cell types)")
savefig(paga_plot(obj, reduction="umap", cluster_col="cell_type",
        min_conn_weight=0.02, title="Lineage Graph — Cell Types"),
    "paga_celltype_graph.png")

# Top marker genes per cluster — violin + UMAP feature plot
step("Top 3 marker genes per cluster")
top3 = combine(groupby(markers, :cluster)) do df
    sort(df, :avg_log2FC, rev=true)[1:min(3, nrow(df)), :]
end
genes_to_plot = unique(top3.gene)
println("  → Top genes: " * join(genes_to_plot, ", "))

savefig(violin_plot(obj, features=genes_to_plot, group="cell_type"),
    "top_markers_violin.png")

if !isempty(genes_to_plot)
    best_gene = genes_to_plot[1]
    step("UMAP feature plot: $best_gene")
    savefig(feature_plot(obj, feature=best_gene,
            title="Top Marker: $best_gene"), "top_marker_umap_$(best_gene).png")
end

println("\n  ✓ Standard plots saved.")

# ==============================================================================
# 4. TRAJECTORY INFERENCE
# ==============================================================================
section("4. TRAJECTORY INFERENCE")

step("Building Diffusion Map (10 components)...")
run_diffusion_map!(obj, n_components=10)

# Auto-select biologically meaningful root for breast cancer
ROOT_CANDIDATES = [
    "Epithelial cells", "Luminal cells", "Basal cells",
    "Progenitor cells", "Stem cells"
]
available_types = unique(obj.meta_data.cell_type)
root_type = let found = filter(t -> t in available_types, ROOT_CANDIDATES)
    isempty(found) ? available_types[1] : first(found)
end
step("Root cell type selected: '$root_type'")

step("Computing graph-based pseudotime from root...")
run_pseudotime!(obj, root_type, cluster_col="cell_type", method="graph")

step("Plotting pseudotime on UMAP...")
savefig(feature_plot(obj, :pseudotime,
        title="Pseudotime from '$root_type'"), "umap_pseudotime.png")



step("Generating pseudotime gene dynamics heatmap...")
obj.meta_data[!, :pseudotime_bin] = map(obj.meta_data.pseudotime) do t
    t < 0.33 ? "Early" : t < 0.66 ? "Mid" : "Late"
end
pb_time = get_pseudobulk!(obj, "pseudotime_bin")
savefig(plot_pseudobulk_heatmap(pb_time, n_genes=40, ranking=:range,
        title="Gene Dynamics Over Pseudotime"), "pseudotime_heatmap.png")

println("\n  ✓ Trajectory plots saved.")

# ==============================================================================
# 5. TRAJECTORY UNCERTAINTY FRAMEWORK (TES + TDS + LPS)
# ==============================================================================
section("5. TRAJECTORY UNCERTAINTY FRAMEWORK")

step("Computing TES + TDS + LPS trajectory uncertainty metrics...")
trajectory_uncertainty!(obj, key_prefix="traj_unc")

# Summary statistics
println("\n  Uncertainty metrics summary across $(size(obj.meta_data, 1)) cells:")
for metric in ["traj_unc_tes", "traj_unc_tds", "traj_unc_lps"]
    if metric in names(obj.meta_data)
        vals = obj.meta_data[!, metric]
        println("    $metric → μ=$(round(mean(vals), digits=3))   max=$(round(maximum(vals), digits=3))   min=$(round(minimum(vals), digits=3))")
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

println("\n=== Biological Interpretation ===")
println("• High TES regions: cells surrounded by mixed developmental stages")
println("• High TDS regions: potential branch points or fate decision zones")
println("• Positive LPS: cells actively progressing along the trajectory")
println("• Negative LPS: cells in early or 'behind' regions relative to neighbors")

step("Saving trajectory uncertainty plots...")
savefig(
    trajectory_uncertainty_plot(obj, title="Trajectory Uncertainty Framework — Breast Cancer"),
    "trajectory_uncertainty.png"
)
println("    ✓ trajectory_uncertainty.png")

step("Exporting full uncertainty metrics table...")
CSV.write("trajectory_uncertainty_metrics.csv",
    select(obj.meta_data, [:cell_type, :pseudotime, :traj_unc_tes, :traj_unc_tds, :traj_unc_lps]))
println("    ✓ trajectory_uncertainty_metrics.csv")



# ==============================================================================
# VALIDATION EXPERIMENTS
# ==============================================================================
section("5b. UNCERTAINTY VALIDATION")

tes_key = "traj_unc_tes"
tds_key = "traj_unc_tds"
lps_key = "traj_unc_lps"

# 1. Top 5% Highlight Plots
step("Generating top 5% highlight plots...")
p_top = highlight_top_uncertainty(obj, top_pct=0.05,
    title="Top 5% Most Uncertain Cells")
savefig(p_top, "top5_uncertainty.png")
println("    ✓ top5_uncertainty.png  (red = high uncertainty)")

# 2. Cell-type Ranking
step("Cell-type level summary...")
cluster_summary = combine(
    groupby(obj.meta_data, :cell_type),
    tes_key => mean => :mean_TES,
    tds_key => mean => :mean_TDS,
    lps_key => mean => :mean_LPS
)
sort!(cluster_summary, :mean_TES, rev=true)
println("\n  Mean Uncertainty by Cell Type (sorted by TES):")
show(cluster_summary, allrows=true)
CSV.write("celltype_uncertainty_summary.csv", cluster_summary)

# 3. Correlation Analysis
step("Correlation between metrics...")
tes = obj.meta_data[!, tes_key]
tds = obj.meta_data[!, tds_key]
lps = obj.meta_data[!, lps_key]

cor_matrix = cor([tes tds lps])
println("\n  Metric Correlations:")
println("    TES ↔ TDS : $(round(cor_matrix[1,2], digits=3))")
println("    TES ↔ LPS : $(round(cor_matrix[1,3], digits=3))")
println("    TDS ↔ LPS : $(round(cor_matrix[2,3], digits=3))")

if abs(cor_matrix[1, 2]) > 0.7
    @warn "TES and TDS are highly correlated. Consider if they capture distinct biology."
end

# 4. Save summary table
summary_table = DataFrame(
    Metric=["TES", "TDS", "LPS"],
    Mean=[mean(tes), mean(tds), mean(lps)],
    Max=[maximum(tes), maximum(tds), maximum(lps)],
    Correlation_TES=[1.0, cor_matrix[1, 2], cor_matrix[1, 3]]
)
CSV.write("uncertainty_correlations.csv", summary_table)


# ==============================================================================
# BIOLOGICAL VALIDATION OF UNCERTAINTY METRICS
# ==============================================================================
section("5c. BIOLOGICAL VALIDATION")

tes_key = "traj_unc_tes"
tds_key = "traj_unc_tds"

# 1. Top 5% TES cells vs Rest (Differential Expression)
step("Finding markers enriched in Top 5% TES cells...")
top5_tes = obj.meta_data[!, tes_key] .>= quantile(obj.meta_data[!, tes_key], 0.95)

# Temporarily add group label
obj.meta_data[!, :top_tes_group] = ifelse.(top5_tes, "Top5_TES", "Background")

markers_tes = find_all_markers(obj, group="top_tes_group", min_pct=0.1, logfc_threshold=0.25)

# Filter to markers of Top5_TES group
top_tes_markers = filter(row -> row.cluster == "Top5_TES", markers_tes)

if nrow(top_tes_markers) > 0
    println("\n  Top 10 markers enriched in high-TES cells:")
    show(top_tes_markers[1:min(10, nrow(top_tes_markers)), :], allrows=true)
    CSV.write("top5_tes_markers.csv", top_tes_markers)
else
    println("  No significant markers found for top TES cells.")
end

# Cleanup temporary column
select!(obj.meta_data, Not(:top_tes_group))

# 2. Cell-type Ranking (already done, just print nicely)
step("Cell-type uncertainty ranking (biological plausibility check)...")
ranking = combine(
    groupby(obj.meta_data, :cell_type),
    tes_key => mean => :mean_TES,
    tds_key => mean => :mean_TDS
)
sort!(ranking, :mean_TES, rev=true)
println("\n  Cell types ranked by mean TES (highest uncertainty first):")
show(ranking, allrows=true)

# 3. Interpretation Summary
println("\n=== Biological Validation Summary ===")
println("• High TES cells enriched for stemness/EMT/activation genes? → Check top5_tes_markers.csv")
println("• Cell type ranking matches known biology?")
println("   Expected high TES: Progenitors, Stem-like, Activated immune, Fibroblasts")
println("   Expected low TES: Mature epithelial, Terminally differentiated cells")
println("• TES-TDS correlation = ", round(cor(obj.meta_data[!, tes_key], obj.meta_data[!, tds_key]), digits=3))

# ==============================================================================
# 6. SUMMARY
# ==============================================================================
section("TEST COMPLETE — OUTPUT FILES")

outputs = [
    ("umap_kmeans.png", "UMAP — K-Means clusters"),
    ("umap_graph_cluster.png", "UMAP — Label Propagation clusters"),
    ("umap_annotated.png", "UMAP — Cell type annotations"),
    ("celltype_heatmap.png", "Pseudobulk cell type signature genes"),
    ("pseudobulk_counts.csv", "Raw pseudobulk counts table"),
    ("volcano_cluster_$(target_cluster).png", "Volcano plot (cluster $target_cluster)"),
    ("paga_celltype_graph.png", "PAGA lineage connectivity graph"),
    ("top_markers_violin.png", "Violin: top 3 markers per cluster"),
    ("top_marker_umap_$(isempty(genes_to_plot) ? "none" : genes_to_plot[1]).png", "UMAP feature plot: top marker gene"),
    ("umap_pseudotime.png", "Pseudotime on UMAP"),
    ("pseudotime_heatmap.png", "Gene dynamics heatmap over pseudotime"),
    ("trajectory_uncertainty.png", "TES + TDS + LPS uncertainty framework"),
    ("trajectory_uncertainty_metrics.csv", "Per-cell TES, TDS, LPS metrics"),
]

for (fname, desc) in outputs
    marker = isfile(fname) ? "✓" : "·"
    println("  $marker  $fname  —  $desc")
end