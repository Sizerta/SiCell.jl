# ==============================================================================
# SiCell — Glioblastoma Pipeline Test
# Covers: QC, Clustering, Annotation, DE, Trajectory Inference,
#         Trajectory Uncertainty Framework (TES + TDS + LPS)
# Dataset: 320k scFFPE Glioblastoma (10x GEM-X FLEX, 16-plex)
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

DATA_PATH = ".\\320k_scFFPE_16-plex_GEM-X_FLEX_Glioblastoma_BC1-2_count_sample_filtered_feature_bc_matrix.h5"

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
run_clustering!(obj, k=12)

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

step("Annotating clusters using cell marker")
db = load_cellmarker(species="Hs")

# Force annotation with relaxed settings
annotate_clusters!(obj, markers, db,
    species="Hs",
    min_score=0.03)


obj.meta_data[!, :cell_type] = replace(
    obj.meta_data.cell_type, "Unknown" => "Hypoxic tumour cells"
)

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

# Auto-select biologically meaningful root for Glioblastoma
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
    trajectory_uncertainty_plot(obj, title="Trajectory Uncertainty Framework — Glioblastoma"),
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

top5_tes = tes .>= quantile(tes, 0.95)

rank2 = combine(
    groupby(DataFrame(
            celltype=obj.meta_data.cell_type,
            top=top5_tes
        ), :celltype),
    :top => mean => :fraction_top5
)

sort!(rank2, :fraction_top5, rev=true)

println("\nTop 5% TES cells:")
show(rank2, allrows=true)

top5_tds = tds .>= quantile(tds, 0.95)

rank3 = combine(
    groupby(DataFrame(
            celltype=obj.meta_data.cell_type,
            top=top5_tds
        ), :celltype),
    :top => mean => :fraction_top5
)

sort!(rank3, :fraction_top5, rev=true)

println("\nTop 5% TdS cells:")
show(rank3, allrows=true)


# ==============================================================================
# NEW ANALYSES: Hypoxia-TES Integration (Analyses 1-3)
# ==============================================================================
section("6. HYPOXIA-TES INTEGRATION ANALYSES")

using HypothesisTests
using EffectSizes
using StatsBase

tes_key = "traj_unc_tes"

# ============================================================================
# Analysis 1: TES in Hypoxic vs Non-Hypoxic Tumour Cells
# ============================================================================
step("Analysis 1: TES comparison (Hypoxic vs Non-hypoxic)")

hypoxic = obj.meta_data.cell_type .== "Hypoxic tumour cells"
non_hypoxic = .!hypoxic

tes_hyp = obj.meta_data[hypoxic, tes_key]
tes_non = obj.meta_data[non_hypoxic, tes_key]

if length(tes_hyp) > 10 && length(tes_non) > 10
    println("\nHypoxic tumour:")
    println("  n = $(length(tes_hyp))")
    println("  mean TES = $(round(mean(tes_hyp), digits=4)) ± $(round(std(tes_hyp), digits=4))")

    println("\nAll other cells:")
    println("  n = $(length(tes_non))")
    println("  mean TES = $(round(mean(tes_non), digits=4)) ± $(round(std(tes_non), digits=4))")

    # Mann–Whitney U test
    mw_test = MannWhitneyUTest(tes_hyp, tes_non)
    pval = pvalue(mw_test)
    p_str = pval < 1e-10 ? "< 1e-10" : string(round(pval, digits=6))
    println("\nMann–Whitney U test: p = $p_str")

    # Cohen's d
    es = CohenD(tes_hyp, tes_non)          # ← fixed
    d_val = effectsize(es)
    effect_interp = abs(d_val) < 0.2 ? "negligible" :
                    abs(d_val) < 0.5 ? "small" :
                    abs(d_val) < 0.8 ? "medium" : "large"
    println("Cohen's d = $(round(d_val, digits=3)) → Effect size: $effect_interp")
else
    println("Insufficient cells in one group.")
end

# ============================================================================
# Analysis 2: TES Marker Genes
# ============================================================================
step("Analysis 2: TES , TDSv marker genes (high-TES cells)")

top5_tes_idx = obj.meta_data[!, tes_key] .>= quantile(obj.meta_data[!, tes_key], 0.95)
obj.meta_data[!, :high_tes] = ifelse.(top5_tes_idx, "High_TES", "Other")

markers_high_tes = find_all_markers(obj, group="high_tes", min_pct=0.1, logfc_threshold=0.25)

high_tes_markers = filter(row -> row.cluster == "High_TES", markers_high_tes)
sort!(high_tes_markers, :avg_log2FC, rev=true)

println("\nTop 15 markers enriched in high-TES cells:")
if nrow(high_tes_markers) > 0
    show(high_tes_markers[1:min(15, nrow(high_tes_markers)), :], allrows=true)
    CSV.write("high_tes_markers.csv", high_tes_markers)
else
    println("No significant markers found.")
end




top5_tds_idx = obj.meta_data[!, tds_key] .>= quantile(obj.meta_data[!, tds_key], 0.95)
obj.meta_data[!, :high_tds] = ifelse.(top5_tds_idx, "High_tds", "Other")

markers_high_tds = find_all_markers(obj, group="high_tds", min_pct=0.1, logfc_threshold=0.25)

high_tds_markers = filter(row -> row.cluster == "High_tds", markers_high_tds)
sort!(high_tds_markers, :avg_log2FC, rev=true)

println("\nTop 15 markers enriched in high-tds cells:")
if nrow(high_tds_markers) > 0
    show(high_tds_markers[1:min(15, nrow(high_tds_markers)), :], allrows=true)
    CSV.write("high_tds_markers.csv", high_tds_markers)
else
    println("No significant markers found.")
end

# ============================================================================
# Analysis 2: Biological Score Decomposition of TES
# ============================================================================

section("6. HYPOXIA-TES INTEGRATION ANALYSES")

step("Analysis 2: Computing biological signature scores")

# ── Robust gene-name accessor ────────────────────────────────────────────
function get_gene_names(obj)
    vd = obj.var_data
    vd_cols = Symbol.(names(vd))

    for col in [:gene, :gene_name, :gene_names, :symbol, :symbols,
        :features, :feature, :feature_name, :name, :gene_id, :index, :x]
        if col in vd_cols
            println("    → Gene names found in column: '$(String(col))'")
            return String.(vd[!, col])
        end
    end

    for col in names(vd)
        if eltype(vd[!, col]) <: AbstractString
            sample = vd[1:min(5, nrow(vd)), col]
            if any(v -> occursin(r"[A-Za-z]", string(v)), sample)
                println("    → Auto-detected gene name column: '$col'")
                return String.(vd[!, col])
            end
        end
    end

    error("Could not find a gene name column in obj.var_data.\nAvailable columns: $(names(vd))")
end

# ── Score function ──────────────────────────────────────────────────────
function signature_score(obj, genes)
    gene_names = get_gene_names(obj)

    X = obj.norm_data
    if isnothing(X)
        @warn "obj.norm_data is nothing — falling back to obj.counts."
        X = obj.counts
    end

    gene_upper = uppercase.(gene_names)
    found_idx = Int[]
    matched_genes = String[]
    missing_genes = String[]

    for g in genes
        g_upper = uppercase(g)
        hits = findall(==(g_upper), gene_upper)
        if !isempty(hits)
            push!(found_idx, hits[1])
            push!(matched_genes, g)
        else
            push!(missing_genes, g)
        end
    end

    if !isempty(missing_genes)
        println("    ⚠ Genes not found: " * join(missing_genes, ", "))
    end

    if isempty(found_idx)
        error("NONE of the requested signature genes were found! Requested: $(genes)")
    end

    println("    ✓ Matched $(length(found_idx))/$(length(genes)) genes: " * join(matched_genes, ", "))

    # The semicolon at the end of the next line prevents printing 33,000 numbers!
    return vec(Array(mean(X[found_idx, :], dims=1)))
end

# ── Signature gene sets ─────────────────────────────────────────────────
hypoxia_genes = ["CA9", "VEGFA", "SLC2A1", "LDHA", "PGK1", "BNIP3", "HK2", "ENO1", "ALDOA"]
stemness_genes = ["CD44", "SOX2", "NES", "PROM1", "OLIG2", "ID1", "ID2", "NANOG"]
cycle_genes = ["MKI67", "TOP2A", "PCNA", "CCNB1", "CDK1", "CENPF"]
emt_genes = ["VIM", "FN1", "CD44", "CHI3L1", "SERPINE1", "COL1A1", "COL1A2"]

# ── Step 1: Compute scores ──────────────────────────────────────────────
step("Step 1 — Computing signature scores per cell")

# Semicolons suppress the 33,000-number dump in the REPL
hypoxia_score = signature_score(obj, hypoxia_genes);
stemness_score = signature_score(obj, stemness_genes);
cycle_score = signature_score(obj, cycle_genes);
emt_score = signature_score(obj, emt_genes);

println("\n  Score distributions:")
for (name, scores) in [("Hypoxia", hypoxia_score),
    ("Stemness", stemness_score),
    ("Cell Cycle", cycle_score),
    ("EMT", emt_score)]
    println("    $name → μ=$(round(mean(scores), digits=4))  σ=$(round(std(scores), digits=4))  " *
            "min=$(round(minimum(scores), digits=4))  max=$(round(maximum(scores), digits=4))")
end

# ── Step 2: Fit linear model ────────────────────────────────────────────
step("Step 2 — Fitting linear model: TES ~ hypoxia + stemness + cellcycle + emt")

using GLM

tes = obj.meta_data[!, tes_key];

df_model = DataFrame(
    TES=tes,
    hypoxia=hypoxia_score,
    stemness=stemness_score,
    cellcycle=cycle_score,
    emt=emt_score
);

model = lm(@formula(TES ~ hypoxia + stemness + cellcycle + emt), df_model);

r2_val = r2(model);
adj_r2 = adjr2(model);
pct_explained = round(r2_val * 100, digits=1);

println("\n  Model fit:")
println("    R²      = $(round(r2_val, digits=4))")
println("    R²_adj  = $(round(adj_r2, digits=4))")
println("    → Known biology explains $(pct_explained)% of TES variance")

println("\n  Coefficient details:")
coeftbl = coeftable(model);
for i in 1:size(coeftbl.cols[1], 1)
    term = coeftbl.rownms[i]
    coef_v = coeftbl.cols[1][i]
    se_val = coeftbl.cols[2][i]
    p_val = coeftbl.cols[4][i]
    sig = p_val < 0.001 ? "***" :
          p_val < 0.01 ? "**" :
          p_val < 0.05 ? "*" : ""
    println("    $(rpad(term, 14)) β=$(round(coef_v, digits=4))  SE=$(round(se_val, digits=4))  p=$(round(p_val, digits=6)) $sig")
end

# ── Step 3: Calculate residual TES ──────────────────────────────────────
step("Step 3 — Calculating residual TES (unexplained uncertainty)")

df_model.TES_pred = predict(model);
df_model.TES_residual = df_model.TES .- df_model.TES_pred;

println("\n  Residual TES distribution:")
println("    μ   = $(round(mean(df_model.TES_residual), digits=4))")
println("    σ   = $(round(std(df_model.TES_residual), digits=4))")
println("    min = $(round(minimum(df_model.TES_residual), digits=4))")
println("    max = $(round(maximum(df_model.TES_residual), digits=4))")
println()
println("  Interpretation:")
println("    Positive residual → Actual TES >> Expected (potential novel biology)")
println("    Negative residual → Actual TES << Expected (well-explained cells)")

# ── Step 4: Identify top 5% residual-high cells ─────────────────────────
step("Step 4 — Identifying top 5% residual-high cells")

cutoff = quantile(df_model.TES_residual, 0.95);
high_residual = df_model.TES_residual .>= cutoff;

println("  Cutoff (95th percentile): $(round(cutoff, digits=4))")
println("  High-residual cells: $(sum(high_residual)) / $(length(high_residual))")

obj.meta_data[!, :TES_residual] = df_model.TES_residual;
obj.meta_data[!, :high_TES_residual] = high_residual;

# ── Step 5: Characterize residual-high cells ─────────────────────────────
step("Step 5 — Characterizing residual-high cells")

step("5a. Cell-type enrichment among high-residual cells")

enrichment = combine(
    groupby(DataFrame(
            celltype=obj.meta_data.cell_type,
            high=high_residual
        ), :celltype),
    :high => mean => :fraction_high_residual,
    nrow => :total_cells
);
sort!(enrichment, :fraction_high_residual, rev=true);

println("\n  Cell-type enrichment in top 5% residual TES:")
show(enrichment, allrows=true)
CSV.write("residual_TES_celltype_enrichment.csv", enrichment)

step("5b. Finding markers enriched in high-residual TES cells")

obj.meta_data[!, :residual_group] = ifelse.(high_residual, "HighResidual_TES", "Background");

markers_residual = find_all_markers(obj, group="residual_group", min_pct=0.1, logfc_threshold=0.25);

high_res_markers = filter(
    row ->
        row.cluster == "HighResidual_TES" &&
            row.p_val_adj < 0.05 &&
            row.avg_log2FC > 0.25,
    markers_residual
)

if nrow(high_res_markers) > 0
    println("\n  Top 15 markers enriched in high-residual TES cells:")
    show(high_res_markers[1:min(15, nrow(high_res_markers)), :], allrows=true)
    CSV.write("high_residual_TES_markers.csv", high_res_markers)

    savefig(
        volcano_plot(high_res_markers, logfc_thresh=0.25, pval_thresh=0.05, top_n=15,
            title="DE: High Residual TES vs Background"),
        "volcano_high_residual_TES.png"
    )
else
    println("  No significant markers found for high-residual TES cells.")
end

select!(obj.meta_data, Not(:residual_group));

step("5c. Signature score profiles of high-residual cells")

for (name, scores) in [("Hypoxia", hypoxia_score),
    ("Stemness", stemness_score),
    ("Cell Cycle", cycle_score),
    ("EMT", emt_score)]
    high_vals = scores[high_residual]
    bg_vals = scores[.!high_residual]
    mw = MannWhitneyUTest(high_vals, bg_vals)
    local p_str = pvalue(mw) < 1e-10 ? "< 1e-10" : string(round(pvalue(mw), digits=6))
    println("    $name → HighRes: μ=$(round(mean(high_vals), digits=4))  " *
            "BG: μ=$(round(mean(bg_vals), digits=4))  p=$p_str")
end

step("5d. Visualizing residual TES on UMAP")

obj.meta_data[!, :TES_residual_bin] = map(df_model.TES_residual) do r
    r >= cutoff ? "Top5% Residual" :
    r <= quantile(df_model.TES_residual, 0.05) ? "Bottom5% Residual" :
    "Middle 90%"
end;

savefig(
    dim_plot(obj, reduction="umap", group="TES_residual_bin",
        title="Residual TES — Cells Beyond Known Biology"),
    "umap_residual_TES.png"
)

select!(obj.meta_data, Not(:TES_residual_bin));

# ── Summary ──────────────────────────────────────────────────────────────
println("\n$( "="^60 )")
println("  RESIDUAL TES ANALYSIS SUMMARY")
println("$( "="^60 )")
println("  Known biology (hypoxia + stemness + cycle + EMT) explains:")
println("    → $(pct_explained)% of TES variance (R² = $(round(r2_val, digits=4)))")
println("    → $(round(100 - r2_val*100, digits=1))% remains unexplained")
println()
println("  Top 5% residual-high cells are enriched in:")

if nrow(enrichment) > 0
    top_enriched = filter(row -> row.fraction_high_residual > mean(enrichment.fraction_high_residual), enrichment)
    for row in eachrow(top_enriched)
        println("    • $(row.celltype) ($(round(row.fraction_high_residual * 100, digits=1))% of type is high-residual)")
    end
end

println()
println("  Key question: Do high-residual cells express novel programs")
println("  not captured by hypoxia, stemness, cell cycle, or EMT?")
println("  → Check high_residual_TES_markers.csv for candidate genes")
println("$( "="^60 )")

model_summary = DataFrame(
    term=coeftbl.rownms,
    estimate=coeftbl.cols[1],
    std_error=coeftbl.cols[2],
    p_value=coeftbl.cols[4],
    r_squared=fill(r2_val, length(coeftbl.rownms)),
    adj_r2=fill(adj_r2, length(coeftbl.rownms))
);
CSV.write("TES_model_coefficients.csv", model_summary);
println("  ✓ TES_model_coefficients.csv saved")

CSV.write("residual_DE.csv", high_res_markers)















# ============================================================================
# TDS Residual Analysis
# ============================================================================

step("Fitting linear model: TDS ~ hypoxia + stemness + cellcycle + emt")

using GLM

tds = obj.meta_data[!, tds_key];

df_model_tds = DataFrame(
    TDS       = tds,
    hypoxia   = hypoxia_score,
    stemness  = stemness_score,
    cellcycle = cycle_score,
    emt       = emt_score
);

model_tds = lm(@formula(TDS ~ hypoxia + stemness + cellcycle + emt), df_model_tds);

r2_val = r2(model_tds);
adj_r2 = adjr2(model_tds);
pct_explained = round(r2_val * 100, digits=1);

println("\n  Model fit:")
println("    R²      = $(round(r2_val, digits=4))")
println("    R²_adj  = $(round(adj_r2, digits=4))")
println("    → Known biology explains $(pct_explained)% of TDS variance")

println("\n  Coefficient details:")
coeftbl = coeftable(model_tds);
for i in 1:size(coeftbl.cols[1], 1)
    term    = coeftbl.rownms[i]
    coef_v  = coeftbl.cols[1][i]
    se_val  = coeftbl.cols[2][i]
    p_val   = coeftbl.cols[4][i]
    sig     = p_val < 0.001 ? "***" :
              p_val < 0.01  ? "**"  :
              p_val < 0.05  ? "*"   : ""
    println("    $(rpad(term, 14)) β=$(round(coef_v, digits=4))  SE=$(round(se_val, digits=4))  p=$(round(p_val, digits=6)) $sig")
end

# ── Step 3: Calculate residual TDS ──────────────────────────────────────
step("Calculating residual TDS (unexplained uncertainty)")

df_model_tds.TDS_pred     = predict(model_tds);
df_model_tds.TDS_residual = df_model_tds.TDS .- df_model_tds.TDS_pred;

println("\n  Residual TDS distribution:")
println("    μ   = $(round(mean(df_model_tds.TDS_residual), digits=4))")
println("    σ   = $(round(std(df_model_tds.TDS_residual), digits=4))")
println("    min = $(round(minimum(df_model_tds.TDS_residual), digits=4))")
println("    max = $(round(maximum(df_model_tds.TDS_residual), digits=4))")
println()
println("  Interpretation:")
println("    Positive residual → Actual TDS >> Expected (potential novel biology)")
println("    Negative residual → Actual TDS << Expected (well-explained cells)")

# ── Step 4: Identify top 5% residual-high cells ─────────────────────────
step("Identifying top 5% residual-high TDS cells")

cutoff_tds = quantile(df_model_tds.TDS_residual, 0.95);
high_residual_tds = df_model_tds.TDS_residual .>= cutoff_tds;

println("  Cutoff (95th percentile): $(round(cutoff_tds, digits=4))")
println("  High-residual cells: $(sum(high_residual_tds)) / $(length(high_residual_tds))")

obj.meta_data[!, :TDS_residual]     = df_model_tds.TDS_residual;
obj.meta_data[!, :high_TDS_residual] = high_residual_tds;

# ── Step 5: Characterize residual-high cells ─────────────────────────────
step("Characterizing residual-high TDS cells")

step("5a. Cell-type enrichment among high-residual TDS cells")

enrichment_tds = combine(
    groupby(DataFrame(
        celltype = obj.meta_data.cell_type,
        high     = high_residual_tds
    ), :celltype),
    :high => mean => :fraction_high_residual,
    nrow => :total_cells
);
sort!(enrichment_tds, :fraction_high_residual, rev=true);

println("\n  Cell-type enrichment in top 5% residual TDS:")
show(enrichment_tds, allrows=true)
CSV.write("residual_TDS_celltype_enrichment.csv", enrichment_tds)

step("5b. Finding markers enriched in high-residual TDS cells")

obj.meta_data[!, :residual_group] = ifelse.(high_residual_tds, "HighResidual_TDS", "Background");

markers_residual_tds = find_all_markers(obj, group="residual_group", min_pct=0.1, logfc_threshold=0.25);

high_res_markers_tds = filter(
    row ->
        row.cluster == "HighResidual_TDS" &&
            row.p_val_adj < 0.05 &&
            row.avg_log2FC > 0.25,
    markers_residual_tds
);

if nrow(high_res_markers_tds) > 0
    println("\n  Top 15 markers enriched in high-residual TDS cells:")
    show(high_res_markers_tds[1:min(15, nrow(high_res_markers_tds)), :], allrows=true)
    CSV.write("high_residual_TDS_markers.csv", high_res_markers_tds)

    savefig(
        volcano_plot(high_res_markers_tds, logfc_thresh=0.25, pval_thresh=0.05, top_n=15,
            title="DE: High Residual TDS vs Background"),
        "volcano_high_residual_TDS.png"
    )
else
    println("  No significant markers found for high-residual TDS cells.")
end

# Clean up temporary column
select!(obj.meta_data, Not(:residual_group));

step("5c. Signature score profiles of high-residual TDS cells")

for (name, scores) in [("Hypoxia", hypoxia_score),
                        ("Stemness", stemness_score),
                        ("Cell Cycle", cycle_score),
                        ("EMT", emt_score)]
    high_vals = scores[high_residual_tds]
    bg_vals   = scores[.!high_residual_tds]
    mw = MannWhitneyUTest(high_vals, bg_vals)
    local p_str = pvalue(mw) < 1e-10 ? "< 1e-10" : string(round(pvalue(mw), digits=6))
    println("    $name → HighRes: μ=$(round(mean(high_vals), digits=4))  " *
            "BG: μ=$(round(mean(bg_vals), digits=4))  p=$p_str")
end

step("5d. Visualizing residual TDS on UMAP")

obj.meta_data[!, :TDS_residual_bin] = map(df_model_tds.TDS_residual) do r
    r >= cutoff_tds ? "Top5% Residual" :
    r <= quantile(df_model_tds.TDS_residual, 0.05) ? "Bottom5% Residual" :
    "Middle 90%"
end;

savefig(
    dim_plot(obj, reduction="umap", group="TDS_residual_bin",
        title="Residual TDS — Cells Beyond Known Biology"),
    "umap_residual_TDS.png"
)

select!(obj.meta_data, Not(:TDS_residual_bin));

# ── Summary ──────────────────────────────────────────────────────────────
println("\n$( "="^60 )")
println("  RESIDUAL TDS ANALYSIS SUMMARY")
println("$( "="^60 )")
println("  Known biology (hypoxia + stemness + cycle + EMT) explains:")
println("    → $(pct_explained)% of TDS variance (R² = $(round(r2_val, digits=4)))")
println("    → $(round(100 - r2_val*100, digits=1))% remains unexplained")
println()
println("  Top 5% residual-high TDS cells are enriched in:")

if nrow(enrichment_tds) > 0
    top_enriched = filter(row -> row.fraction_high_residual > mean(enrichment_tds.fraction_high_residual), enrichment_tds)
    for row in eachrow(top_enriched)
        println("    • $(row.celltype) ($(round(row.fraction_high_residual * 100, digits=1))% of type is high-residual)")
    end
end

println()
println("  Key question: Do high-residual cells express novel programs")
println("  not captured by hypoxia, stemness, cell cycle, or EMT?")
println("  → Check high_residual_TDS_markers.csv for candidate genes")
println("$( "="^60 )")

model_summary_tds = DataFrame(
    term       = coeftbl.rownms,
    estimate   = coeftbl.cols[1],
    std_error  = coeftbl.cols[2],
    p_value    = coeftbl.cols[4],
    r_squared  = fill(r2_val, length(coeftbl.rownms)),
    adj_r2     = fill(adj_r2, length(coeftbl.rownms))
);
CSV.write("TDS_model_coefficients.csv", model_summary_tds);
println("  ✓ TDS_model_coefficients.csv saved")

# Fixed: Save the actual TDS residual markers, not the old 'markers' variable
CSV.write("residual_TDS_DE.csv", high_res_markers_tds);
println("  ✓ residual_TDS_DE.csv saved")





































# ==============================================================================
# STEP 1 — ISOLATE THE UNKNOWN CLUSTER
# ==============================================================================
section("1. ISOLATE UNKNOWN CLUSTER")

# Boolean mask for Unknown cells
unknown_mask = obj.meta_data.cell_type .== "Unknown"
n_unknown = sum(unknown_mask)
n_total = nrow(obj.meta_data)

println("  Unknown cells : $n_unknown / $n_total  ($(round(100*n_unknown/n_total, digits=1))%)")

# Pull TES scores for Unknown vs rest
tes_unknown = obj.meta_data[unknown_mask, :traj_unc_tes]
tes_rest = obj.meta_data[.!unknown_mask, :traj_unc_tes]
pseudotime_unk = obj.meta_data[unknown_mask, :pseudotime]

println("  TES  Unknown : μ=$(round(mean(tes_unknown),    digits=3))  median=$(round(median(tes_unknown),    digits=3))")
println("  TES  Rest    : μ=$(round(mean(tes_rest),       digits=3))  median=$(round(median(tes_rest),       digits=3))")
println("  Pseudotime Unknown: μ=$(round(mean(pseudotime_unk), digits=3))  range=$(round(minimum(pseudotime_unk),digits=2))–$(round(maximum(pseudotime_unk),digits=2))")

# ==============================================================================
# STEP 2 — DIFFERENTIAL EXPRESSION: Unknown vs All Other Cell Types
# ==============================================================================
section("2. DIFFERENTIAL EXPRESSION — Unknown vs Rest")

step("Running DE: Unknown vs all other cell types...")

# Add temporary binary label
obj.meta_data[!, :unknown_group] = ifelse.(unknown_mask, "Unknown", "Other")

de_unknown = find_all_markers(
    obj,
    group="unknown_group",
    min_pct=0.10,
    logfc_threshold=0.25
)

# Keep only Unknown-enriched genes
de_up = filter(row -> row.cluster == "Unknown", de_unknown)
sort!(de_up, :avg_log2FC, rev=true)

println("\n  Top 20 genes upregulated in Unknown cells:")
show(de_up[1:min(20, nrow(de_up)), :], allrows=true)

CSV.write("unknown_cluster_DE.csv", de_up)
println("\n  ✓ Full DE table → unknown_cluster_DE.csv")

# Clean up temp column
select!(obj.meta_data, Not(:unknown_group))

# ==============================================================================
# STEP 3 — SCORE AGAINST KNOWN GSC MARKER PANELS
# ==============================================================================

# ==============================================================================
# STEP 4 — CHECK INDIVIDUAL KEY MARKERS IN UNKNOWN CELLS
# ==============================================================================
section("4. KEY MARKER EXPRESSION IN UNKNOWN CELLS")

key_markers = ["CD44", "SOX2", "NES", "PROM1", "OLIG2", "ID4", "EGFR", "PDGFRA", "VIM", "ZEB1"]

step("Checking expression of key GSC/GBM markers in Unknown vs rest...")

# Get available genes (not all panels genes may be in the dataset)
available_genes = filter(g -> g in obj.var_names, key_markers)
missing_genes = filter(g -> g ∉ obj.var_names, key_markers)

if !isempty(missing_genes)
    println("  ⚠ Genes not in dataset: $(join(missing_genes, ", "))")
end
println("  ✓ Genes found: $(join(available_genes, ", "))")

if !isempty(available_genes)
    # Expression per cell type for each key marker
    expr_df = DataFrame(cell_type=obj.meta_data.cell_type)
    for gene in available_genes
        idx = findfirst(==(gene), obj.var_names)
        expr_df[!, gene] = Vector(obj.counts[idx, :])
    end

    marker_by_type = combine(
        groupby(expr_df, :cell_type),
        [g => mean => g for g in available_genes]...
    )
    sort!(marker_by_type, :CD44, rev=true)

    println("\n  Mean expression of key markers by cell type:")
    show(marker_by_type, allrows=true)
    CSV.write("unknown_marker_expression.csv", marker_by_type)
end

# ==============================================================================
# STEP 5 — VISUALIZATIONS
# ==============================================================================
section("5. VISUALIZATIONS")

# 5a. UMAP highlighting Unknown cells
step("UMAP: Unknown cells highlighted...")
obj.meta_data[!, :is_unknown] = ifelse.(unknown_mask, "Unknown", "Other")
savefig(
    dim_plot(obj, reduction="umap", group="is_unknown",
        title="Unknown Cluster Location on UMAP",
        colors=["#E84040", "#CCCCCC"]),
    "unknown_umap_highlight.png"
)
println("    ✓ unknown_umap_highlight.png")

# 5b. GSC score on UMAP
step("UMAP: GSC score overlay...")
savefig(
    feature_plot(obj, :GSC_score,
        title="GSC Score (CD44, SOX2, NES, PROM1, OLIG2...)"),
    "unknown_gsc_score_umap.png"
)
println("    ✓ unknown_gsc_score_umap.png")

# 5c. TES score on UMAP
step("UMAP: TES score overlay...")
savefig(
    feature_plot(obj, :traj_unc_tes,
        title="TES Score — Trajectory Uncertainty"),
    "unknown_tes_umap.png"
)
println("    ✓ unknown_tes_umap.png")

# 5d. Feature plots for top individual markers
step("Feature plots: CD44, SOX2, ID4, EGFR...")
for gene in filter(g -> g in obj.var_names, ["CD44", "SOX2", "ID4", "EGFR", "PDGFRA", "NES"])
    savefig(
        feature_plot(obj, gene, title="$gene expression"),
        "unknown_marker_$(gene).png"
    )
    println("    ✓ unknown_marker_$(gene).png")
end

# 5e. Violin: GSC score by cell type
step("Violin: GSC score across cell types...")
savefig(
    violin_plot(obj, features=["GSC_score", "MES_score", "Proneural_score"],
        group="cell_type",
        title="Stem Cell Program Scores by Cell Type"),
    "unknown_score_violin.png"
)
println("    ✓ unknown_score_violin.png")

# 5f. Dot plot: key markers in Unknown vs rest
step("Dot plot: key markers in Unknown vs other cell types...")
if !isempty(available_genes)
    savefig(
        dot_plot(obj, features=available_genes, group="cell_type",
            title="Key GBM Markers: Expression × % Cells"),
        "unknown_dotplot.png"
    )
    println("    ✓ unknown_dotplot.png")
end

# 5g. Volcano plot for Unknown DE
step("Volcano: Unknown vs rest...")
savefig(
    volcano_plot(de_up, logfc_thresh=0.5, pval_thresh=0.05, top_n=20,
        title="Unknown Cluster — Upregulated Genes"),
    "unknown_volcano.png"
)
println("    ✓ unknown_volcano.png")

# ==============================================================================
# STEP 6 — TES CORRELATION WITH GSC SCORE
# ==============================================================================
section("6. TES ↔ GSC SCORE CORRELATION")

step("Is high TES driven by stemness?")

tes_vals = obj.meta_data[!, :traj_unc_tes]
gsc_vals = obj.meta_data[!, :GSC_score]
mes_vals = obj.meta_data[!, :MES_score]

cor_tes_gsc = cor(tes_vals, gsc_vals)
cor_tes_mes = cor(tes_vals, mes_vals)

println("\n  TES ↔ GSC score correlation : $(round(cor_tes_gsc, digits=3))")
println("  TES ↔ MES score correlation : $(round(cor_tes_mes, digits=3))")

if cor_tes_gsc > 0.3
    println("\n  ✅ TES correlates with GSC stemness — supports GSC interpretation")
elseif cor_tes_gsc > 0.15
    println("\n  ⚠  Moderate TES-GSC correlation — partial stemness signal")
else
    println("\n  ✗  Weak TES-GSC correlation — Unknown cluster may not be GSCs")
end

# Scatter plot: TES vs GSC score, colored by cell type
step("Scatter: TES vs GSC score...")
scatter_df = DataFrame(
    TES=tes_vals,
    GSC_score=gsc_vals,
    cell_type=obj.meta_data.cell_type
)
CSV.write("tes_gsc_scatter.csv", scatter_df)

# ==============================================================================
# STEP 7 — PSEUDOTIME DISTRIBUTION OF UNKNOWN CELLS
# ==============================================================================
section("7. PSEUDOTIME POSITION OF UNKNOWN CELLS")

step("Where do Unknown cells sit in pseudotime?")

pt_summary = combine(
    groupby(obj.meta_data, :cell_type),
    :pseudotime => mean => :mean_PT,
    :pseudotime => median => :median_PT,
    :pseudotime => minimum => :min_PT,
    :pseudotime => maximum => :max_PT
)
sort!(pt_summary, :mean_PT)

println("\n  Pseudotime position by cell type (early → late):")
show(pt_summary, allrows=true)
CSV.write("unknown_pseudotime_position.csv", pt_summary)

# Interpretation
unk_pt = pt_summary[pt_summary.cell_type.=="Unknown", :mean_PT]
if !isempty(unk_pt)
    overall_mean_pt = mean(obj.meta_data.pseudotime)
    if unk_pt[1] < overall_mean_pt
        println("\n  → Unknown cells sit EARLY in pseudotime (mean PT=$(round(unk_pt[1],digits=3)) vs overall $(round(overall_mean_pt,digits=3)))")
        println("    Consistent with stem/progenitor identity (GSC hypothesis supported)")
    else
        println("\n  → Unknown cells sit LATE in pseudotime")
        println("    May represent a terminal or aberrant state rather than progenitors")
    end
end

# ==============================================================================
# STEP 8 — FINAL VERDICT
# ==============================================================================
section("8. INTERPRETATION SUMMARY")

println("""
  ┌─────────────────────────────────────────────────────────┐
  │  UNKNOWN CLUSTER INVESTIGATION — RESULTS SUMMARY        │
  └─────────────────────────────────────────────────────────┘

  Check the following outputs to reach your conclusion:

  EVIDENCE FOR GSC IDENTITY:
  ✓ unknown_cluster_DE.csv        → top upregulated genes
    → Look for: CD44, SOX2, NES, PROM1, OLIG2, ID4, EGFR
  ✓ unknown_gsc_scores.csv        → GSC score vs other cell types
    → Unknown should rank #1 if GSCs
  ✓ unknown_marker_expression.csv → per-gene expression by cell type
  ✓ unknown_pseudotime_position.csv → should be early if stem-like

  VISUALIZATIONS:
  ✓ unknown_umap_highlight.png    → spatial location on UMAP
  ✓ unknown_gsc_score_umap.png   → does GSC score overlap Unknown?
  ✓ unknown_tes_umap.png         → does TES overlap Unknown?
  ✓ unknown_marker_CD44.png      → CD44 (key GSC marker)
  ✓ unknown_marker_SOX2.png      → SOX2 (key GSC marker)
  ✓ unknown_score_violin.png     → GSC/MES/Proneural by cell type
  ✓ unknown_volcano.png          → Unknown vs rest DE

  DECISION LOGIC:
  If GSC score highest in Unknown + CD44/SOX2/ID4 upregulated
  + early pseudotime → RELABEL as "Glioma Stem Cell (GSC)"

  If mesenchymal markers dominant (VIM, ZEB1, FN1)
  → RELABEL as "Mesenchymal GSC" or "Therapy-Resistant GSC"

  If mixed lineage (neural + immune + other)
  → Likely DOUBLETS → filter out before final analysis
""")

# ==============================================================================
# OPTIONAL: RELABEL IF EVIDENCE IS CLEAR
# ==============================================================================
section("9. RELABELING (run after reviewing evidence)")

# Uncomment and run after you've checked the outputs above:

# -- Option A: relabel all Unknown as GSC
# obj.meta_data[!, :cell_type_refined] = replace(
#     obj.meta_data.cell_type, "Unknown" => "Glioma Stem Cell (GSC)"
# )

# -- Option B: relabel only high-GSC-score Unknown cells as GSC
# gsc_threshold = quantile(obj.meta_data.GSC_score, 0.75)
# obj.meta_data[!, :cell_type_refined] = map(
#     (ct, gs) -> (ct == "Unknown" && gs >= gsc_threshold) ? "Glioma Stem Cell (GSC)" : ct,
#     obj.meta_data.cell_type, obj.meta_data.GSC_score
# )

# -- After relabeling, rerun TES summary:
# rerun_summary = combine(
#     groupby(obj.meta_data, :cell_type_refined),
#     :traj_unc_tes => mean => :mean_TES,
#     :traj_unc_tds => mean => :mean_TDS,
#     :GSC_score    => mean => :mean_GSC
# )
# sort!(rerun_summary, :mean_TES, rev=true)
# println("Updated TES ranking after relabeling:")
# show(rerun_summary, allrows=true)

println("\n  ✓ Investigation complete. Review outputs then decide on relabeling.")





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