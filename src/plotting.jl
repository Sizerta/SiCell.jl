# ==============================================================================
# SiCell Plotting Module
# Contains: Dimensional Reduction, Feature Plots, Volcano, Violin, and PAGA
# ==============================================================================

# Publication-quality defaults
default(
    fontfamily="Helvetica",
    grid=false,
    framestyle=:box,
    legendfontsize=9,
    guidefontsize=11,
    tickfontsize=9,
    titlefontsize=12,
    dpi=300,
    markersize=1.5,
    alpha=0.7
)

# Extended categorical palette — 28 perceptually distinct colors.
# Derived from the Zeileis/Okabe-Ito + tableau sets, optimized for
# color-blind accessibility and print reproduction.
const SCANPY_COLORS = [
    "#1f78b4", "#ff7f00", "#33a02c", "#e31a1c",
    "#00bcd4", "#f4d03f", "#9c27b0", "#ff69b4",
    "#795548", "#607d8b", "#00897b", "#43a047",
    "#e91e63", "#f06292", "#29b6f6", "#ffa726",
    "#66bb6a", "#d4ac0d", "#ffca28", "#546e7a",
    "#a1887f", "#ce93d8", "#4dd0e1", "#1565c0",
    "#7e57c2", "#8d6e63", "#f48fb1", "#90a4ae"
]

# ------------------------------------------------------------------------------
# 1. Internal Style Helpers
# ------------------------------------------------------------------------------

function _embedding_style(n::Int)
    if n < 800
        return (ms=4.2, α=0.95, stroke=0.3)
    elseif n < 5_000
        return (ms=2.0, α=0.75, stroke=0.0)
    else
        return (ms=1.1, α=0.55, stroke=0.0)
    end
end

function _clip_limits(v; q=0.005, pad=0.08)
    lo, hi = quantile(v, (q, 1 - q))
    span = hi - lo
    return (lo - pad * span, hi + pad * span)
end

# ------------------------------------------------------------------------------
# 2. Dimensionality Reduction Plot (UMAP/t-SNE)
# ------------------------------------------------------------------------------

"""
    dim_plot(obj; reduction="umap", group="cluster", title="")

Visualizes cells in reduced dimensional space, colored by a metadata column.

Points are shuffled before plotting so no group systematically paints over
another. Axis limits are clipped at the 0.5th/99.5th percentile to prevent
outliers from compressing the main cloud.
"""
function dim_plot(
    obj::SingleCellObject;
    reduction::String="umap",
    group::String="cluster",
    title::String=""
)
    !(reduction in keys(obj.reductions)) &&
        error("Reduction '$reduction' not found.")

    coords = obj.reductions[reduction]
    x, y = coords[:, 1], coords[:, 2]
    n = length(x)
    style = _embedding_style(n)

    xlims = _clip_limits(x)
    ylims = _clip_limits(y)
    ax1 = uppercase(reduction) * " 1"
    ax2 = uppercase(reduction) * " 2"

    groups = group in names(obj.meta_data) ? obj.meta_data[!, group] : nothing

    if isnothing(groups)
        # ── Ungrouped ──────────────────────────────────────────────
        return scatter(x, y,
            color=SCANPY_COLORS[1],
            markersize=style.ms,
            alpha=style.α,
            markerstrokewidth=0,
            xlabel=ax1, ylabel=ax2,
            title=title,
            legend=false,
            xlims=xlims, ylims=ylims,
            size=(820, 640))
    end

    # ── Grouped ────────────────────────────────────────────────────
    unique_labels = sort(unique(groups), by=string)
    n_labels = length(unique_labels)
    color_map = Dict(lbl => SCANPY_COLORS[mod1(i, length(SCANPY_COLORS))]
                     for (i, lbl) in enumerate(unique_labels))

    # Shuffle draw order so no group consistently occludes another
    draw_order = randperm(n)
    xs = x[draw_order]
    ys = y[draw_order]
    gs = groups[draw_order]

    p = plot(size=(820, 640),
        title=title,
        xlabel=ax1, ylabel=ax2,
        xlims=xlims, ylims=ylims,
        framestyle=:box,
        grid=false)

    # Grey background halo — all cells in neutral color first,
    # gives visual context when a group is highlighted
    scatter!(p, xs, ys,
        color="#e8e8e8",
        markersize=style.ms * 0.9,
        alpha=0.35,
        markerstrokewidth=0,
        label="")

    # Colored pass — each group on top of the grey halo
    for lbl in unique_labels
        mask = gs .== lbl
        scatter!(p, xs[mask], ys[mask],
            label=string(lbl),
            color=color_map[lbl],
            markersize=style.ms,
            alpha=style.α,
            markerstrokewidth=0)
    end

    # Legend: outside right for >6 groups, inside top-right otherwise
    legend_pos = n_labels > 6 ? :outertopright : :topright
    plot!(p, legend=legend_pos,
        legendfontsize=8,
        foreground_color_legend=nothing,   # no legend border
        background_color_legend=nothing)

    return p
end

# --------------------------------------------------------------------------
# 3. Feature Plotting
# --------------------------------------------------------------------------

function feature_plot(
    obj::SingleCellObject;
    feature::Union{AbstractString,Symbol},
    reduction::String="umap",
    title::String=""
)
    feature_plot(obj, feature; reduction=reduction, title=title)
end

function feature_plot(
    obj::SingleCellObject,
    feature::Union{AbstractString,Symbol};
    reduction::String="umap",
    title::String=""
)

    feat_name = String(feature)

    # ----------------------------------------------------------------------
    # Check whether feature is metadata or gene
    # ----------------------------------------------------------------------
    if feat_name in String.(names(obj.meta_data))

        col = obj.meta_data[!, Symbol(feat_name)]

        try
            expression = Float64.(col)
        catch
            error("Metadata column '$feat_name' cannot be converted to numeric values.")
        end

    else


        gene_idx = findfirst(obj.var_data.gene_name .== feat_name)
        if isnothing(gene_idx)
            gene_idx = findfirst(obj.var_data.gene_id .== feat_name)
        end
        isnothing(gene_idx) &&
            error("Feature '$feat_name' not found in metadata or gene list.")


        expression = vec(obj.norm_data[gene_idx, :])
        expression = Float64.(expression)

    end

    # ----------------------------------------------------------------------
    # Validate reduction
    # ----------------------------------------------------------------------
    if !(reduction in keys(obj.reductions))
        error("Reduction '$reduction' not found.")
    end

    coords = obj.reductions[reduction]

    x = Vector(coords[:, 1])
    y = Vector(coords[:, 2])

    axis_name = uppercase(reduction)

    # ----------------------------------------------------------------------
    # Color scaling
    # ----------------------------------------------------------------------
    clims = quantile(expression, (0.01, 0.99))

    # ----------------------------------------------------------------------
    # Axis clipping
    # ----------------------------------------------------------------------
    xlims = _clip_limits(x)
    ylims = _clip_limits(y)

    display_title = isempty(title) ? feat_name : title

    # ----------------------------------------------------------------------
    # Adaptive point styling
    # ----------------------------------------------------------------------
    n = length(x)
    style = _embedding_style(n)

    # ----------------------------------------------------------------------
    # 3D reduction
    # ----------------------------------------------------------------------
    if size(coords, 2) >= 3

        z = Vector(coords[:, 3])
        zlims = _clip_limits(z)

        p = scatter(
            x,
            y,
            z,
            zcolor=expression,
            clims=clims,
            xlims=xlims,
            ylims=ylims,
            zlims=zlims,
            markersize=style.ms * 0.85,
            alpha=style.α,
            markerstrokewidth=style.stroke,
            xlabel="$axis_name 1",
            ylabel="$axis_name 2",
            zlabel="$axis_name 3",
            title=display_title,
            legend=false,
            size=(800, 600),
            colorbar_title=feat_name,
            palette=:viridis
        )

    else

        p = scatter(
            x,
            y,
            zcolor=expression,
            clims=clims,
            xlims=xlims,
            ylims=ylims,
            aspect_ratio=1,
            markersize=style.ms,
            alpha=style.α,
            markerstrokewidth=style.stroke,
            xlabel="$axis_name 1",
            ylabel="$axis_name 2",
            title=display_title,
            colorbar=true,
            legend=false,
            size=(800, 600),
            colorbar_title=feat_name,
            palette=:viridis
        )

    end

    return p
end
# ------------------------------------------------------------------------------
# 4. Volcano Plot
# ------------------------------------------------------------------------------



function volcano_plot(markers::DataFrame;
    logfc_thresh::Float64=0.25,
    pval_thresh::Float64=0.05,
    top_n::Int=10,
    title::String="Volcano Plot")

    df = copy(markers)
    df[!, :neg_log_p] = -log10.(max.(df.p_val_adj, 1e-300))
    df[!, :sig] = (df.p_val_adj .< pval_thresh) .&
                  (df.avg_log2FC .> logfc_thresh)

    ns = df[.!df.sig, :]
    sig = df[df.sig, :]

    # Axis limits with padding
    all_fc = df.avg_log2FC
    all_lp = df.neg_log_p
    x_pad = 0.08 * (maximum(all_fc) - minimum(all_fc))
    y_pad = 0.04 * maximum(all_lp)
    xlims = (minimum(all_fc) - x_pad, maximum(all_fc) + x_pad)
    ylims = (0.0, maximum(all_lp) + 4 * y_pad)

    p = plot(
        size=(820, 640),
        title=title,
        xlabel="Log2 Fold Change",
        ylabel="-Log10(adj. p-value)",
        xlims=xlims, ylims=ylims,
        framestyle=:box,
        grid=false,
        legend=:topright,
        foreground_color_legend=nothing,
        background_color_legend=nothing,
        legendfontsize=9
    )

    # NS points — grey, small, semi-transparent
    scatter!(p, ns.avg_log2FC, ns.neg_log_p,
        label="NS",
        color="#aaaaaa",
        markersize=2.5,
        alpha=0.35,
        markerstrokewidth=0)

    # Significant points — deep blue, slightly larger
    scatter!(p, sig.avg_log2FC, sig.neg_log_p,
        label="Significant",
        color="#1565c0",
        markersize=3.5,
        alpha=0.80,
        markerstrokewidth=0)

    # Threshold lines — dashed, subtle
    vline!(p, [logfc_thresh],
        color="#555555", linestyle=:dash, linewidth=1.0, label="")
    hline!(p, [-log10(pval_thresh)],
        color="#555555", linestyle=:dash, linewidth=1.0, label="")

    # Label top significant genes by p-value
    if nrow(sig) > 0
        top = sort(sig, :p_val_adj)[1:min(top_n, nrow(sig)), :]
        for r in eachrow(top)
            annotate!(p,
                r.avg_log2FC,
                r.neg_log_p + 0.6 * y_pad,
                text(r.gene, 7, :bottom, :black))
        end
    end

    return p
end

# ------------------------------------------------------------------------------
# 5. Violin Plot
# ------------------------------------------------------------------------------
function violin_plot(
    obj::SingleCellObject;
    feature::Union{Nothing,String}=nothing,
    features::Union{Nothing,Vector{String}}=nothing,
    group::String="cluster"
)

    if isnothing(feature) && isnothing(features)
        error("Provide either feature=\"GENE\" or features=[\"GENE1\",\"GENE2\"]")
    end

    genes = isnothing(features) ? [feature] : features

    plots = []

    for gene in genes

        gene_idx = findfirst(==(gene), obj.var_data.gene_name)
        isnothing(gene_idx) && error("Gene '$gene' not found.")

        df = DataFrame(
            grp=obj.meta_data[!, group],
            val=vec(obj.norm_data[gene_idx, :])
        )

        sort!(df, :grp)

        p = @df df violin(
            :grp,
            :val,
            group=:grp,
            fillalpha=0.6,
            title=gene,
            ylabel="Expression",
            legend=false
        )

        @df df dotplot!(
            :grp,
            :val,
            color=:black,
            markersize=1,
            alpha=0.2
        )

        push!(plots, p)
    end

    return length(plots) == 1 ?
           plots[1] :
           plot(
        plots...,
        layout=(length(plots), 1),
        size=(900, 350 * length(plots))
    )
end

# ------------------------------------------------------------------------------
# 6.PseudoBulk & PAGA Plot
# ------------------------------------------------------------------------------
function plot_pseudobulk_heatmap(
    pb::DataFrame;
    n_genes::Int=40,
    title::String="Pseudobulk Heatmap",
    ranking::Symbol=:range
)

    gene_names = String.(pb.gene)

    sample_cols = names(pb)[2:end]

    expr = Matrix{Float64}(pb[:, sample_cols])

    # ---------------------------------------------------------
    # Remove genes with no signal
    # ---------------------------------------------------------
    keep_nonzero = vec(sum(expr, dims=2)) .> 0

    expr = expr[keep_nonzero, :]
    gene_names = gene_names[keep_nonzero]

    # ---------------------------------------------------------
    # Gene ranking
    # ---------------------------------------------------------
    scores =
        if ranking == :variance

            vec(var(expr, dims=2))

        elseif ranking == :cv

            vec(std(expr, dims=2)) ./
            (vec(mean(expr, dims=2)) .+ eps())

        elseif ranking == :range

            vec(maximum(expr, dims=2) .-
                minimum(expr, dims=2))

        else

            error(
                "Unknown ranking method '$ranking'. " *
                "Choose :variance, :cv, or :range."
            )

        end

    # ---------------------------------------------------------
    # Select top genes
    # ---------------------------------------------------------
    n_select = min(n_genes, length(scores))

    top_idx = sortperm(scores, rev=true)[1:n_select]

    expr_top = expr[top_idx, :]
    genes_top = gene_names[top_idx]

    # ---------------------------------------------------------
    # Row-wise z-score scaling
    # ---------------------------------------------------------
    expr_scaled = copy(expr_top)

    for i in 1:size(expr_scaled, 1)

        μ = mean(expr_scaled[i, :])
        σ = std(expr_scaled[i, :])

        if σ > 0
            expr_scaled[i, :] .=
                (expr_scaled[i, :] .- μ) ./ σ
        end

    end

    p = heatmap(
        sample_cols,
        reverse(genes_top),
        reverse(expr_scaled, dims=1),
        xlabel="Cell Type",
        ylabel="Genes",
        title=title,
        colorbar=true,
        size=(900, 700)
    )

    return p
end

function paga_plot(
    obj::SingleCellObject;
    reduction::String="umap",
    cluster_col::String="cluster",
    min_conn_weight::Float64=0.02,
    title::String="PAGA Connectivity"
)

    reduction ∉ keys(obj.reductions) &&
        error("Reduction '$reduction' not found.")

    !haskey(obj.graphs, "neighbors") &&
        error("Neighbor graph not found. Run find_neighbors! first.")

    coords = obj.reductions[reduction]
    clusters = obj.meta_data[!, cluster_col]

    u_clusters = sort(unique(clusters))
    n_clusters = length(u_clusters)

    # ---------------------------------------------------------
    # Cluster centroids (median is more robust than mean)
    # ---------------------------------------------------------
    centroids = Dict(
        c => vec(median(coords[clusters.==c, :], dims=1))
        for c in u_clusters
    )

    # ---------------------------------------------------------
    # Cluster sizes
    # ---------------------------------------------------------
    cluster_sizes = [
        sum(clusters .== c)
        for c in u_clusters
    ]

    println("\nCluster sizes:")
    for (c, sz) in zip(u_clusters, cluster_sizes)
        println("Cluster $c => $sz cells")
    end

    # ---------------------------------------------------------
    # Connectivity matrix
    # ---------------------------------------------------------
    connectivity = zeros(Float64, n_clusters, n_clusters)

    cluster_index = Dict(
        c => i for (i, c) in enumerate(u_clusters)
    )

    g = obj.graphs["neighbors"]

    for e in edges(g)

        s = src(e)
        d = dst(e)

        c1 = clusters[s]
        c2 = clusters[d]

        c1 == c2 && continue

        i = cluster_index[c1]
        j = cluster_index[c2]

        connectivity[i, j] += 1
        connectivity[j, i] += 1
    end

    # ---------------------------------------------------------
    # Normalize by cluster size
    # ---------------------------------------------------------
    for i in 1:n_clusters
        for j in 1:n_clusters

            if connectivity[i, j] > 0

                connectivity[i, j] /=
                    sqrt(cluster_sizes[i] * cluster_sizes[j])

            end
        end
    end

    # ---------------------------------------------------------
    # Scale to [0,1]
    # ---------------------------------------------------------
    max_conn = maximum(connectivity)

    if max_conn > 0
        connectivity ./= max_conn
    end

    println("\nConnectivity matrix:")
    display(round.(connectivity, digits=3))

    # ---------------------------------------------------------
    # Plot
    # ---------------------------------------------------------
    p = plot(
        aspect_ratio=1,
        legend=false,
        ticks=false,
        grid=false,
        title=title,
        size=(800, 800)
    )

    # ---------------------------------------------------------
    # Draw edges
    # ---------------------------------------------------------
    for i in 1:n_clusters
        for j in (i+1):n_clusters

            w = connectivity[i, j]

            if w >= min_conn_weight

                p1 = centroids[u_clusters[i]]
                p2 = centroids[u_clusters[j]]

                plot!(
                    p,
                    [p1[1], p2[1]],
                    [p1[2], p2[2]],
                    color=:gray50,
                    alpha=0.7,
                    linewidth=1 + 12 * w,
                    label=""
                )
            end
        end
    end

    # ---------------------------------------------------------
    # Draw nodes
    # ---------------------------------------------------------
    for (i, c) in enumerate(u_clusters)

        xy = centroids[c]

        scatter!(
            p,
            [xy[1]],
            [xy[2]],
            markersize=14,
            color=SCANPY_COLORS[mod1(i, length(SCANPY_COLORS))],
            markerstrokewidth=0.8,
            label=""
        )

        annotate!(
            p,
            xy[1],
            xy[2],
            text(string(c), 9, :black)
        )
    end

    return p
end