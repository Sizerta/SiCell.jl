# ==============================================================================
# SiCell Trajectory Inference Module
# Methods: Diffusion Maps, Graph-based Pseudotime,
#          Trajectory Uncertainty Framework (TES + TDS + LPS)
# ==============================================================================

"""
    run_diffusion_map!(obj; n_components=15, k=20, alpha=0.5)

Computes a Diffusion Map embedding.
"""
function run_diffusion_map!(
    obj::SingleCellObject;
    n_components::Int=15,
    k::Int=20,
    reduction::String="pca",
    key::String="diffusion",
    alpha::Float64=0.5
)
    !(reduction in keys(obj.reductions)) &&
        error("PCA required for Diffusion Maps.")

    X = obj.reductions[reduction]
    n_cells = size(X, 1)

    println("Constructing Diffusion Map from $n_cells cells...")

    if haskey(obj.graphs, "neighbor_indices")
        println("Using existing neighbor graph...")
        idxs = obj.graphs["neighbor_indices"]
        dists = obj.graphs["neighbor_distances"]
    else
        println("Building KNN search...")
        tree = KDTree(transpose(X))
        idxs, dists = knn(tree, transpose(X), min(k + 1, n_cells), true)
        obj.graphs["neighbor_indices"] = idxs
        obj.graphs["neighbor_distances"] = dists
    end

    sigmas = [d[min(k + 1, length(d))] for d in dists]

    I, J, V = Int[], Int[], Float64[]
    sizehint!(I, n_cells * k)
    sizehint!(J, n_cells * k)
    sizehint!(V, n_cells * k)

    for i in 1:n_cells
        for (pos, nbr) in enumerate(idxs[i])
            nbr == i && continue
            dist_sq = dists[i][pos]^2
            s_ij = sigmas[i] * sigmas[nbr]
            if s_ij > 1e-12
                push!(I, i)
                push!(J, nbr)
                push!(V, exp(-dist_sq / s_ij))
            end
        end
    end

    W = sparse(I, J, V, n_cells, n_cells)
    W = max.(W, W')

    d_raw = vec(sum(W, dims=2))
    scale = d_raw .^ (-alpha)
    scale[.!isfinite.(scale)] .= 0.0
    D_alpha = sparse(Diagonal(scale))
    W_alpha = D_alpha * W * D_alpha

    d_alpha = vec(sum(W_alpha, dims=2))
    inv_sqrt_deg = 1.0 ./ sqrt.(d_alpha)
    inv_sqrt_deg[.!isfinite.(inv_sqrt_deg)] .= 0.0
    D_inv_sqrt = sparse(Diagonal(inv_sqrt_deg))
    L_sym = D_inv_sqrt * W_alpha * D_inv_sqrt

    vals, vecs = eigs(L_sym; nev=n_components + 1, which=:LM)
    vals, vecs = real.(vals), real.(vecs)
    λ = vals[2:end]
    embedding = (D_inv_sqrt*vecs)[:, 2:end]

    for i in 1:length(λ)
        denom = 1 - λ[i]
        if abs(denom) > 1e-12
            embedding[:, i] .*= λ[i] / denom
        end
    end

    obj.reductions[key] = embedding
    println("Diffusion Map complete. Components: $n_components")
    return nothing
end

# ==============================================================================
# Graph-based Pseudotime 
# ==============================================================================

"""
    run_pseudotime!(obj, root_cell_idx; method="graph", key="pseudotime")

Compute pseudotime from a root cell index.

Methods:
  "graph"      — Dijkstra shortest paths on KNN graph (geodesic, recommended)
  "euclidean"  — Euclidean distance from root cell in reduction space
  "dc1"        — First diffusion component normalized to [0,1] (continuous,
                 no root cell needed but root_cell_idx used to set direction)


"""
function run_pseudotime!(obj::SingleCellObject, root_cell_idx::Int;
    reduction::String="diffusion",
    key::String="pseudotime",
    method::String="dc1")

    if method == "graph"
        !haskey(obj.graphs, "neighbors") && error("KNN graph required for graph pseudotime. Run find_neighbors!() first.")
        g = obj.graphs["neighbors"]
        println("Computing graph-based pseudotime (Dijkstra) from cell $root_cell_idx...")
        ds = dijkstra_shortest_paths(g, root_cell_idx)
        dists = ds.dists
        max_finite = maximum(filter(isfinite, dists))
        dists[.!isfinite.(dists)] .= max_finite
        obj.meta_data[!, key] = dists ./ maximum(dists)

    elseif method == "euclidean"
        !(reduction in keys(obj.reductions)) && error("Reduction '$reduction' not found. Run run_diffusion_map!() first.")
        X = obj.reductions[reduction]
        root_vec = X[root_cell_idx, :]
        dists = [sqrt(sum((X[i, :] .- root_vec) .^ 2)) for i in 1:size(X, 1)]
        obj.meta_data[!, key] = dists ./ maximum(dists)

    elseif method == "dc1"
        # Continuous pseudotime from first diffusion component
        # DC1 is the smoothest function on the graph manifold and
        # captures the dominant axis of transcriptional progression
        !(reduction in keys(obj.reductions)) && error("Reduction '$reduction' not found. Run run_diffusion_map!() first.")
        X = obj.reductions[reduction]
        size(X, 2) < 1 && error("Diffusion reduction has no components.")

        dc1 = X[:, 1]

        # Normalize to [0, 1]
        dc1_min = minimum(dc1)
        dc1_max = maximum(dc1)
        dc1_max ≈ dc1_min && error("DC1 has zero variance — diffusion map may not have converged.")
        pt = (dc1 .- dc1_min) ./ (dc1_max - dc1_min)

        # Orient so root cell has low pseudotime (flip if needed)
        if pt[root_cell_idx] > 0.5
            pt = 1.0 .- pt
            println("  DC1 flipped: root cell was at $(round(1-pt[root_cell_idx]+pt[root_cell_idx], digits=3)), now at $(round(pt[root_cell_idx], digits=3))")
        end

        println("Computing DC1-based continuous pseudotime...")
        println("  Root cell $(root_cell_idx) pseudotime: $(round(pt[root_cell_idx], digits=4))")
        println("  Unique values: $(length(unique(pt)))  (continuous = good)")
        println("  Range: $(round(minimum(pt), digits=4)) – $(round(maximum(pt), digits=4))")

        obj.meta_data[!, key] = pt

    else
        error("Unknown method: '$method'. Use 'graph', 'euclidean', or 'dc1'.")
    end

    println("Pseudotime stored in meta_data['$key'] (method=$method).")
end

"""
    run_pseudotime!(obj, root_group; cluster_col, key, method)

Compute pseudotime using the first cell of root_group as the root.

"""
function run_pseudotime!(
    obj::SingleCellObject, root_group::String;
    cluster_col::String="cell_type",
    key::String="pseudotime",
    method::String="graph",
    reduction::String="diffusion"
)
    !(cluster_col in names(obj.meta_data)) && error("Column '$cluster_col' not found.")
    idx = findfirst(obj.meta_data[!, cluster_col] .== root_group)
    isnothing(idx) && error("Root group '$root_group' not found in column '$cluster_col'.")
    run_pseudotime!(obj, idx; key=key, method=method, reduction=reduction)
end
"""
    pseudotime_uncertainty!(obj; ...)

**Deprecated**: Please use `trajectory_uncertainty!(obj)` instead.

This function is kept only for compatibility with existing scripts.
"""
function pseudotime_uncertainty!(obj::SingleCellObject;
    pseudotime_key::String="pseudotime",
    key::String="pseudotime_uncertainty")

    @warn "pseudotime_uncertainty! is deprecated. Use trajectory_uncertainty! instead."

    !(pseudotime_key in names(obj.meta_data)) &&
        error("'$pseudotime_key' not found.")

    !haskey(obj.graphs, "neighbors") &&
        error("Neighbor graph not found. Run find_neighbors! first.")

    pt = obj.meta_data[!, pseudotime_key]
    g = obj.graphs["neighbors"]
    n_cells = length(pt)

    scores = zeros(Float64, n_cells)
    for i in 1:n_cells
        nbrs = neighbors(g, i)
        isempty(nbrs) && continue
        scores[i] = std(pt[nbrs])
    end

    max_score = maximum(scores)
    scores = max_score > 0 ? scores ./ max_score : scores

    obj.meta_data[!, key] = scores
    println("Legacy pseudotime uncertainty stored in meta_data['$key'].")
end
# ==============================================================================
# Trajectory Uncertainty Framework (TES + TDS + LPS)
# ==============================================================================

"""
    trajectory_uncertainty!(obj; n_bins=10, key_prefix="traj_unc")

Computes three complementary metrics:

- **TES**: Temporal Entropy (mixing of developmental stages)
- **TDS**: Trajectory Divergence (branching / directional ambiguity)
- **LPS**: Local Progression Score (raw: positive = neighbors ahead)

** recommended**: Use graph-based pseudotime.
"""

function trajectory_uncertainty!(
    obj::SingleCellObject;
    pseudotime_key::String="pseudotime",
    reduction::String="diffusion",
    n_bins::Int=10,
    key_prefix::String="traj_unc"
)
    !(pseudotime_key in names(obj.meta_data)) && error("Run run_pseudotime! first.")
    !haskey(obj.graphs, "neighbors") && error("Run find_neighbors! first.")
    !(reduction in keys(obj.reductions)) && error("Reduction '$reduction' not found.")

    pt = obj.meta_data[!, pseudotime_key]

    # Bounds check — TES is only bounded [0,1] when pseudotime is in [0,1]
    minimum(pt) < 0.0 && error("Pseudotime contains negative values. Pseudotime must be in [0,1].")
    maximum(pt) > 1.0 && error("Pseudotime contains values > 1. Pseudotime must be in [0,1].")

    g = obj.graphs["neighbors"]
    X = obj.reductions[reduction]
    n_cells = length(pt)

    tes = zeros(Float64, n_cells)
    tds = zeros(Float64, n_cells)
    lps_raw = zeros(Float64, n_cells)

    println("Computing TES + Forward-TDS + LPS (n_bins=$n_bins)...")

    for i in 1:n_cells
        nbrs = neighbors(g, i)
        length(nbrs) < 5 && continue

        nbr_pt = pt[nbrs]
        self_pt = pt[i]

        # ── TES — Pairwise Temporal Disagreement, O(k log k) ──────────────
        # Uses the sorted identity:
        #   sum_{i<j} (x_j - x_i) = sum_j (2j - k - 1) * x_j   (1-indexed)
        # Normalized by number of pairs: k*(k-1)/2
        # Result is bounded [0,1] when pseudotime is in [0,1]
        k = length(nbr_pt)
        if k >= 2
            sorted_pt = sort(nbr_pt)
            weighted_sum = 0.0
            for j in 1:k
                weighted_sum += (2j - k - 1) * sorted_pt[j]
            end
            tes[i] = weighted_sum / (k * (k - 1) / 2)
        else
            tes[i] = 0.0
        end

        # ── LPS — Local Progression (unchanged) ───────────────────────────
        lps_raw[i] = mean(nbr_pt) - self_pt

        # ── TDS — Forward Directional Divergence (unchanged) ──────────────
        forward_mask = nbr_pt .> self_pt
        if sum(forward_mask) >= 3 && size(X, 2) >= 2
            forward_nbrs = nbrs[forward_mask]
            center = X[i, :]
            vecs = X[forward_nbrs, :] .- center'
            norms = sqrt.(sum(vecs .^ 2, dims=2))
            norms[norms.<1e-12] .= 1.0
            vecs ./= norms
            mean_vec = vec(mean(vecs, dims=1))
            tds[i] = 1.0 - norm(mean_vec)
        else
            tds[i] = 0.0
        end
    end

    obj.meta_data[!, "$(key_prefix)_tes"] = tes
    obj.meta_data[!, "$(key_prefix)_tds"] = tds
    obj.meta_data[!, "$(key_prefix)_lps"] = lps_raw

    println("trajectory uncertainty metrics stored (prefix: $key_prefix)")
    println("   → TES (temporal mixing)")
    println("   → TDS (forward directional divergence / branching)")
    println("   → LPS (local progression — positive = advancing)")
end
# ==============================================================================
# Visualization
# ==============================================================================

"""
    trajectory_uncertainty_plot(obj; reduction="umap", key_prefix="traj_unc")
"""
function trajectory_uncertainty_plot(obj::SingleCellObject;
    reduction::String="umap",
    key_prefix::String="traj_unc",
    title::String="Trajectory Uncertainty Framework"
)
    !(reduction in keys(obj.reductions)) && error("Reduction not found.")

    tes_key = "$(key_prefix)_tes"
    tds_key = "$(key_prefix)_tds"
    lps_key = "$(key_prefix)_lps"

    if !(tes_key in names(obj.meta_data))
        error("Run trajectory_uncertainty! first.")
    end

    coords = obj.reductions[reduction]
    x, y = coords[:, 1], coords[:, 2]
    n = length(x)
    style = _embedding_style(n)

    # Slightly larger markers + higher alpha for better visibility on large datasets
    ms = style.ms * 1.15
    α = min(0.92, style.α * 1.3)

    p1 = scatter(x, y, zcolor=obj.meta_data[!, tes_key],
        title="TES (Temporal Mixing)", colormap=:viridis,
        colorbar_title="TES", markersize=ms, alpha=α, markerstrokewidth=0)

    p2 = scatter(x, y, zcolor=obj.meta_data[!, tds_key],
        title="TDS (Branching / Divergence)", colormap=:plasma,
        colorbar_title="TDS", markersize=ms, alpha=α, markerstrokewidth=0)

    p3 = scatter(x, y, zcolor=obj.meta_data[!, lps_key],
        title="LPS (Local Progression)", colormap=:RdBu,
        colorbar_title="LPS (raw)", markersize=ms, alpha=α, markerstrokewidth=0,
        clims=(-maximum(abs.(obj.meta_data[!, lps_key])), maximum(abs.(obj.meta_data[!, lps_key]))))

    return plot(p1, p2, p3, layout=(1, 3), size=(1850, 580), title=title)
end

"""
    uncertainty_plot(...) — Enhanced TCS Plot with top 10% outline
"""
function uncertainty_plot(obj::SingleCellObject;
    reduction::String="umap",
    pseudotime_key::String="pseudotime",
    uncertainty_key::String="pseudotime_uncertainty",
    title::String="Trajectory Confidence Score (TCS)"
)
    !(reduction in keys(obj.reductions)) && error("Reduction '$reduction' not found.")
    !(pseudotime_key in names(obj.meta_data)) && error("'$pseudotime_key' not found.")
    !(uncertainty_key in names(obj.meta_data)) && error("'$uncertainty_key' not found.")

    coords = obj.reductions[reduction]
    x, y = coords[:, 1], coords[:, 2]
    style = _embedding_style(length(x))

    pt = obj.meta_data[!, pseudotime_key]
    unc = obj.meta_data[!, uncertainty_key]

    top10_threshold = quantile(unc, 0.9)
    highlight = unc .>= top10_threshold

    xlims = _clip_limits(x)
    ylims = _clip_limits(y)

    p = plot(size=(850, 660), title=title,
        xlabel=uppercase(reduction) * " 1", ylabel=uppercase(reduction) * " 2",
        colorbar_title="TCS", legend=false, xlims=xlims, ylims=ylims, dpi=300)

    # Background
    scatter!(p, x, y, zcolor=pt, alpha=0.65, markersize=style.ms * 0.95,
        markerstrokewidth=0, colormap=:viridis, clims=(0.0, 1.0))

    # Top 10% highlighted
    scatter!(p, x[highlight], y[highlight], zcolor=pt[highlight],
        alpha=0.95, markersize=style.ms * 1.5, markerstrokewidth=1.4,
        markerstrokecolor=:black, colormap=:viridis, clims=(0.0, 1.0))

    return p
end

"""
    highlight_top_uncertainty(obj; top_pct=0.05)
"""
function highlight_top_uncertainty(obj::SingleCellObject;
    reduction::String="umap",
    key_prefix::String="traj_unc",
    top_pct::Float64=0.05,
    title::String="Top 5% Uncertainty Cells"
)
    coords = obj.reductions[reduction]
    x, y = coords[:, 1], coords[:, 2]
    style = _embedding_style(length(x))

    tes_key = "$(key_prefix)_tes"
    tds_key = "$(key_prefix)_tds"

    top_tes = obj.meta_data[!, tes_key] .>= quantile(obj.meta_data[!, tes_key], 1 - top_pct)
    top_tds = obj.meta_data[!, tds_key] .>= quantile(obj.meta_data[!, tds_key], 1 - top_pct)

    p1 = scatter(x, y, color="#eeeeee", alpha=0.45, markersize=style.ms * 0.9, markerstrokewidth=0, label="All cells")
    scatter!(p1, x[top_tes], y[top_tes], color="#d32f2f", alpha=0.95, markersize=style.ms * 2.0,
        markerstrokewidth=1.0, markerstrokecolor=:black, label="Top $(Int(top_pct*100))% TES")

    p2 = scatter(x, y, color="#eeeeee", alpha=0.45, markersize=style.ms * 0.9, markerstrokewidth=0, label="All cells")
    scatter!(p2, x[top_tds], y[top_tds], color="#1976d2", alpha=0.95, markersize=style.ms * 2.0,
        markerstrokewidth=1.0, markerstrokecolor=:black, label="Top $(Int(top_pct*100))% TDS")

    return plot(p1, p2, layout=(1, 2), size=(1450, 650), title=title, legend=:topright)
end