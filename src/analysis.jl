# ==============================================================================
# SiCell Analysis Module
# Contains: Preprocessing, PCA, Clustering (KMeans/Graph), UMAP, and Pseudobulk
# ==============================================================================

# Note: This file assumes the following standard libraries are loaded in the parent module:
# Statistics, SparseArrays, DataFrames, LinearAlgebra, Random, Arpack, NearestNeighbors, Graphs

# ==============================================================================
# 1. Preprocessing & Normalization
# ==============================================================================
using LinearAlgebra, SparseArrays

"""
    normalize_data!(obj::SingleCellObject; scale_factor::Float64=10000.0)

Log-normalizes count data: log((count / cell_total * scale_factor) + 1).
Stores result in `obj.norm_data`.
"""
function normalize_data!(obj::SingleCellObject; scale_factor::Float64=10000.0)
    # Input validation
    if size(obj.counts, 1) == 0 || size(obj.counts, 2) == 0
        error("Counts matrix is empty. Cannot normalize.")
    end
    if scale_factor <= 0 || !isfinite(scale_factor)
        error("scale_factor must be positive and finite, got $scale_factor")
    end

    println("Normalizing data (scale_factor=$scale_factor)...")

    # Pre-calculate cell totals
    total_counts = vec(sum(obj.counts, dims=1))

    # Copy structure for sparse matrix optimization
    norm_matrix = copy(obj.counts)
    vals = nonzeros(norm_matrix)
    m, n = size(norm_matrix)

    # Threaded normalization
    Threads.@threads for j in 1:n
        cell_total = total_counts[j]
        if cell_total > 0
            # Iterate only over non-zero values for this column (cell)
            for k in nzrange(norm_matrix, j)
                vals[k] = log((vals[k] / cell_total * scale_factor) + 1.0)
            end
        end
    end

    obj.norm_data = norm_matrix
    println("Normalization complete.")
end

"""
    find_variable_features!(obj::SingleCellObject; n_features::Int=2000)

Identifies highly variable genes based on dispersion. 
Updates `obj.var_data` with `:mean`, `:dispersion`, and `:highly_variable` columns.
"""
function find_variable_features!(obj::SingleCellObject; n_features::Int=2000)
    if isnothing(obj.counts) || size(obj.counts, 1) == 0
        error("Counts matrix is empty. Cannot find variable features.")
    end
    if n_features <= 0 || n_features > size(obj.counts, 1)
        error("Invalid n_features: $n_features (Total genes: $(size(obj.counts, 1)))")
    end

    println("Finding variable features...")
    raw = obj.counts

    # Calculate stats
    gene_means = vec(mean(raw, dims=2))
    gene_vars = vec(var(raw, dims=2))

    # Calculate dispersion
    dispersion = gene_vars ./ (gene_means .+ 1e-9)
    dispersion_norm = log.(dispersion)

    # Select top features
    sorted_indices = sortperm(dispersion_norm, rev=true)
    top_indices = sorted_indices[1:n_features]

    # Update Metadata
    obj.var_data[!, :mean] = gene_means
    obj.var_data[!, :dispersion] = dispersion_norm
    obj.var_data[!, :highly_variable] .= false
    obj.var_data[top_indices, :highly_variable] .= true

    println("Found $n_features highly variable genes.")
end

"""
    scale_data!(obj::SingleCellObject; features=nothing, clip_value=10.0)

Performs gene-wise Z-score scaling on normalized data.
Stores result in `obj.scaled_data` (Dense Matrix).
"""
function scale_data!(obj::SingleCellObject; features=nothing, clip_value=10.0)
    if isnothing(obj.norm_data)
        error("Normalized data not found. Run `normalize_data!` first.")
    end

    # Determine Feature Indices
    features_idx::Vector{Int} = if isnothing(features)
        if !hasproperty(obj.var_data, :highly_variable)
            error("Highly variable genes not found. Run `find_variable_features!` first.")
        end
        findall(obj.var_data.highly_variable)
    else
        # Handle user-provided features (String names or Integer indices)
        if eltype(features) <: AbstractString
            if !hasproperty(obj.var_data, :gene_names)
                error("`:gene_names` column not found in `obj.var_data`.")
            end
            indices = [findfirst(==(name), obj.var_data.gene_names) for name in features]
            if any(isnothing, indices)
                missing_names = features[isnothing.(indices)]
                error("Features not found: $(join(missing_names, ", "))")
            end
            collect(Int, indices)
        elseif eltype(features) <: Integer
            collect(Int, features)
        else
            error("`features` must be gene names (String) or indices (Integer).")
        end
    end

    # Extract data (Genes x Cells)
    X = Matrix(obj.norm_data[features_idx, :])

    # Calculate stats (dims=2 operates on rows/genes)
    gene_means = mean(X, dims=2)
    gene_sds = std(X, dims=2) .+ 1e-8 # Epsilon for stability

    # Z-score scaling
    X_scaled = (X .- gene_means) ./ gene_sds

    # Clipping
    if isfinite(clip_value)
        X_scaled = clamp.(X_scaled, -clip_value, clip_value)
    end

    obj.scaled_data = X_scaled
    println("Scaling complete for $(length(features_idx)) features.")
    return obj
end

# ==============================================================================
# 2. Dimensionality Reduction (PCA & UMAP)
# ==============================================================================

"""
    run_pca!(obj; ndims=50, key="pca")

Compute PCA. Uses iterative SVD (Arpack) for sparse data (normalized) 
or standard SVD for dense data (scaled).
"""
function run_pca!(obj::SingleCellObject; ndims::Int=50, key::String="pca")
    n_genes, n_cells = size(obj.counts)
    max_possible = min(n_genes, n_cells) - 1

    if ndims > max_possible
        @warn "Requested ndims ($ndims) exceeds limits. Adjusting to $max_possible."
        ndims = max_possible
    end

    if !isnothing(obj.scaled_data)
        println("Running PCA on scaled (dense) data...")
        # SVD on genes x cells
        res = svd(obj.scaled_data)
        # embedding = Cells x Components
        k = min(ndims, size(res.V, 2))
        embedding = res.V[:, 1:k]
    else
        println("Running PCA on normalized (sparse) data...")
        if !hasproperty(obj.var_data, :highly_variable)
            error("Highly variable genes not found. Run find_variable_features! first.")
        end

        hvg_mask = Vector{Bool}(obj.var_data[!, :highly_variable])
        n_hvg = sum(hvg_mask)

        if ndims >= n_hvg
            @warn "ndims ($ndims) must be less than HVGs ($n_hvg). Adjusting."
            ndims = n_hvg - 1
        end

        embedding = _sparse_pca(obj.norm_data, hvg_mask, ndims)
    end

    obj.reductions[key] = embedding
    println("PCA complete. Stored $(size(embedding, 2)) dimensions in '$key'.")
    return nothing
end

"""
    _sparse_pca(data, hvg_mask, ndims)

Internal helper: Sparse-safe PCA using iterative SVD on transpose.
"""
function _sparse_pca(data::SparseMatrixCSC, hvg_mask::Vector{Bool}, ndims::Int)
    # Subset HVGs (Genes x Cells)
    X = data[hvg_mask, :]

    # Iterative SVD on Transpose (Cells x Genes)
    # Arpack.svds returns tuple (svd_obj, ...)
    svd_tuple = Arpack.svds(X'; nsv=ndims)
    svd_result = svd_tuple[1]

    U = svd_result.U
    S = svd_result.S

    # Embedding = U * Diagonal(S) -> (Cells x Dims)
    return U * Diagonal(S)
end

"""
    run_umap!(obj; dims=320, n_components=2, n_neighbors=15, n_epochs=50, ...)

Runs UMAP on a specific reduction (default: PCA).
"""
function run_umap!(
    obj::SingleCellObject;
    dims::Int=20,
    n_components::Int=2,
    n_neighbors::Int=15,
    n_epochs::Int=20,
    key::String="umap",
    reduction::String="pca"
)
    # Validation
    if !(reduction in keys(obj.reductions))
        error("Reduction '$reduction' not found.")
    end
    if n_components > 3
        error("n_components must be <= 3")
    end

    _check_dependency(:UMAP, "UMAP.jl")

    println("Running UMAP (n_dims=$n_components, k=$n_neighbors, epochs=$n_epochs)...")

    # Get data view
    avail_dims = size(obj.reductions[reduction], 2)
    pca_data = @view obj.reductions[reduction][:, 1:min(dims, avail_dims)]
    data_t = transpose(pca_data) # UMAP expects (Features x Samples)

    # Dynamic dispatch to UMAP package
    parent_mod = parentmodule(@__MODULE__)
    UMAP_mod = getfield(parent_mod, :UMAP)

    umap_embedding = UMAP_mod.umap(
        data_t,
        n_components;
        n_neighbors=n_neighbors,
        n_epochs=n_epochs
    )

    obj.reductions[key] = transpose(umap_embedding)
    println("UMAP complete. Stored in obj.reductions[\"$key\"].")
end

# ==============================================================================
# 3. Clustering (K-Means & Graph-Based)
# ==============================================================================

"""
    find_neighbors!(obj; k=15, dims=nothing, key="pca")

Builds a k-Nearest Neighbors graph from the specified reduction.
Stores the graph in `obj.graphs["neighbors"]`.
"""
function find_neighbors!(
    obj::SingleCellObject;
    k::Int=15,
    dims::Union{Int,Nothing}=nothing,
    key::String="pca"
)

    if !(key in keys(obj.reductions))
        error("Reduction '$key' not found.")
    end

    n_cells, total_dims = size(obj.reductions[key])

    if k >= n_cells
        error("k ($k) must be less than n_cells ($n_cells)")
    end

    # Dimensions to use
    use_dims = isnothing(dims) ? min(30, total_dims) : dims

    println("Building KNN graph (k=$k, dims=$use_dims)...")

    pca_data = @view obj.reductions[key][:, 1:use_dims]
    data_t = transpose(pca_data)  # KDTree expects (Dims × Samples)

    try
        # ------------------------------------------------------------
        # KDTree + KNN search
        # ------------------------------------------------------------
        kdtree = KDTree(data_t)

        # k+1 because each cell finds itself
        idxs, dists = knn(
            kdtree,
            data_t,
            min(k + 1, n_cells),
            true
        )

        # ------------------------------------------------------------
        # Cache neighbors for later reuse
        # (Diffusion Maps, PAGA, trajectory inference, etc.)
        # ------------------------------------------------------------
        obj.graphs["neighbor_indices"] = idxs
        obj.graphs["neighbor_distances"] = dists

        # ------------------------------------------------------------
        # Build graph directly (NO dense adjacency matrix)
        # ------------------------------------------------------------
        g = SimpleGraph(n_cells)

        for i in 1:n_cells
            neighbor_count = 0

            for neighbor_idx in idxs[i]

                if neighbor_idx == i
                    continue
                end

                add_edge!(g, i, neighbor_idx)

                neighbor_count += 1

                if neighbor_count >= k
                    break
                end
            end
        end

        obj.graphs["neighbors"] = g

        println(
            "KNN graph built for $(nv(g)) cells ",
            "with $(ne(g)) edges."
        )

    catch e
        error("Failed to build KNN graph: $e")
    end
end

"""
    run_clustering!(obj; k=10, key="cluster", reduction="pca")

Runs K-Means clustering on the reduction embedding.
"""
function run_clustering!(obj::SingleCellObject; k::Int=10, key::String="cluster", reduction::String="pca")
    if !(reduction in keys(obj.reductions))
        error("Reduction '$reduction' not found.")
    end

    _check_dependency(:Clustering, "Clustering.jl")

    println("Running K-Means clustering (k=$k)...")

    data = transpose(obj.reductions[reduction])
    parent_mod = parentmodule(@__MODULE__)
    Clustering = getfield(parent_mod, :Clustering)

    result = Clustering.kmeans(data, k; maxiter=200, display=:none)

    obj.meta_data[!, key] = assignments(result)
    println("Clustering complete. Found $(length(unique(obj.meta_data[!, key]))) clusters.")
end

"""
    run_graph_clustering!(obj; method="louvain", resolution=1.0, ...)

Runs clustering on the constructed KNN graph. 
Methods: "louvain" (native), "label_propagation", "connected_components".
"""
function run_graph_clustering!(obj::SingleCellObject;
    method::String="louvain",
    resolution::Float64=1.0,
    max_iterations::Int=10,
    key::String="graph_cluster")

    if !haskey(obj.graphs, "neighbors")
        error("KNN graph not found. Run find_neighbors! first.")
    end

    g = obj.graphs["neighbors"]
    if nv(g) == 0
        error("Graph is empty.")
    end

    println("Running graph clustering (method=$method)...")
    clustering = Int[]

    if method == "louvain"
        adj = adjacency_matrix(g)
        clustering = NativeLouvain.run_louvain(adj; resolution=resolution)
    elseif method == "label_propagation"
        clustering = _label_propagation_clustering(g, max_iterations)
    elseif method == "connected_components"
        clustering = _connected_components_clustering(g)
    else
        error("Unknown method '$method'.")
    end

    obj.meta_data[!, key] = clustering
    println("Clustering complete. Found $(length(unique(clustering))) clusters.")
end

# ==============================================================================
# 4. Internal Algorithms & Utilities
# ==============================================================================

"""
Internal Module: NativeLouvain
Implements the Louvain algorithm (Phase 1) for community detection.
"""
module NativeLouvain
using LinearAlgebra
using SparseArrays
using Random

function run_louvain(adj::SparseMatrixCSC{T,Int}; resolution::Float64=1.0) where T
    W = Symmetric(adj)
    n = size(W, 1)
    k = vec(sum(W, dims=1))
    m2 = sum(k) # Total weight (2m)
    communities = collect(1:n)
    comm_k_tot = copy(k)

    improved = true
    cur_iter = 0
    nodes = collect(1:n)

    while improved && cur_iter < 100
        improved = false
        cur_iter += 1
        shuffle!(nodes)

        for i in nodes
            c_old = communities[i]
            ki = k[i]
            neighbor_range = nzrange(adj, i)
            isempty(neighbor_range) && continue

            neighbor_rows = rowvals(adj)
            neighbor_communities = Dict{Int,Float64}()

            for idx in neighbor_range
                neighbor = neighbor_rows[idx]
                w = nonzeros(adj)[idx]
                c_nbr = communities[neighbor]
                neighbor_communities[c_nbr] = get(neighbor_communities, c_nbr, 0.0) + w
            end

            best_c = c_old
            best_gain = 0.0

            for (c_candidate, k_i_in) in neighbor_communities
                c_candidate == c_old && continue
                sigma_tot = comm_k_tot[c_candidate]
                gain = k_i_in - resolution * (ki * sigma_tot / m2)

                if gain > best_gain + 1e-10
                    best_gain = gain
                    best_c = c_candidate
                end
            end

            if best_c != c_old
                comm_k_tot[c_old] -= ki
                comm_k_tot[best_c] += ki
                communities[i] = best_c
                improved = true
            end
        end
    end

    # Renumber 1..K
    unique_comms = sort(unique(communities))
    mapper = Dict(c => i for (i, c) in enumerate(unique_comms))
    return [mapper[c] for c in communities]
end
end # End NativeLouvain

function _label_propagation_clustering(g::AbstractGraph, max_iterations::Int)
    n = nv(g)
    n == 0 && return Int[]
    communities = collect(1:n)

    for _ in 1:max_iterations
        changed = false
        for node in shuffle(1:n)
            nbrs = neighbors(g, node)
            isempty(nbrs) && continue

            # Count neighbor labels
            label_counts = Dict{Int,Int}()
            for nbr in nbrs
                lbl = communities[nbr]
                label_counts[lbl] = get(label_counts, lbl, 0) + 1
            end

            # Find max
            max_c = maximum(values(label_counts))
            candidates = [l for (l, c) in label_counts if c == max_c]
            new_label = candidates[1]

            if communities[node] != new_label
                communities[node] = new_label
                changed = true
            end
        end
        !changed && break
    end

    unique_comms = sort(unique(communities))
    comm_map = Dict(val => idx for (idx, val) in enumerate(unique_comms))
    return [comm_map[c] for c in communities]
end

function _connected_components_clustering(g::AbstractGraph)
    n = nv(g)
    n == 0 && return Int[]
    components = Graphs.connected_components(g)
    assignments = zeros(Int, n)
    for (id, component) in enumerate(components)
        assignments[component] .= id
    end
    return assignments
end

"""
    get_pseudobulk!(obj, group)

Aggregates counts by a metadata group. Returns DataFrame.
"""
function get_pseudobulk!(obj::SingleCellObject, group::String)
    if !(group in names(obj.meta_data))
        error("Group '$group' not found in metadata.")
    end

    groups = obj.meta_data[!, group]
    unique_groups = sort(unique(groups))
    n_genes = size(obj.counts, 1)

    pb_matrix = zeros(Float64, n_genes, length(unique_groups))

    for (i, g) in enumerate(unique_groups)
        indices = findall(==(g), groups)
        pb_matrix[:, i] = vec(sum(obj.counts[:, indices], dims=2))
    end

    df = DataFrame(pb_matrix, Symbol.(unique_groups))
    insertcols!(df, 1, :gene => obj.var_data.gene_name)
    return df
end
