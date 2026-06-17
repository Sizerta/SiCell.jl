# ==============================================================================
# SiCell Batch Integration Module
# Contains: Dataset Merging, Harmony Correction, BBKNN, and Mixing Metrics
# ==============================================================================

# Note: This file assumes the following standard libraries are loaded in the parent module:
# DataFrames, SparseArrays, Statistics, LinearAlgebra, Random, Clustering, NearestNeighbors, Graphs

# ==============================================================================
# 1. Dataset Merging
# ==============================================================================

"""
    merge_batches!(obj1, obj2; batch1_name="Batch1", batch2_name="Batch2")

Merges two `SingleCellObject`s by finding the intersection of genes and concatenating cells.
Returns a **new** `SingleCellObject` containing the combined data.



[Image of Venn diagram set intersection]


# Arguments
- `obj1`, `obj2`: The objects to merge.
- `batch1_name`, `batch2_name`: Labels for the 'batch' metadata column.

# Returns
- A new `SingleCellObject`.
"""
function merge_batches!(obj1::SingleCellObject, obj2::SingleCellObject;
    batch1_name::String="Batch1", batch2_name::String="Batch2")

    println("Merging datasets...")

    # 1. Find common genes (Intersection)
    genes1 = obj1.var_data.gene_name
    genes2 = obj2.var_data.gene_name

    common_genes = intersect(genes1, genes2)
    println("Found $(length(common_genes)) common genes.")

    if length(common_genes) < 100
        error("Very low gene overlap (<100). Check if datasets use the same gene nomenclature.")
    end

    # 2. Get indices for common genes
    # indexin returns indices in the target array (genes1/genes2) corresponding to common_genes
    idx1 = filter(!isnothing, indexin(common_genes, genes1))
    idx2 = filter(!isnothing, indexin(common_genes, genes2))

    # 3. Subset counts matrices to common genes
    counts1_sub = obj1.counts[idx1, :]
    counts2_sub = obj2.counts[idx2, :]

    # 4. Concatenate horizontally (combine cells)
    merged_counts = hcat(counts1_sub, counts2_sub)

    # 5. Prepare Metadata
    meta1 = copy(obj1.meta_data)
    meta2 = copy(obj2.meta_data)

    # Add Batch Identifiers
    meta1[!, :batch] .= batch1_name
    meta2[!, :batch] .= batch2_name

    # Concatenate metadata
    # flexible=true allows for column mismatches (fills missing with missing)
    merged_meta = vcat(meta1, meta2, cols=:union)

    # 6. Prepare Variable Data (Gene info)
    # Use gene info from the first object, subsetted to common genes
    merged_var = obj1.var_data[idx1, :]

    # 7. Create new SingleCellObject
    merged_obj = SingleCellObject(merged_counts, merged_meta, merged_var)

    println("Merge complete. Total cells: $(size(merged_counts, 2))")
    return merged_obj
end

# ==============================================================================
# 2. Harmony Integration (Iterative Correction)
# ==============================================================================

"""
    run_harmony!(obj, batch_key; max_iter=10, n_clusters=10, key="harmony")

Applies the Harmony batch correction algorithm to the PCA reduction.
Aligns cells from different batches in a shared embedding space using iterative clustering.



# Arguments
- `obj`: The `SingleCellObject`.
- `batch_key`: Column name in `obj.meta_data` for batch IDs.
- `key`: Key to store the corrected embedding (default: "harmony").
"""
function run_harmony!(obj::SingleCellObject, batch_key::String;
    max_iter::Int=10, n_clusters::Int=10, key::String="harmony")

    # 1. Validation
    if !("pca" in keys(obj.reductions))
        error("PCA reduction not found. Run `run_pca!` first.")
    end
    if !(batch_key in names(obj.meta_data))
        error("Batch key '$batch_key' not found in metadata.")
    end

    _check_dependency(:Clustering, "Clustering.jl")

    # Dynamic dispatch to Clustering module
    parent_mod = parentmodule(@__MODULE__)
    Clustering = getfield(parent_mod, :Clustering)

    println("Running Harmony Integration (k=$n_clusters, iter=$max_iter)...")

    # Prepare Data
    X = copy(obj.reductions["pca"]) # Cells x Dims

    # Global centering (Standardize)
    X .-= mean(X, dims=1)

    # Map batch IDs to integers
    unique_batches = unique(obj.meta_data[!, batch_key])
    batch_map = Dict(unique_batches .=> 1:length(unique_batches))
    batch_idx = [batch_map[b] for b in obj.meta_data[!, batch_key]]
    n_batches = length(unique_batches)

    # Initialize Z (Embedding)
    Z = copy(X)

    # 2. Harmony Iteration Loop
    for iter in 1:max_iter
        # Step A: Clustering
        # kmeans expects (Dims x Cells), so we transpose Z
        kmeans_result = Clustering.kmeans(transpose(Z), n_clusters; maxiter=100, init=:kmpp, display=:none)
        assignments = kmeans_result.assignments

        # Get global cluster centroids (n_clusters x n_dims)
        cluster_centroids = transpose(kmeans_result.centers)

        # Step B: Evaluate Batch Mixing
        batch_centroids = zeros(n_batches, n_clusters, size(Z, 2))
        batch_counts = zeros(n_batches, n_clusters)

        # Aggregate data by Batch and Cluster
        for i in 1:size(Z, 1)
            c = assignments[i]
            b = batch_idx[i]
            batch_centroids[b, c, :] .+= Z[i, :]
            batch_counts[b, c] += 1
        end

        # Compute averages
        for d in 1:size(Z, 2)
            batch_centroids[:, :, d] ./= max.(batch_counts, 1.0)
        end

        # Step C: Apply Correction
        Z_new = copy(Z)
        for i in 1:size(Z, 1)
            c = assignments[i]
            b = batch_idx[i]
            if batch_counts[b, c] > 0
                # Calculate vector pointing from Batch Centroid to Global Centroid
                correction_vector = batch_centroids[b, c, :] - cluster_centroids[c, :]
                Z_new[i, :] .-= correction_vector
            end
        end

        # Step D: Orthogonalization / Re-scaling
        Z_new .-= mean(Z_new, dims=1)
        norm_initial = norm(Z)
        norm_new = norm(Z_new)
        if norm_new > 0
            Z_new .*= (norm_initial / norm_new)
        end

        Z = Z_new
    end

    obj.reductions[key] = Z
    println("Harmony complete. Corrected embedding stored in reductions[\"$key\"].")
    return nothing
end

# ==============================================================================
# 3. BBKNN (Batch Balanced KNN)
# ==============================================================================

"""
    run_bbknn!(obj, batch_key; reduction="harmony", neighbors_within_batch=3, trim=nothing)

Batch Balanced KNN. Constructs a graph where each cell connects to neighbors 
within its own batch AND neighbors in other batches, ensuring connectivity across datasets.



# Arguments
- `neighbors_within_batch`: How many neighbors to find per batch. 
  (Total neighbors ≈ neighbors_within_batch * n_batches)
- `trim`: Maximum total neighbors to keep per cell.
"""
function run_bbknn!(obj::SingleCellObject, batch_key::String;
    reduction::String="harmony", neighbors_within_batch::Int=3,
    trim::Union{Int,Nothing}=nothing)

    # 1. Validation
    if !(reduction in keys(obj.reductions))
        error("Reduction '$reduction' not found.")
    end
    if !(batch_key in names(obj.meta_data))
        error("Batch key '$batch_key' not found in metadata.")
    end

    println("Running BBKNN graph construction...")
    X = obj.reductions[reduction] # Cells x Dims
    n_cells = size(X, 1)

    unique_batches = unique(obj.meta_data[!, batch_key])
    n_batches = length(unique_batches)
    batch_map = Dict(unique_batches .=> 1:length(unique_batches))
    batch_idx = [batch_map[b] for b in obj.meta_data[!, batch_key]]

    if isnothing(trim)
        trim = 10 * neighbors_within_batch * n_batches
    end

    # Map batch ID -> Vector of cell indices
    batch_dict = Dict(b => findall(==(b), batch_idx) for b in 1:n_batches)

    # Storage for neighbors
    all_neighbors = [Int[] for _ in 1:n_cells]

    # 2. Find Neighbors using KDTree per batch
    data_t = transpose(X) # Dims x Cells (for KDTree)

    # Pre-build trees
    batch_trees = Dict{Int,KDTree}()
    for (b_id, indices) in batch_dict
        if length(indices) > 1
            batch_trees[b_id] = KDTree(data_t[:, indices])
        end
    end

    # Iterate cells
    for i in 1:n_cells
        cell_vec = data_t[:, i]
        temp_neighbors = Tuple{Float64,Int}[]

        # Search in every batch
        for (b_id, indices) in batch_dict
            if length(indices) <= 1
                continue
            end

            tree = batch_trees[b_id]
            k_search = min(neighbors_within_batch + 1, length(indices))

            # Find k neighbors in this specific batch
            idxs, dists = knn(tree, cell_vec, k_search, true)

            for (k, local_idx) in enumerate(idxs)
                global_idx = indices[local_idx]
                if global_idx != i
                    push!(temp_neighbors, (dists[k], global_idx))
                end
            end
        end

        # Sort combined neighbors by distance
        sort!(temp_neighbors)

        # Trim to top K
        keep_count = min(length(temp_neighbors), trim)
        neighbor_indices = unique([n[2] for n in temp_neighbors[1:keep_count]])

        all_neighbors[i] = neighbor_indices
    end

    # 3. Build Graph
    adj_matrix = zeros(Bool, n_cells, n_cells)
    for i in 1:n_cells
        for nbr in all_neighbors[i]
            adj_matrix[i, nbr] = true
            adj_matrix[nbr, i] = true
        end
    end

    g = SimpleGraph(adj_matrix)
    obj.graphs["neighbors"] = g

    # 4. Report Quality
    score = calculate_mixing_score(obj; graph_key="neighbors", batch_key=batch_key)
    perfect_score = 1.0 - (1.0 / n_batches)

    println("BBKNN complete. Graph stored in obj.graphs[\"neighbors\"].")
    println("  🔍 Mixing Score: $(round(score, digits=4)) (Target ≈ $(round(perfect_score, digits=4)))")

    return nothing
end

"""
    calculate_mixing_score(obj; graph_key="neighbors", batch_key="batch")

Calculates a Batch Mixing Score (Entropy-like metric) to evaluate integration quality.
Score ≈ (1 - 1/n_batches) indicates perfect mixing.
"""
function calculate_mixing_score(obj::SingleCellObject; graph_key::String="neighbors", batch_key::String="batch")
    if !(graph_key in keys(obj.graphs))
        error("Graph '$graph_key' not found.")
    end
    if !(batch_key in names(obj.meta_data))
        error("Batch key '$batch_key' not found.")
    end

    g = obj.graphs[graph_key]
    batches = obj.meta_data[!, batch_key]
    n_batches = length(unique(batches))

    if n_batches <= 1
        return 0.0
    end

    total_mixing_ratio = 0.0
    n_valid_cells = 0

    for i in 1:nv(g)
        cell_neighbors = neighbors(g, i)
        if isempty(cell_neighbors)
            continue
        end

        my_batch = batches[i]

        # Count neighbors from DIFFERENT batches
        foreign_neighbors = count(nbr -> batches[nbr] != my_batch, cell_neighbors)

        total_mixing_ratio += foreign_neighbors / length(cell_neighbors)
        n_valid_cells += 1
    end

    return n_valid_cells == 0 ? 0.0 : total_mixing_ratio / n_valid_cells
end
