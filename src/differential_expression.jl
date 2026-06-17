# ==============================================================================
# SiCell Differential Expression Module
# Contains: Marker Discovery (Wilcoxon Rank-Sum), P-value Adjustment
# ==============================================================================

# Note: This file assumes the following standard libraries are loaded in the parent module:
# Statistics, DataFrames, SparseArrays, MultipleTesting, Distributions, StatsBase, LinearAlgebra

"""
    find_all_markers(obj; group="cluster", min_pct=0.1, logfc_threshold=0.25)

Finds marker genes for each group (cluster) using the Wilcoxon Rank-Sum test (Mann-Whitney U).
This implementation is highly optimized for sparse matrices and uses multi-threading.

# Arguments
- `obj`: SingleCellObject.
- `group`: Metadata column to group by (e.g., "cluster" or "cell_type").
- `min_pct`: Minimum fraction of cells expressing the gene in the cluster.
- `logfc_threshold`: Minimum log2 fold-change.

# Returns
- A `DataFrame` containing genes, p-values, adjusted p-values, and fold-changes.
"""
function find_all_markers(
    obj::SingleCellObject;
    group::String="cluster",
    min_pct::Float64=0.1,
    logfc_threshold::Float64=0.25
)
    # 1. Smart Auto-Detection
    # Fallback to 'graph_cluster' if 'cluster' is missing but graph clustering was performed
    if group == "cluster" && !("cluster" in names(obj.meta_data)) && ("graph_cluster" in names(obj.meta_data))
        group = "graph_cluster"
        println("Notice: 'cluster' column not found. Using 'graph_cluster' for marker discovery.")
    end

    # 2. Pre-flight Checks
    if isnothing(obj.norm_data)
        error("Normalized data not found. Please run `normalize_data!(obj)` first.")
    end

    if !(group in names(obj.meta_data))
        error("Group column '$group' not found in metadata.")
    end

    # 3. Setup Data
    norm_sparse = obj.norm_data
    n_genes, n_cells = size(norm_sparse)
    genes = obj.var_data.gene_name

    clusters = sort(unique(obj.meta_data[!, group]))
    n_clusters = length(clusters)
    cluster_labels = obj.meta_data[!, group]

    cluster_to_index = Dict(c => i for (i, c) in enumerate(clusters))
    cluster_cell_counts = [count(==(c), cluster_labels) for c in clusters]

    # Pre-allocate accumulators (Genes x Clusters)
    SumExpr = zeros(Float64, n_genes, n_clusters)
    SumCounts = zeros(Float64, n_genes, n_clusters)
    SumRanks = zeros(Float64, n_genes, n_clusters)

    # Transpose for faster row-wise (gene-wise) access
    data_T = sparse(norm_sparse')

    println("Calculating statistics for $n_genes genes across $n_clusters clusters using $(Threads.nthreads()) threads...")

    # 4. Single-Pass Statistics Calculation
    Threads.@threads for i in 1:n_genes
        nz_range = nzrange(data_T, i)
        isempty(nz_range) && continue

        gene_vals = view(nonzeros(data_T), nz_range)
        cell_idx = view(rowvals(data_T), nz_range)
        local_labels_for_gene = cluster_labels[cell_idx]

        # A. Sum Expression and Counts (Percentage expressed)
        @inbounds for (val, lbl) in zip(gene_vals, local_labels_for_gene)
            k = cluster_to_index[lbl]
            SumExpr[i, k] += val
            SumCounts[i, k] += 1
        end

        # B. Sum Ranks for Wilcoxon
        # We account for the "zeros" by calculating the average rank of all zero entries
        n_zeros = n_cells - length(gene_vals)
        avg_zero_rank = (1 + n_zeros) / 2
        ranks = tiedrank(gene_vals) .+ n_zeros # Shift ranks by number of zeros

        for k in 1:n_clusters
            sum_nz_ranks = 0.0
            n_nz = 0

            # Manual loop to avoid intermediate allocations
            @inbounds for j in eachindex(local_labels_for_gene)
                if local_labels_for_gene[j] == clusters[k]
                    sum_nz_ranks += ranks[j]
                    n_nz += 1
                end
            end

            n_nz == 0 && continue

            # Rank sum = (Sum of ranks of non-zeros) + (count of zeros in cluster * avg_zero_rank)
            n_zero_in_cluster = cluster_cell_counts[k] - n_nz
            SumRanks[i, k] = sum_nz_ranks + (n_zero_in_cluster * avg_zero_rank)
        end
    end

    # 5. Comparative Analysis (Cluster vs. All Others)
    results = DataFrame()
    total_sum_expr = vec(sum(SumExpr, dims=2))
    total_sum_cnts = vec(sum(SumCounts, dims=2))
    normal_dist = Normal(0, 1)

    for (k, clus) in enumerate(clusters)
        n1 = cluster_cell_counts[k]
        n2 = n_cells - n1
        (n1 == 0 || n2 == 0) && continue

        # Mean Expression and LogFC
        mean1 = SumExpr[:, k] ./ n1
        mean2 = (total_sum_expr .- SumExpr[:, k]) ./ n2
        logfc = log2.(mean1 .+ 1e-9) .- log2.(mean2 .+ 1e-9)

        # Percent Expressed
        pct1 = SumCounts[:, k] ./ n1
        pct2 = (total_sum_cnts .- SumCounts[:, k]) ./ n2

        # Wilcoxon U Statistic & P-value
        # U = R1 - (n1(n1+1)/2)
        U1 = SumRanks[:, k] .- (n1 * (n1 + 1) / 2)
        mu = n1 * n2 / 2
        sigma = sqrt(n1 * n2 * (n_cells + 1) / 12)

        # Standard normal approximation
        Z = [sigma > 0 ? (u - mu) / sigma : 0.0 for u in U1]
        p_vals = 2.0 .* ccdf.(normal_dist, abs.(Z))

        # Filter markers based on thresholds
        mask = (logfc .> logfc_threshold) .& (pct1 .>= min_pct)

        if any(mask)
            append!(results, DataFrame(
                gene=genes[mask],
                cluster=fill(clus, sum(mask)),
                avg_log2FC=logfc[mask],
                p_val=p_vals[mask],
                pct_1=pct1[mask],
                pct_2=pct2[mask]
            ))
        end
    end

    # 6. Post-Processing: Multiple Testing Correction
    if nrow(results) > 0
        println("Adjusting p-values using Benjamini-Hochberg...")
        results.p_val_adj = adjust(results.p_val, BenjaminiHochberg())

        # Final cleanup: Significance filter and sorting
        filter!(row -> row.p_val_adj < 0.05, results)
        sort!(results, [:cluster, :avg_log2FC], rev=[false, true])
    else
        @warn "No marker genes found with current thresholds."
    end

    return results
end