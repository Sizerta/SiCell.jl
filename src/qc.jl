# ==============================================================================
# SiCell Quality Control Module
# Contains: QC Metric Calculation, Cell Filtering
# ==============================================================================

"""
    calculate_qc_metrics!(obj::SingleCellObject; mito_prefix::String="MT-")

Calculates standard single-cell QC metrics and updates the object's metadata:
1. `n_counts`: Total UMI counts per cell (Library Size).
2. `n_genes`: Number of unique genes detected per cell.
3. `percent_mito`: Percentage of reads mapping to mitochondrial genes.

# Arguments
- `obj`: SingleCellObject.
- `mito_prefix`: Prefix used to identify mitochondrial genes (e.g., "MT-" for human, "mt-" for mouse).
"""
function calculate_qc_metrics!(obj::SingleCellObject; mito_prefix::String="MT-")
    # 1. Validation
    n_features, n_cells = size(obj.counts)
    if n_cells == 0
        error("Counts matrix is empty. Cannot calculate QC metrics.")
    end

    # 2. Total Counts (Column Sums)
    # Using vec(sum(...)) on a SparseMatrix is efficient
    n_counts = vec(sum(obj.counts, dims=1))

    # 3. Detected Genes
    # Optimized: Instead of creating a full binary matrix with .> 0, 
    # we count non-zeros in each column of the CSC matrix
    n_genes = zeros(Int, n_cells)
    colptr = obj.counts.colptr
    for j in 1:n_cells
        n_genes[j] = colptr[j+1] - colptr[j]
    end

    # 4. Mitochondrial Content
    gene_names = obj.var_data.gene_name
    mito_indices = findall(x -> startswith(lowercase(x), lowercase(mito_prefix)), gene_names)

    if isempty(mito_indices)
        @warn "No mitochondrial genes found with prefix '$mito_prefix'. Check your gene naming convention."
        percent_mito = zeros(Float64, n_cells)
    else
        # Sum counts for mito-genes only
        mito_counts = vec(sum(obj.counts[mito_indices, :], dims=1))
        percent_mito = (mito_counts ./ n_counts) .* 100.0
    end

    # 5. Update Metadata
    obj.meta_data[!, :n_counts] = n_counts
    obj.meta_data[!, :n_genes] = n_genes
    obj.meta_data[!, :percent_mito] = percent_mito

    println("QC Metrics calculated: avg $(round(mean(n_genes), digits=0)) genes/cell.")
    return nothing
end

"""
    filter_cells!(obj::SingleCellObject; 
                 min_genes=200, max_genes=Inf, 
                 min_counts=0, max_mito=5.0)

Filters cells based on calculated QC metrics. Updates the object in-place.
"""
function filter_cells!(obj::SingleCellObject;
    min_genes::Number=200,
    max_genes::Number=Inf,
    min_counts::Number=0,
    max_mito::Number=5.0
)
    # Check if metrics exist
    if !all(col -> col in names(obj.meta_data), ["n_genes", "percent_mito", "n_counts"])
        error("QC metrics not found in metadata. Run `calculate_qc_metrics!(obj)` first.")
    end

    # Create filter mask
    keep_mask = (obj.meta_data.n_genes .>= min_genes) .&
                (obj.meta_data.n_genes .<= max_genes) .&
                (obj.meta_data.n_counts .>= min_counts) .&
                (obj.meta_data.percent_mito .<= max_mito)

    n_before = size(obj.counts, 2)
    n_after = sum(keep_mask)

    if n_after == 0
        @error "No cells passed the filtering criteria! Object was not modified."
        return nothing
    end

    # Subset Data
    obj.counts = obj.counts[:, keep_mask]
    obj.meta_data = obj.meta_data[keep_mask, :]

    # Subset dependent data if present
    if !isnothing(obj.norm_data)
        obj.norm_data = obj.norm_data[:, keep_mask]
    end

    # Subset dimensionality reductions (UMAP, PCA, etc.)
    for key in keys(obj.reductions)
        # Note: Reductions are stored as (Dims x Cells) internally
        obj.reductions[key] = obj.reductions[key][:, keep_mask]
    end

    println("Filtering complete: Kept $n_after / $n_before cells ($(round(n_after/n_before*100, digits=1))%).")
    return nothing
end