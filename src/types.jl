# ==============================================================================
# SiCell Core Module
# Defines the SingleCellObject and basic manipulation methods

function _check_dependency(pkg::Symbol, func_name::String)
    # This checks if the package was successfully loaded into the SiCell module
    if !isdefined(SiCell, pkg)
        error("$func_name requires the $pkg package. \n" *
              "Please run: using Pkg; Pkg.add(\"$pkg\") \n" *
              "Then ensure 'using $pkg' is at the top of your SiCell.jl file.")
    end
end
"""
    SingleCellObject

The primary data structure for the SiCell.jl ecosystem. It stores raw counts,
normalized data, cell/gene metadata, and dimensionality reductions.

# Fields
- `counts`: Raw UMI counts (Genes x Cells).
- `norm_data`: Log-normalized or SCTransformed data.
- `scaled_data`: Z-scored data (usually only for highly variable genes).
- `meta_data`: Cell-level information (QC, clusters, pseudotime).
- `var_data`: Gene-level information (Names, IDs, HVG status).
- `reductions`: Dictionary of embeddings (e.g., "pca", "umap", "diffusion").
- `graphs`: Dictionary for adjacency matrices (KNN, SNN).
"""
mutable struct SingleCellObject
    counts::SparseMatrixCSC{Float64,Int64}
    norm_data::Union{SparseMatrixCSC{Float64,Int64},Nothing}
    scaled_data::Union{Matrix{Float64},Nothing}
    meta_data::DataFrame
    var_data::DataFrame
    reductions::Dict{String,Matrix{Float64}}
    graphs::Dict{String,Any}

    function SingleCellObject(raw_counts::AbstractMatrix, meta::DataFrame, var::DataFrame)
        # Ensure counts are sparse and Float64 for mathematical operations
        sparse_counts = isa(raw_counts, SparseMatrixCSC) ?
                        Float64.(raw_counts) :
                        sparse(Float64.(raw_counts))

        # Dimensionality Validation
        n_genes, n_cells = size(sparse_counts)

        if n_cells != nrow(meta)
            error("Dimension mismatch: Matrix has $n_cells columns, but metadata has $(nrow(meta)) rows.")
        end

        if n_genes != nrow(var)
            error("Dimension mismatch: Matrix has $n_genes rows, but var_data has $(nrow(var)) rows.")
        end

        new(sparse_counts, nothing, nothing, meta, var,
            Dict{String,Matrix{Float64}}(), Dict{String,Any}())
    end
end

# ------------------------------------------------------------------------------
# Utilities
# ------------------------------------------------------------------------------

function Base.show(io::IO, obj::SingleCellObject)
    n_genes, n_cells = size(obj.counts)
    println(io, "SingleCellObject with $n_cells cells and $n_genes genes.")
    println(io, "  - Features: ", names(obj.meta_data))
    println(io, "  - Normalized: ", !isnothing(obj.norm_data))
    println(io, "  - Scaled: ", !isnothing(obj.scaled_data))

    if !isempty(obj.reductions)
        println(io, "  - Reductions: ", join(keys(obj.reductions), ", "))
    end

    if !isempty(obj.graphs)
        println(io, "  - Graphs: ", join(keys(obj.graphs), ", "))
    end
end

"""
    subset_object(obj, mask)

Returns a new SingleCellObject containing only the cells selected by the `mask`.
The mask can be a `Vector{Bool}` or a `Vector{Int}`.
"""
function subset_object(obj::SingleCellObject, mask::AbstractVector)
    n_cells = size(obj.counts, 2)

    if eltype(mask) == Bool
        length(mask) != n_cells &&
            error("Boolean mask length must match cell count ($n_cells).")
    else
        (any(mask .< 1) || any(mask .> n_cells)) &&
            error("Cell indices out of bounds. Valid range: 1 to $n_cells.")
    end
    # Create new components
    new_counts = obj.counts[:, mask]
    new_meta = obj.meta_data[mask, :]

    # Initialize the new object (reusing var_data as genes don't change)
    new_obj = SingleCellObject(new_counts, new_meta, copy(obj.var_data))

    # Transfer normalized data if it exists
    if !isnothing(obj.norm_data)
        new_obj.norm_data = obj.norm_data[:, mask]
    end

    # Optional: Reductions are usually invalidated by subsetting unless 
    # the user specifically wants to keep coordinates. 
    # Here we reset them to maintain biological consistency.
    println("Successfully subsetted to $(size(new_obj.counts, 2)) cells.")
    return new_obj
end