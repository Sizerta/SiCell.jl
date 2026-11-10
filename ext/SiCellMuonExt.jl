module SiCellMuonExt

using SiCell
using Muon
using DataFrames
using SparseArrays
using Graphs

import SiCell: from_anndata, to_anndata

# Keep DataFrame annotations and numeric arrays losslessly in SiCell's
# passthrough payload. Nullable numeric arrays may be present in that payload,
# but Muon.AnnData's documented constructor accepts numeric arrays or
# DataFrames, so those arrays are converted to DataFrames on export.
const AnnDataAnnotationArray = AbstractArray{<:Union{Missing,Number}}
const AnnDataAnnotationValue = Union{AnnDataAnnotationArray,DataFrame}
const AnnDataConstructorAnnotationValue = Union{AbstractArray{<:Number},DataFrame}
const AnnDataMatrixValue = AbstractMatrix{<:Number}

function _copy_annotation_mapping(mapping)
    result = Dict{String,AnnDataAnnotationValue}()
    for (key, value) in pairs(mapping)
        if value isa DataFrame || value isa AnnDataAnnotationArray
            result[String(key)] = deepcopy(value)
        end
    end
    return result
end

function _constructor_annotation_value(value::AnnDataAnnotationArray, key::String)
    # Muon.AnnData's documented constructor accepts numeric arrays and
    # DataFrames (unlike some MuData aligned mappings, it does not advertise
    # nullable numeric arrays). Represent nullable 1-D/2-D annotations as a
    # DataFrame so missing values are retained without relying on an
    # unsupported AnnData array type.
    if ndims(value) == 1
        return DataFrame(annotation=collect(value))
    elseif ndims(value) == 2
        return DataFrame(value, :auto)
    end
    throw(ArgumentError(
        "Cannot export nullable annotation '$key' with $(ndims(value)) dimensions to Muon.AnnData. " *
        "Muon.AnnData supports numeric arrays or DataFrames; convert this annotation to a 1-D or 2-D DataFrame first."
    ))
end

function _split_constructor_annotations(mapping)
    accepted = Dict{String,AnnDataConstructorAnnotationValue}()
    for (key, value) in pairs(mapping)
        name = String(key)
        if value isa DataFrame || value isa AbstractArray{<:Number}
            accepted[name] = deepcopy(value)
        elseif value isa AnnDataAnnotationArray
            accepted[name] = _constructor_annotation_value(value, name)
        end
    end
    return accepted
end

function _copy_matrix_mapping(mapping)
    result = Dict{String,AnnDataMatrixValue}()
    for (key, value) in pairs(mapping)
        value isa AbstractMatrix{<:Number} || continue
        result[String(key)] = deepcopy(value)
    end
    return result
end

function _copy_uns(mapping)
    result = Dict{String,Any}()
    for (key, value) in pairs(mapping)
        result[String(key)] = deepcopy(value)
    end
    return result
end

function _to_string_names(names, expected_length, axis_name)
    result = String.(collect(names))
    length(result) == expected_length ||
        throw(DimensionMismatch("AnnData $axis_name names do not match the matrix dimensions."))
    return result
end

"""
    from_anndata(adata::Muon.AnnData; layer=nothing)

Convert a Muon.jl `AnnData` object into a SiCell `SingleCellObject`.

AnnData stores observations × variables, while SiCell stores genes × cells,
so the selected matrix is transposed during import. By default `adata.X` is
used. Pass `layer="counts"` (for example) to use a named AnnData layer instead.

Observation/variable metadata, names, layers, embeddings, pairwise annotations,
and unstructured metadata are retained for export with `to_anndata`.
"""
function from_anndata(adata::Muon.AnnData; layer::Union{Nothing,AbstractString,Symbol}=nothing)
    Muon.isbacked(adata) &&
        throw(ArgumentError("Backed AnnData objects are not supported yet. Read with `Muon.readh5ad(path; backed=false)`."))

    obs_names = _to_string_names(adata.obs_names, size(adata, 1), "observation")
    var_names = _to_string_names(adata.var_names, size(adata, 2), "variable")

    selected = if isnothing(layer)
        adata.X
    else
        layer_name = String(layer)
        haskey(adata.layers, layer_name) ||
            throw(ArgumentError("AnnData does not contain a layer named '$layer_name'. Available layers: $(join(string.(collect(keys(adata.layers))), ", "))."))
        adata.layers[layer_name]
    end

    size(selected) == (length(obs_names), length(var_names)) ||
        throw(DimensionMismatch("Selected AnnData matrix has size $(size(selected)); expected ($(length(obs_names)), $(length(var_names)))."))

    # Muon/AnnData: cells × genes. SiCell: genes × cells.
    counts = sparse(Float64.(transpose(selected)))

    obs = isnothing(adata.obs) ? DataFrame() : deepcopy(adata.obs)
    var = isnothing(adata.var) ? DataFrame() : deepcopy(adata.var)

    # SiCell's QC functions expect a gene_name column. Keep all original var
    # columns and add the canonical name only when it is absent.
    if !("gene_name" in names(var))
        var[!, :gene_name] = copy(var_names)
    end

    obj = SiCell.SingleCellObject(counts, obs, var)

    # Store all observation-aligned embeddings in the passthrough payload and
    # expose 2-D numeric obsm entries through SiCell's reduction API.
    original_obsm = _copy_annotation_mapping(adata.obsm)
    for (key, value) in pairs(original_obsm)
        if value isa AbstractMatrix{<:Number} && size(value, 1) == length(obs_names)
            reduction_key = startswith(key, "X_") ? key[3:end] : key
            obj.reductions[reduction_key] = Matrix{Float64}(value)
        end
    end

    obj.anndata = Dict{String,Any}(
        "obs_names" => obs_names,
        "var_names" => var_names,
        "obsm" => original_obsm,
        "varm" => _copy_annotation_mapping(adata.varm),
        "layers" => _copy_matrix_mapping(adata.layers),
        "obsp" => _copy_matrix_mapping(adata.obsp),
        "varp" => _copy_matrix_mapping(adata.varp),
        "uns" => _copy_uns(adata.uns),
    )

    return obj
end

function _axis_names(obj, axis::Symbol, n::Int)
    key = axis === :obs ? "obs_names" : "var_names"
    payload_names = get(obj.anndata, key, nothing)
    if payload_names isa AbstractVector && length(payload_names) == n
        return String.(payload_names)
    end

    df = axis === :obs ? obj.meta_data : obj.var_data
    fallback = axis === :obs ? "barcode" : "gene_name"
    if fallback in names(df)
        return string.(df[!, fallback])
    end

    prefix = axis === :obs ? "cell" : "gene"
    return ["$(prefix)_$(i)" for i in 1:n]
end

function _add_sicell_graphs!(obsp, obj::SiCell.SingleCellObject, n_cells::Int)
    if haskey(obj.graphs, "neighbors")
        graph = obj.graphs["neighbors"]
        if graph isa Graphs.AbstractGraph && nv(graph) == n_cells
            obsp["connectivities"] = sparse(Float64.(adjacency_matrix(graph)))
        else
            @warn "Skipping SiCell neighbor graph during AnnData export because its type or size is incompatible with the current cell count."
        end
    end

    if haskey(obj.graphs, "neighbor_indices") && haskey(obj.graphs, "neighbor_distances")
        indices = obj.graphs["neighbor_indices"]
        distances = obj.graphs["neighbor_distances"]
        if length(indices) == n_cells && length(distances) == n_cells
            rows = Int[]
            cols = Int[]
            values = Float64[]
            for i in 1:n_cells
                for (j, distance) in zip(indices[i], distances[i])
                    j == i && continue
                    1 <= j <= n_cells || continue
                    isfinite(distance) || continue
                    push!(rows, i)
                    push!(cols, Int(j))
                    push!(values, Float64(distance))
                end
            end
            obsp["distances"] = sparse(rows, cols, values, n_cells, n_cells)
        else
            @warn "Skipping SiCell neighbor distances during AnnData export because their size is incompatible with the current cell count."
        end
    end

    return obsp
end

function _reduction_obsm_key(key::String, original_obsm)
    # Preserve an existing key exactly (e.g. `spatial`), while using the
    # conventional `X_` prefix for reductions created by SiCell.
    haskey(original_obsm, key) && return key
    prefixed = "X_" * key
    haskey(original_obsm, prefixed) && return prefixed
    return prefixed
end

function _selected_output_matrix(obj::SiCell.SingleCellObject, X::Symbol)
    matrix = if X === :counts
        obj.counts
    elseif X === :normalized
        isnothing(obj.norm_data) &&
            throw(ArgumentError("No normalized matrix is available. Run `normalize_data!` first or use X=:counts."))
        obj.norm_data
    elseif X === :scaled
        isnothing(obj.scaled_data) &&
            throw(ArgumentError("No scaled matrix is available. Run `scale_data!` first or use another X option."))
        obj.scaled_data
    else
        throw(ArgumentError("X must be :counts, :normalized, or :scaled; got :$X."))
    end

    size(matrix) == size(obj.counts) ||
        throw(DimensionMismatch("The selected $X matrix has size $(size(matrix)), but AnnData export requires the full gene × cell shape $(size(obj.counts)). Scaled data may contain only highly variable genes."))
    return matrix
end

"""
    to_anndata(obj::SingleCellObject; X=:counts)

Convert a SiCell object into a Muon.jl `AnnData` object.

`X` selects the exported expression matrix: `:counts` (default),
`:normalized`, or `:scaled`. Scaled data are exportable only when they contain
all genes. When `X` is not `:counts`, raw counts are also placed in
`layers["counts"]`.
"""
function to_anndata(obj::SiCell.SingleCellObject; X::Symbol=:counts)
    n_genes, n_cells = size(obj.counts)
    matrix = _selected_output_matrix(obj, X)

    # SiCell: genes × cells. AnnData: cells × genes.
    X_out = matrix isa SparseMatrixCSC ? copy(transpose(matrix)) : permutedims(matrix)

    obs_names = _axis_names(obj, :obs, n_cells)
    var_names = _axis_names(obj, :var, n_genes)
    obs = deepcopy(obj.meta_data)
    var = deepcopy(obj.var_data)

    # Reuse annotations loaded from AnnData, then overwrite embeddings with
    # the current SiCell reductions so newly computed results are exported.
    obsm = get(obj.anndata, "obsm", Dict{String,AnnDataAnnotationValue}())
    obsm = _copy_annotation_mapping(obsm)
    for (key, embedding) in obj.reductions
        size(embedding, 1) == n_cells ||
            throw(DimensionMismatch("Reduction '$key' has $(size(embedding, 1)) rows; expected $n_cells cells."))
        out_key = _reduction_obsm_key(key, obsm)
        obsm[out_key] = Matrix{Float64}(embedding)
    end

    varm = _copy_annotation_mapping(get(obj.anndata, "varm", Dict{String,AnnDataAnnotationValue}()))
    obsm_constructor = _split_constructor_annotations(obsm)
    varm_constructor = _split_constructor_annotations(varm)
    layers = _copy_matrix_mapping(get(obj.anndata, "layers", Dict{String,AnnDataMatrixValue}()))
    obsp = _copy_matrix_mapping(get(obj.anndata, "obsp", Dict{String,AnnDataMatrixValue}()))
    _add_sicell_graphs!(obsp, obj, n_cells)
    varp = _copy_matrix_mapping(get(obj.anndata, "varp", Dict{String,AnnDataMatrixValue}()))
    uns = _copy_uns(get(obj.anndata, "uns", Dict{String,Any}()))

    if X !== :counts
        layers["counts"] = copy(transpose(obj.counts))
    end

    adata_out = Muon.AnnData(
        X=X_out,
        obs=obs,
        obs_names=obs_names,
        var=var,
        var_names=var_names,
        obsm=obsm_constructor,
        varm=varm_constructor,
        layers=layers,
        obsp=obsp,
        varp=varp,
        uns=uns,
    )

    return adata_out
end

end # module SiCellMuonExt
