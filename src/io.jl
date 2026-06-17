# =============================================================================
# io.jl - Data Input/Output for SiCell.jl
# =============================================================================

using MatrixMarket
using CSV
using DataFrames
using GZip
using HDF5
using SparseArrays

# =============================================================================
# 10x Genomics Format (Directory)
# =============================================================================

"""
    read_10x(data_dir::String)

Reads a standard 10x Genomics directory containing matrix.mtx, barcodes.tsv, and features.tsv (or genes.tsv).
Returns a `SingleCellObject`.
"""
function read_10x(data_dir::String)
    # Input validation
    if !isdir(data_dir)
        error("Directory '$data_dir' does not exist or is not a directory.")
    end
    if isempty(data_dir)
        error("Data directory path cannot be empty.")
    end

    # Define expected file names (handle .gz versions too)
    mtx_path = joinpath(data_dir, "matrix.mtx")
    barcodes_path = joinpath(data_dir, "barcodes.tsv")
    features_path = joinpath(data_dir, "features.tsv")

    # Check for gzipped versions if plain files don't exist
    if !isfile(mtx_path) && isfile(mtx_path * ".gz")
        mtx_path *= ".gz"
    end
    if !isfile(barcodes_path) && isfile(barcodes_path * ".gz")
        barcodes_path *= ".gz"
    end
    if !isfile(features_path) && isfile(features_path * ".gz")
        features_path *= ".gz"
    elseif !isfile(features_path) && isfile(joinpath(data_dir, "genes.tsv"))
        features_path = joinpath(data_dir, "genes.tsv")
    end

    println("Loading 10x data from: $data_dir")

    # Validate files exist
    if !isfile(mtx_path) && !isfile(mtx_path * ".gz")
        error("Matrix file not found. Expected 'matrix.mtx' or 'matrix.mtx.gz'")
    end
    if !isfile(barcodes_path) && !isfile(barcodes_path * ".gz")
        error("Barcodes file not found.")
    end
    if !isfile(features_path) && !isfile(features_path * ".gz") && !isfile(joinpath(data_dir, "genes.tsv"))
        error("Features file not found.")
    end

    # 1. Load Matrix
    raw_counts = MatrixMarket.mmread(mtx_path)
    if size(raw_counts, 1) == 0 || size(raw_counts, 2) == 0
        error("Loaded matrix is empty. Check matrix.mtx file.")
    end
    if any(x -> x < 0, nonzeros(raw_counts))
        error("Matrix contains negative values. Counts should be non-negative.")
    end

    # 2. Load Barcodes (Cells)
    barcodes_df = CSV.read(barcodes_path, DataFrame; header=false)
    if nrow(barcodes_df) == 0
        error("Barcodes file is empty.")
    end
    if nrow(barcodes_df) != size(raw_counts, 2)
        error("Number of barcodes ($(nrow(barcodes_df))) does not match matrix columns ($(size(raw_counts, 2))).")
    end
    rename!(barcodes_df, 1 => "barcode")

    # 3. Load Features (Genes)
    features_df = CSV.read(features_path, DataFrame; header=false)
    if nrow(features_df) == 0
        error("Features file is empty.")
    end
    if nrow(features_df) != size(raw_counts, 1)
        error("Number of features ($(nrow(features_df))) does not match matrix rows ($(size(raw_counts, 1))).")
    end
    if ncol(features_df) < 2
        rename!(features_df, 1 => "gene_name")
        features_df.gene_id = features_df.gene_name
    else
        rename!(features_df, 1 => "gene_id", 2 => "gene_name")
    end

    return SingleCellObject(raw_counts, barcodes_df, features_df)
end


# =============================================================================
# AnnData (.h5ad) Format
# =============================================================================

function _read_h5ad_dataframe(file::HDF5.File, group_name::String)

    if !haskey(file, group_name)
        return DataFrame()
    end

    g = file[group_name]
    df = DataFrame()

    # --------------------------------------------------
    # Index column
    # --------------------------------------------------

    if haskey(g, "_index")
        df[!, "_index"] = vec(read(g["_index"]))
    elseif haskey(g, "index")
        df[!, "_index"] = vec(read(g["index"]))
    end

    # --------------------------------------------------
    # Old Scanpy categories
    # obs/__categories/*
    # --------------------------------------------------

    categories_group = nothing

    if haskey(g, "__categories")
        categories_group = g["__categories"]
    elseif haskey(g, "__categories__")
        categories_group = g["__categories__"]
    end

    # --------------------------------------------------
    # Read all columns
    # --------------------------------------------------

    for key in keys(g)

        if key in ["_index", "index", "__categories", "__categories__"]
            continue
        end

        ds = g[key]

        if !(ds isa HDF5.Dataset)
            continue
        end

        try

            data = read(ds)

            # Flatten vectors
            if isa(data, Matrix) &&
               (size(data, 1) == 1 || size(data, 2) == 1)

                data = vec(data)
            end

            # ==================================================
            # FORMAT 1:
            # old Scanpy categorical
            # ==================================================

            if !isnothing(categories_group) &&
               haskey(categories_group, key)

                labels = read(categories_group[key])

                mapped = [
                    x >= 0 ?
                    String(labels[Int(x)+1]) :
                    missing
                    for x in data
                ]

                df[!, key] = mapped
                continue
            end

            # ==================================================
            # FORMAT 2:
            # AnnData categorical references
            # ==================================================

            attrs = attributes(ds)

            if haskey(attrs, "categories")

                cat_ref = read(attrs["categories"])

                cat_labels = nothing

                try

                    if isa(cat_ref, HDF5.Reference)

                        ref_obj = HDF5.open(file, cat_ref)
                        cat_labels = read(ref_obj)
                        close(ref_obj)

                    elseif isa(cat_ref, AbstractString)

                        cat_labels = read(file[String(cat_ref)])

                    end

                catch e
                    @warn "Could not resolve categories for $key" exception = e
                end

                if !isnothing(cat_labels)

                    mapped = [
                        x >= 0 ?
                        String(cat_labels[Int(x)+1]) :
                        missing
                        for x in data
                    ]

                    df[!, key] = mapped

                else
                    df[!, key] = data
                end

            else

                df[!, key] = data

            end

        catch e

            @warn "Skipping column $key" exception = e

        end

    end

    return df
end

"""
    read_h5ad(filepath::String; raw::Bool=true)

Reads an AnnData (.h5ad) file - a standard format for single-cell RNA-seq data.
"""
function read_h5ad(filepath::String; raw::Bool=true)
    # Input Validation
    if !isfile(filepath)
        error("File '$filepath' does not exist.")
    end
    if !endswith(filepath, ".h5ad")
        println("Warning: File '$filepath' does not have a .h5ad extension.")
    end

    println("Loading .h5ad from: $filepath")

    h5open(filepath, "r") do file
        # 1. Load Counts Matrix
        group_path = "X"
        if raw && haskey(file, "raw") && haskey(file["raw"], "X")
            group_path = "raw/X"
            println("  Using raw.X for counts...")
        end

        counts_obj = file[group_path]
        counts = nothing

        if isa(counts_obj, HDF5.Group)
            # Sparse Matrix Handling (CSR or CSC)
            attrs = attributes(counts_obj)
            encoding = haskey(attrs, "encoding-type") ? read(attrs["encoding-type"]) : "csr_matrix"
            shape = read(attrs["shape"]) # (N_cells, M_genes)

            data = read(counts_obj["data"])
            indices = read(counts_obj["indices"]) .+ 1 # 0-based to 1-based
            indptr = read(counts_obj["indptr"]) .+ 1   # 0-based to 1-based

            if encoding == "csr_matrix"
                n_cells = Int(shape[1])
                n_genes = Int(shape[2])
                counts = SparseMatrixCSC(n_genes, n_cells, indptr, indices, data)
            elseif encoding == "csc_matrix"
                n_cells = Int(shape[1])
                n_genes = Int(shape[2])
                mat_t = SparseMatrixCSC(n_cells, n_genes, indptr, indices, data)
                counts = sparse(transpose(mat_t))
            else
                error("Unsupported sparse encoding: $encoding")
            end
        elseif isa(counts_obj, HDF5.Dataset)
            # Dense Matrix
            matrix = read(counts_obj)
            counts = sparse(matrix)
        else
            error("Could not read matrix from $group_path")
        end

        # 2. Load Metadata
        obs_df = _read_h5ad_dataframe(file, "obs")
        var_group = (group_path == "raw/X") ? "raw/var" : "var"
        var_df = _read_h5ad_dataframe(file, var_group)

        # === STANDARDIZE METADATA ===
        # Cell barcodes
        if "_index" in names(obs_df)
            rename!(obs_df, "_index" => "barcode")

        elseif "index" in names(obs_df)
            rename!(obs_df, "index" => "barcode")

        else
            obs_df.barcode = ["cell_$i" for i in 1:nrow(obs_df)]
        end

        # Gene names (most common failure point)
        if "_index" in names(var_df)
            rename!(var_df, "_index" => "gene_name")
        elseif hasproperty(var_df, :gene_names)
            rename!(var_df, :gene_names => "gene_name")
        elseif hasproperty(var_df, :index)
            rename!(var_df, :index => "gene_name")
        elseif ncol(var_df) > 0
            # Fallback: take first column as gene names
            rename!(var_df, 1 => "gene_name")
        else
            error("Could not find any gene names in var group.")
        end

        # Ensure gene_id exists
        if !("gene_id" in names(var_df))
            var_df[!, "gene_id"] = var_df.gene_name
        end

        # Final safety check
        if !("gene_name" in names(var_df)) || nrow(var_df) == 0
            error("Could not determine gene names from .h5ad var group.")
        end

        # 4. Create Object
        obj = SingleCellObject(counts, obs_df, var_df)

        # 5. Load Embeddings (obsm) if available
        if haskey(file, "obsm")
            for key in keys(file["obsm"])
                clean_key = replace(key, "X_" => "")
                embedding = read(file["obsm"][key])

                if size(embedding, 2) == size(obj.counts, 2)
                    obj.reductions[clean_key] = embedding
                elseif size(embedding, 1) == size(obj.counts, 2)
                    obj.reductions[clean_key] = copy(transpose(embedding))
                else
                    println("Warning: Skipping reduction '$key' due to dimension mismatch: $(size(embedding)). Expected $(size(obj.counts, 2)) cells.")
                end
            end
        end

        println("Successfully loaded $(size(obj.counts, 2)) cells and $(size(obj.counts, 1)) genes.")
        return obj
    end
end


# =============================================================================
# 10x Genomics Format (H5)
# =============================================================================

"""
    read_h5(filepath::String)

Reads a 10x Genomics HDF5 file (.h5) using native HDF5.jl.
Supports both Cell Ranger v2 and v3 format specifications.
"""
function read_h5(filepath::String)
    # --- 1. Input Validation ---
    if !isfile(filepath)
        error("File '$filepath' does not exist.")
    end
    if !endswith(lowercase(filepath), ".h5")
        println("Warning: File '$filepath' does not have a .h5 extension.")
    end

    println("Loading .h5 from: $filepath")

    h5open(filepath, "r") do file
        # --- 2. Locate the Matrix Group ---
        # 10x v3 files use a "matrix" group. v2 uses the genome name (e.g. "GRCh38").
        root = nothing
        if haskey(file, "matrix")
            root = file["matrix"]
        else
            # Search for a group that looks like a genome (contains "barcodes" and "data")
            for key in keys(file)
                if isa(file[key], HDF5.Group) && haskey(file[key], "barcodes") && haskey(file[key], "data")
                    root = file[key]
                    println("  Detected genome group: $key")
                    break
                end
            end
        end

        if isnothing(root)
            error("Invalid 10x H5 format. Could not find 'matrix' group or valid genome group.")
        end

        # --- 3. Read Sparse Matrix Components ---
        # 10x stores data in CSC format (Compressed Sparse Column)
        data = read(root["data"])
        indices = read(root["indices"]) .+ 1 # Convert 0-based (Python/C) to 1-based (Julia)
        indptr = read(root["indptr"]) .+ 1   # Convert 0-based to 1-based
        shape = read(root["shape"])          # [rows (genes), cols (cells)]

        n_genes = Int(shape[1])
        n_cells = Int(shape[2])

        # Construct the sparse matrix
        counts_matrix = SparseMatrixCSC(n_genes, n_cells, indptr, indices, data)

        # --- 4. Read Metadata ---

        # A. Barcodes (Cells)
        barcodes = read(root["barcodes"])
        cell_meta = DataFrame(barcode=barcodes)

        # B. Features (Genes)
        var_meta = DataFrame()

        # Handle v3 structure (inside "features" group)
        if haskey(root, "features")
            feat_group = root["features"]
            # Read datasets directly
            var_meta[!, :gene_id] = read(feat_group["id"])
            var_meta[!, :gene_name] = read(feat_group["name"])

            # Feature type (e.g., "Gene Expression", "Antibody Capture")
            if haskey(feat_group, "feature_type")
                var_meta[!, :feature_type] = read(feat_group["feature_type"])
            else
                var_meta[!, :feature_type] .= "Gene Expression"
            end

            # Handle v2 structure (datasets "genes" and "gene_names" at root)
        elseif haskey(root, "genes")
            # In v2, 'genes' is usually the ID, 'gene_names' is the symbol
            var_meta[!, :gene_id] = read(root["genes"])
            if haskey(root, "gene_names")
                var_meta[!, :gene_name] = read(root["gene_names"])
            else
                var_meta[!, :gene_name] = var_meta[!, :gene_id]
            end
            var_meta[!, :feature_type] .= "Gene Expression"
        else
            error("Could not find feature/gene definitions in H5 file.")
        end

        # --- 5. Return Object ---
        println("Successfully loaded $n_cells cells and $n_genes genes.")
        return SingleCellObject(counts_matrix, cell_meta, var_meta)
    end
end