# ==============================================================================
# SiCell Annotation Module
# Contains: Database Loading, Automatic Cluster Annotation, and Visualization
# ==============================================================================

# Note: This file assumes the following standard libraries are loaded in the parent module:
# DataFrames, CSV, Statistics, Plots, StatsPlots

# ==============================================================================
# 1. Database Loading
# ==============================================================================

"""
    load_panglaodb([filepath::String])

Reads the PanglaoDB file. 
If `filepath` is not provided, it looks for "panglao.tsv" in the package's `data/` directory.
Returns a DataFrame with normalized column names.
"""
function load_panglaodb(filepath::String="")
    # 1. Handle auto-discovery of the data file
    if isempty(filepath)
        root_dir = pkgdir(@__MODULE__)
        if isnothing(root_dir)
            # Fallback for scripts/REPL where pkgdir might be undefined
            filepath = "panglao.tsv"
        else
            filepath = joinpath(root_dir, "data", "panglao.tsv")
        end
        println("Using default database path: $filepath")
    end

    # 2. Validation
    if !isfile(filepath)
        error("""
        Database file not found at '$filepath'. 
        Please create a 'data' folder in your package root and place 'panglao.tsv' there.
        """)
    end

    # 3. Read and Normalize
    df = CSV.read(filepath, DataFrame, silencewarnings=true)

    # Normalize column names: lowercase and replace spaces with underscores
    rename!(df, names(df) .=> lowercase.(replace.(names(df), " " => "_")))

    return df
end


"""
    load_cellmarker([filepath::String]; species="Hs", cancer_type=nothing, tissue_type=nothing)

Loads the CellMarker 2.0 database and normalizes it to the same schema
used by `annotate_clusters!`, making it a drop-in replacement for
`load_panglaodb()`.

Download from: http://117.50.127.228/CellMarker/CellMarker_download.html
File: All_Cell_marker.tsv or All_Cell_marker.csv

# Arguments
- `species`: `"Hs"` for Human (default), `"Mm"` for Mouse, `"all"` for both.
- `cancer_type`: Optionally narrow to a cancer type e.g. `"Breast Cancer"`, `"Normal"`.
- `tissue_type`: Optionally narrow to a tissue e.g. `"Bone marrow"`, `"Blood"`.

Specifying only `species` gives researchers full freedom across all tissues
and cancer types for that organism — the most common use case.

# Returns
A DataFrame with normalized columns compatible with `annotate_clusters!`:
- `species`              — "Hs" or "Mm"
- `official_gene_symbol` — gene symbol (uppercase)
- `cell_type`            — cell name
- `organ`                — tissue type

# Example
```julia
# Human (default) — all tissues and cancer types
db = load_cellmarker()

# Mouse — all tissues
db = load_cellmarker(species="Mm")

# Mouse bone marrow specifically
db = load_cellmarker(species="Mm", tissue_type="Bone marrow")

# Human breast cancer cells
db = load_cellmarker(species="Hs", cancer_type="Breast Cancer")

# All species combined
db = load_cellmarker(species="all")
```
"""
function load_cellmarker(filepath::String="";
    species::String="Hs",
    cancer_type::Union{String,Nothing}=nothing,
    tissue_type::Union{String,Nothing}=nothing)

    # ── Auto-discover file ────────────────────────────────────────────────
    if isempty(filepath)
        root_dir = pkgdir(@__MODULE__)
        candidates = [
            joinpath(isnothing(root_dir) ? "." : root_dir, "Data", "Cell_marker_All.csv"),
            joinpath(isnothing(root_dir) ? "." : root_dir, "data", "Cell_marker_All.csv"),
            joinpath(isnothing(root_dir) ? "." : root_dir, "Data", "Cell_marker_All.tsv"),
            joinpath(isnothing(root_dir) ? "." : root_dir, "data", "Cell_marker_All.tsv"),
            joinpath(isnothing(root_dir) ? "." : root_dir, "Data", "cellmarker2.tsv"),
            joinpath(isnothing(root_dir) ? "." : root_dir, "data", "cellmarker2.tsv"),
            "Cell_marker_All.csv",
            "Cell_marker_All.tsv",
            "cellmarker2.tsv"
        ]
        filepath = something(findfirst(isfile, candidates) |>
                             i -> isnothing(i) ? nothing : candidates[i], "")
        isempty(filepath) &&
            error("CellMarker 2.0 file not found. Download Cell_marker_All.tsv " *
                  "from http://117.50.127.228/CellMarker/CellMarker_download.html " *
                  "and place it in your package data/ directory.")
        println("Using CellMarker 2.0 database: $filepath")
    end

    isfile(filepath) || error("File not found: $filepath")

    # ── Read ─────────────────────────────────────────────────────────────
    delim = endswith(lowercase(filepath), ".csv") ? ',' : '\t'
    df = CSV.read(filepath, DataFrame, delim=delim, silencewarnings=true)
    rename!(df, names(df) .=> lowercase.(replace.(String.(names(df)), " " => "_")))

    # ── Normalize species column to "Hs" / "Mm" ──────────────────────────
    # CellMarker uses "Human" / "Mouse" — map to PanglaoDB convention
    if "species" in names(df)
        df[!, :species] = map(df.species) do s
            s = string(s)
            occursin("Human", s) || occursin("human", s) ? "Hs" :
            occursin("Mouse", s) || occursin("mouse", s) ? "Mm" : s
        end
    end

    # ── Species filter ───────────────────────────────────────────────────
    if lowercase(species) != "all" && "species" in names(df)
        target = uppercase(species) == "MM" ? "Mm" : "Hs"
        df = filter(row -> row.species == target, df)
        nrow(df) == 0 &&
            error("No entries found for species='$species'. Use \"Hs\", \"Mm\", or \"all\".")
    end

    # ── Optional filters ─────────────────────────────────────────────────
    if !isnothing(cancer_type)
        col = "cancer_type" in names(df) ? :cancer_type : nothing
        if !isnothing(col)
            df = filter(row -> occursin(cancer_type,
                    coalesce(string(row[col]), "")), df)
        end
    end

    if !isnothing(tissue_type)
        col = "tissue_type" in names(df) ? :tissue_type : nothing
        if !isnothing(col)
            df = filter(row -> occursin(tissue_type,
                    coalesce(string(row[col]), "")), df)
        end
    end

    nrow(df) == 0 &&
        error("CellMarker database is empty after filtering. " *
              "Check your cancer_type and tissue_type arguments.")

    # ── Normalize to PanglaoDB-compatible schema ──────────────────────────
    # PanglaoDB schema: species | official_gene_symbol | cell_type | organ
    # CellMarker cols:  species | Symbol/marker         | cell_name  | tissue_type
    gene_col = "symbol" in names(df) ? :symbol :
               "marker" in names(df) ? :marker : nothing
    name_col = "cell_name" in names(df) ? :cell_name :
               "cell_type" in names(df) ? :cell_type : nothing
    organ_col = "tissue_type" in names(df) ? :tissue_type :
                "tissue_class" in names(df) ? :tissue_class : nothing

    isnothing(gene_col) && error("Cannot find gene symbol column in CellMarker file.")
    isnothing(name_col) && error("Cannot find cell name column in CellMarker file.")

    normalized = DataFrame(
        species=df.species,
        official_gene_symbol=uppercase.(coalesce.(string.(df[!, gene_col]), "")),
        cell_type=coalesce.(string.(df[!, name_col]), "Unknown"),
        organ=isnothing(organ_col) ? fill("", nrow(df)) :
              coalesce.(string.(df[!, organ_col]), "")
    )

    # Drop rows with missing gene symbols
    filter!(row -> !isempty(row.official_gene_symbol), normalized)

    println("CellMarker 2.0 loaded: $(nrow(normalized)) entries, " *
            "$(length(unique(normalized.cell_type))) cell types, " *
            "$(length(unique(normalized.species))) species.")

    return normalized
end

# ==============================================================================
# 2. Annotation Logic
# ==============================================================================

"""
    annotate_clusters!(obj, markers, db; ...)

Annotates clusters based on marker gene overlap with a reference database (PanglaoDB).
Updates `obj.meta_data` with a new column (default: "cell_type").

# Arguments
- `obj`: SingleCellObject.
- `markers`: DataFrame of marker genes (must contain :cluster and :gene).
- `db`: Reference database DataFrame.
- `species`: Filter DB by species (e.g., "Hs", "Mm").
- `organ`: Filter DB by organ (optional).
- `min_score`: Minimum Jaccard score required to assign a label.
"""
function annotate_clusters!(obj::SingleCellObject, markers::DataFrame, db::DataFrame;
    species::String="Hs",
    organ::Union{String,Nothing}=nothing,
    min_score::Float64=0.05,
    cluster_col::String="cluster",
    new_label_col::String="cell_type")

    @warn "Automatic annotation is based on healthy references. Tumor/disease datasets may be misannotated."

    # 1. Calculate Scores
    # Returns Dict: Cluster ID -> Vector of (Label, Score)
    scores_dict = _calculate_annotation_scores(markers, db, species, organ)

    # 2. Assign Best Labels
    cluster_ids = unique(markers.cluster)
    annotation_map = Dict{Any,String}()

    for clus in cluster_ids
        # Get the best match tuple (Label, Score)
        if isempty(scores_dict[clus])
            annotation_map[clus] = "Unknown"
            continue
        end

        top_match = scores_dict[clus][1]
        label = top_match[2] >= min_score ? top_match[1] : "Unknown"
        annotation_map[clus] = label
    end

    # 3. Update Metadata
    if !(cluster_col in names(obj.meta_data))
        error("Cluster column '$cluster_col' not found in metadata.")
    end

    current_clusters = obj.meta_data[!, cluster_col]
    new_labels = [get(annotation_map, c, "Unknown") for c in current_clusters]
    obj.meta_data[!, new_label_col] = new_labels

    println("Annotation complete. Labels stored in '$new_label_col'.")
    return annotation_map
end

# ==============================================================================
# 3. Visualization
# ==============================================================================

"""
    plot_annotation_heatmap(markers, db; top_n=5)

Visualizes the confidence of cell type annotations. 
Creates a heatmap of Clusters (x) vs. Potential Cell Types (y), colored by Jaccard Score.
Only shows cell types that appear in the top `top_n` candidates for at least one cluster.
"""
function plot_annotation_heatmap(markers::DataFrame, db::DataFrame;
    species::String="Hs",
    organ::Union{String,Nothing}=nothing,
    top_n::Int=3)

    println("Generating annotation confidence heatmap...")

    # 1. Calculate all scores
    scores_dict = _calculate_annotation_scores(markers, db, species, organ)
    clusters = sort(collect(keys(scores_dict)))

    # 2. Identify "Relevant" Cell Types to plot
    # We filter to keep the plot readable (avoiding 1000+ rows)
    relevant_types = Set{String}()

    for clus in clusters
        # Get top N candidates for this cluster
        candidates = scores_dict[clus][1:min(top_n, length(scores_dict[clus]))]
        for (label, score) in candidates
            if score > 0.01 # Filter out zero/noise matches
                push!(relevant_types, label)
            end
        end
    end

    sorted_types = sort(collect(relevant_types))

    if isempty(sorted_types)
        error("No overlapping cell types found to plot.")
    end

    # 3. Build Matrix for Heatmap (Cell Types x Clusters)
    score_matrix = zeros(Float64, length(sorted_types), length(clusters))

    for (j, clus) in enumerate(clusters)
        clus_scores = Dict(scores_dict[clus]) # Convert vector of tuples to Dict for lookup
        for (i, type_label) in enumerate(sorted_types)
            score_matrix[i, j] = get(clus_scores, type_label, 0.0)
        end
    end

    # 4. Plot
    p = heatmap(
        string.(clusters),
        sorted_types,
        score_matrix,
        xlabel="Cluster",
        ylabel="Cell Type",
        title="Annotation Confidence (Jaccard Index)",
        color=:viridis,
        aspect_ratio=:auto,
        xrotation=45,
        size=(600, max(400, length(sorted_types) * 20)), # Dynamic height based on rows
        margins=5Plots.mm
    )

    return p
end

# ==============================================================================
# 4. Internal Helper Functions
# ==============================================================================

function _calculate_annotation_scores(markers, db, species, organ)
    # Filter Database
    species_mask = occursin.(species, db.species)
    db_filtered = db[species_mask, :]

    if !isnothing(organ)
        # coalesce handles missing values in 'organ' column
        organ_mask = occursin.(organ, coalesce.(db_filtered.organ, ""))
        db_filtered = db_filtered[organ_mask, :]
    end

    if nrow(db_filtered) == 0
        error("Reference database is empty after filtering for species='$species' and organ='$organ'.")
    end

    # Build Reference Sets (Cell Type -> Set of Genes)
    ref_sets = Dict{String,Set{String}}()
    for row in eachrow(db_filtered)
        ctype = row.cell_type
        gene = uppercase(row.official_gene_symbol)
        if !haskey(ref_sets, ctype)
            ref_sets[ctype] = Set{String}()
        end
        push!(ref_sets[ctype], gene)
    end

    # Prepare for scoring
    cluster_ids = unique(markers.cluster)
    results = Dict{Any,Vector{Tuple{String,Float64}}}()

    # Thread-safe lock for writing to 'results'
    dict_lock = ReentrantLock()

    println("Calculating annotation scores for $(length(cluster_ids)) clusters using $(Threads.nthreads()) threads...")

    Threads.@threads for clus in cluster_ids
        # Extract top 50 markers for this cluster
        clus_markers = filter(row -> row.cluster == clus, markers)
        top_genes = uppercase.(clus_markers.gene[1:min(50, nrow(clus_markers))])
        top_genes_set = Set(top_genes)

        scores = Vector{Tuple{String,Float64}}()

        # Calculate Jaccard Index against all reference cell types
        for (ref_label, ref_genes) in ref_sets
            intersection = length(intersect(top_genes_set, ref_genes))
            union_len = length(union(top_genes_set, ref_genes))
            score = union_len > 0 ? intersection / union_len : 0.0
            push!(scores, (ref_label, score))
        end

        # Sort descending by score
        sort!(scores, by=x -> x[2], rev=true)

        # Thread-safe write
        lock(dict_lock) do
            results[clus] = scores
        end
    end

    return results
end