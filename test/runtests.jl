using Test
using SiCell
using SparseArrays
using DataFrames
using Random
using Graphs
using Statistics
using Plots # For smoke testing plotting


# Set headless mode for plotting tests to avoid display errors in CI
ENV["GKSwstype"] = "100"

# Set seed for reproducibility
Random.seed!(42)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================
"""
    create_test_object(n_genes=100, n_cells=50; prefix="cell")

Creates a synthetic SingleCellObject for testing.
"""
function create_test_object(n_genes=100, n_cells=50; prefix="cell")
    # Create sparse count matrix with valid structure (integers)
    counts = sprand(n_genes, n_cells, 0.3) * 100
    counts = round.(Int, counts)

    # Create metadata
    meta = DataFrame(
        barcode=["$(prefix)_$i" for i in 1:n_cells],
        batch=fill("Batch1", n_cells) # Default batch
    )

    # Create variable data
    var_data = DataFrame(
        gene_id=["gene_$i" for i in 1:n_genes],
        gene_name=["GENE_$i" for i in 1:n_genes]
    )

    # Add some "mitochondrial" genes
    mito_indices = [1, 5, 10]
    for i in mito_indices
        if i <= n_genes
            var_data.gene_name[i] = "MT-GENE_$i"
        end
    end

    return SingleCellObject(counts, meta, var_data)
end

# ==============================================================================
# MAIN TEST SUITE
# ==============================================================================
@testset "SiCell.jl Comprehensive Suite" begin

    # --- 1. CORE DATA STRUCTURES ---
    @testset "Types & Validation" begin
        obj = create_test_object(10, 5)
        @test obj isa SingleCellObject
        @test size(obj.counts) == (10, 5)

        # Test Subset
        mask = [true, false, true, false, true]
        sub_obj = subset_object(obj, mask)
        @test size(sub_obj.counts, 2) == 3
        @test nrow(sub_obj.meta_data) == 3

        # Error handling
        @test_throws ErrorException subset_object(obj, [1, 100]) # Index out of bounds
    end

    # --- 2. PRE-PROCESSING ---
    @testset "Preprocessing (QC -> Norm -> PCA)" begin
        obj = create_test_object(100, 50)

        # QC
        calculate_qc_metrics!(obj)
        @test "percent_mito" in names(obj.meta_data)

        filter_cells!(obj, min_genes=0, max_mito=100.0) # Permissive filter

        # Norm
        normalize_data!(obj)
        @test !isnothing(obj.norm_data)

        # Features
        find_variable_features!(obj, n_features=20)
        @test sum(obj.var_data.highly_variable) == 20

        # PCA
        run_pca!(obj, ndims=5)
        @test size(obj.reductions["pca"]) == (50, 5)

        # Test Dense path
        scale_data!(obj)
        run_pca!(obj, ndims=5, key="pca_dense")
        @test size(obj.reductions["pca_dense"]) == (50, 5)
    end

    # --- 3. CLUSTERING & GRAPHS ---
    @testset "Clustering & Graph Algorithms" begin
        obj = create_test_object(100, 100)
        normalize_data!(obj)
        find_variable_features!(obj, n_features=50)
        run_pca!(obj, ndims=10)

        # Neighbors
        find_neighbors!(obj, k=5)
        @test haskey(obj.graphs, "neighbors")

        # Method 1: K-Means (Standard)
        run_clustering!(obj, k=3, key="kmeans")
        @test "kmeans" in names(obj.meta_data)

        # Method 2: Louvain (Graph) - Native Implementation
        run_graph_clustering!(obj, method="louvain", resolution=1.0, key="louvain")
        @test "louvain" in names(obj.meta_data)
        @test eltype(obj.meta_data.louvain) <: Integer

        # Method 3: Label Propagation
        run_graph_clustering!(obj, method="label_propagation", key="label_prop")
        @test "label_prop" in names(obj.meta_data)

        # UMAP
        run_umap!(obj, dims=5, n_epochs=10) # Low epochs for speed
        @test size(obj.reductions["umap"]) == (100, 2)
    end

    # --- 4. BATCH CORRECTION ---
    @testset "Batch Integration" begin
        obj1 = create_test_object(100, 50; prefix="A")
        obj1.meta_data.batch .= "BatchA"

        obj2 = create_test_object(100, 50; prefix="B")
        obj2.meta_data.batch .= "BatchB"

        merged = merge_batches!(obj1, obj2, batch1_name="BatchA", batch2_name="BatchB")
        @test size(merged.counts, 2) == 100

        normalize_data!(merged)
        find_variable_features!(merged, n_features=50)
        run_pca!(merged, ndims=10)

        run_harmony!(merged, "batch", max_iter=2, key="harmony")
        @test "harmony" in keys(merged.reductions)

        run_bbknn!(merged, "batch", reduction="harmony", neighbors_within_batch=3)
        @test haskey(merged.graphs, "neighbors")
    end

    # --- 5. TRAJECTORY INFERENCE ---
    @testset "Trajectory (Diffusion & Pseudotime)" begin
        obj = create_test_object(100, 50)
        normalize_data!(obj)
        find_variable_features!(obj, n_features=50)
        run_pca!(obj, ndims=10)
        find_neighbors!(obj, k=10)

        run_diffusion_map!(obj, n_components=5, key="diff")
        @test "diff" in keys(obj.reductions)

        run_pseudotime!(obj, 1, reduction="diff", key="dpt")
        @test "dpt" in names(obj.meta_data)
    end

    # --- 7. Differential expression & Annotation ---
    @testset "Differential Expression" begin
        obj = create_test_object(100, 100)
        obj.meta_data.cluster = rand(1:3, 100)
        normalize_data!(obj)

        markers = find_all_markers(obj, group="cluster")
        @test markers isa DataFrame
    end

    @testset "Annotation" begin

        # Load lightweight test CellMarker database
        test_db_path = joinpath(@__DIR__, "data", "CellMarker_test.csv")

        @test isfile(test_db_path)

        # Create synthetic data
        obj = create_test_object(100, 100)
        obj.meta_data.cluster = rand(1:3, 100)

        normalize_data!(obj)

        # Find markers
        markers = find_all_markers(obj, group="cluster")

        # Load test CellMarker database
        db = load_cellmarker(test_db_path)

        @test db isa DataFrame
        @test nrow(db) > 0

        # Annotation should emit the expected warning
        @test_logs (:warn,) annotate_clusters!(
            obj,
            markers,
            db,
            cluster_col="cluster"
        )

        # Ensure annotation column was added
        @test "cell_type" in names(obj.meta_data)

    end

    # --- 8. PLOTTING (Smoke Tests) ---
    @testset "Plotting Functions" begin
        obj = create_test_object(100, 50)

        normalize_data!(obj)
        find_variable_features!(obj, n_features=50)
        run_pca!(obj, ndims=10)
        run_umap!(obj, n_epochs=5, dims=10)
        obj.meta_data.cluster = rand(1:3, 50)

        @test_nowarn dim_plot(obj, reduction="pca", group="cluster")

        find_neighbors!(obj, k=15, key="pca")
        @test_nowarn paga_plot(obj, cluster_col="cluster")
    end
end