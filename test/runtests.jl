using Test
using SiCell
using SparseArrays
using DataFrames
using Random
using Graphs
using Statistics
using Plots # For smoke testing plotting
using Muon  # For h5ad round-trip regression tests


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

        run_centroid_alignment!(merged, "batch", max_iter=2, key="centroid")
        @test "centroid" in keys(merged.reductions)

        # Deprecated alias still works and emits a warning
        # The deprecated alias keeps its historical default reduction key.
        @test_warn "deprecated" run_harmony!(merged, "batch", max_iter=1)
        @test "harmony" in keys(merged.reductions)

        # Explicit keys remain supported for callers migrating gradually.
        @test_warn "deprecated" run_harmony!(merged, "batch", max_iter=1, key="legacy")
        @test "legacy" in keys(merged.reductions)

        run_bbknn!(merged, "batch", reduction="centroid", neighbors_within_batch=3)
        @test haskey(merged.graphs, "neighbors")
        @test ne(merged.graphs["neighbors"]) > 0
        ms = calculate_mixing_score(merged; graph_key="neighbors", batch_key="batch")
        @test 0.0 <= ms <= 1.0
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

    # --- 6. TRAJECTORY UNCERTAINTY FRAMEWORK (TUF) ---
    @testset "Trajectory Uncertainty (TES/TDS/LPS)" begin
        obj = create_test_object(100, 50)
        normalize_data!(obj)
        find_variable_features!(obj, n_features=50)
        run_pca!(obj, ndims=10)
        find_neighbors!(obj, k=10)
        run_diffusion_map!(obj, n_components=5)
        run_pseudotime!(obj, 1)

        trajectory_uncertainty!(obj)
        @test "traj_unc_tes" in names(obj.meta_data)
        @test "traj_unc_tds" in names(obj.meta_data)
        @test "traj_unc_lps" in names(obj.meta_data)

        # TES and TDS are both derived quantities bounded in [0,1] when
        # pseudotime is in [0,1] (see docstring / software paper Section 2.9)
        @test all(0.0 .<= obj.meta_data.traj_unc_tes .<= 1.0)
        @test all(0.0 .<= obj.meta_data.traj_unc_tds .<= 1.0)
        @test eltype(obj.meta_data.traj_unc_lps) == Float64

        # Bounds-check guard: pseudotime outside [0,1] must error, not silently compute
        obj.meta_data.bad_pt = obj.meta_data.pseudotime .+ 10.0
        @test_throws ErrorException trajectory_uncertainty!(obj, pseudotime_key="bad_pt")
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

        # Annotation should emit the expected warnings: the healthy-reference
        # disclaimer plus the empty-marker guard (with per-cluster BH over all
        # genes and the default thresholds, this synthetic data yields no
        # significant markers, so the guard path is taken deterministically
        # with the fixed test seed).
        @test_logs (:warn,) (:warn,) annotate_clusters!(
            obj,
            markers,
            db,
            cluster_col="cluster"
        )

        # Ensure annotation column was added (all cells labeled Unknown)
        @test "cell_type" in names(obj.meta_data)
        @test all(==( "Unknown"), obj.meta_data.cell_type)

        # A typo in the destination metadata key must be reported even when
        # the marker table is empty; it must not be masked by the empty-marker
        # fallback.
        empty_markers = DataFrame(gene=String[], cluster=String[])
        @test_throws ArgumentError annotate_clusters!(
            obj, empty_markers, db; cluster_col="missing_cluster_key"
        )
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

    # --- 9. REGRESSION: differential expression statistics ---
    # avg_log2FC must follow the Scanpy linear-scale convention, and p_val_adj
    # must be a per-cluster BH over ALL genes (independent of thresholds).
    @testset "Regression: DE statistics" begin
        # Controlled 3-gene, 16-cell object (8 cells per group). Gene 1 is
        # high in A, gene 2 is high in B, gene 3 is identical in both.
        g1A = collect(15:22);  g2A = 30 .- g1A
        g1B = collect(3:10);   g2B = 30 .- g1B
        g3A = collect(5:12);   g3B = collect(5:12)
        # 3 genes x 16 cells: rows are genes, columns are cells (A cells first)
        counts = permutedims([g1A g2A g3A; g1B g2B g3B])
        meta = DataFrame(
            barcode=["cell_$i" for i in 1:16],
            cluster=[fill("A", 8); fill("B", 8)]
        )
        var = DataFrame(gene_id=["g1", "g2", "g3"], gene_name=["GENE1", "GENE2", "GENE3"])
        obj = SingleCellObject(counts, meta, var)
        normalize_data!(obj)

        m1 = find_all_markers(obj, group="cluster", logfc_threshold=0.0, min_pct=0.0)
        @test "p_val_adj" in names(m1)

        # Scanpy-convention log2FC for GENE1 in cluster A, recomputed
        # independently from the normalized data: means of the log1p values
        # first, then expm1 of the means (this is what Scanpy does)
        row = only(filter(r -> r.gene == "GENE1" && r.cluster == "A", m1))
        m_log_A = mean(vec(obj.norm_data[1, 1:8]))
        m_log_B = mean(vec(obj.norm_data[1, 9:16]))
        expected = log2((expm1(m_log_A) + 1e-9) / (expm1(m_log_B) + 1e-9))
        @test isapprox(row.avg_log2FC, expected; rtol=1e-8)
        @test row.avg_log2FC > 0

        # p_val_adj must not depend on the filtering thresholds: with a
        # per-cluster BH over all genes, the adjusted value for a gene that
        # passes both threshold sets is identical across runs. (The previous
        # pooled, post-filter adjustment changed with the thresholds.)
        m2 = find_all_markers(obj, group="cluster", logfc_threshold=1.0, min_pct=0.0)
        row2 = only(filter(r -> r.gene == "GENE1" && r.cluster == "A", m2))
        @test row.p_val == row2.p_val
        @test row.p_val_adj == row2.p_val_adj

        # Adjusted p-values stay in (0, 1] and never undercut the raw p-value
        @test all(r -> 0.0 < r.p_val_adj <= 1.0, eachrow(m1))
        @test all(r -> r.p_val_adj >= r.p_val, eachrow(m1))

        # GENE3 is identical in both groups: it must NOT come out significant
        @test !("GENE3" in m1.gene)
    end

    # --- 10. REGRESSION: h5ad reader orientation ---
    # AnnData obsm embeddings are (cells x dims) and must be stored that way;
    # AnnData dense X is (cells x genes) and must be transposed to genes x cells.
    @testset "Regression: h5ad orientation" begin
        X = sparse(Float64[1 0 3 0; 0 2 0 4; 5 0 0 1])  # 3 cells x 4 genes
        obs = DataFrame(batch=["b1", "b1", "b2"])
        var = DataFrame(gene_id=["id$i" for i in 1:4])
        umap = [1.0 2.0; 3.0 4.0; 5.0 6.0]              # 3 cells x 2 dims
        adata = Muon.AnnData(
            X=X, obs=obs, obs_names=["c1", "c2", "c3"],
            var=var, var_names=["g$i" for i in 1:4],
            obsm=Dict{String,Matrix{Float64}}("X_umap" => umap),
        )
        path = joinpath(tempdir(), "sicell_orientation_test.h5ad")
        Muon.writeh5ad(path, adata)

        obj = read_h5ad(path)
        # Embeddings must come back in cells x dims orientation, not transposed
        @test size(obj.reductions["umap"]) == (3, 2)
        @test obj.reductions["umap"] == umap
        # Sparse X path: genes x cells
        @test size(obj.counts) == (4, 3)

        # Float32 obsm entries (common in Scanpy exports) must convert cleanly
        pca32 = reshape(Float32.(collect(1.0f0:12.0f0)), 3, 4)  # 3 cells x 4 dims
        adata32 = Muon.AnnData(
            X=X, obs=obs, obs_names=["c1", "c2", "c3"],
            var=var, var_names=["g$i" for i in 1:4],
            obsm=Dict{String,Matrix{Float32}}("X_pca" => pca32),
        )
        path32 = joinpath(tempdir(), "sicell_orientation_test_f32.h5ad")
        Muon.writeh5ad(path32, adata32)
        obj32 = read_h5ad(path32)
        @test size(obj32.reductions["pca"]) == (3, 4)
        @test obj32.reductions["pca"] == Float64.(pca32)

        # Dense X path must also be transposed to genes x cells
        dense_adata = Muon.AnnData(
            X=Matrix(X), obs=obs, obs_names=["c1", "c2", "c3"],
            var=var, var_names=["g$i" for i in 1:4],
        )
        dense_path = joinpath(tempdir(), "sicell_dense_test.h5ad")
        Muon.writeh5ad(dense_path, dense_adata)
        dense_obj = read_h5ad(dense_path)
        @test size(dense_obj.counts) == (4, 3)
        @test Matrix(dense_obj.counts) == Matrix(transpose(X))

        # A square dense matrix is dimensionally ambiguous. Require the caller
        # to specify the orientation returned by HDF5.jl. Muon.jl writes the
        # var × obs layout used by the second branch above.
        square_X = Float64[1 2 3; 4 5 6; 7 8 9]
        square_adata = Muon.AnnData(
            X=square_X,
            obs=DataFrame(batch=["b1", "b1", "b2"]),
            obs_names=["s1", "s2", "s3"],
            var=DataFrame(gene_id=["g1", "g2", "g3"]),
            var_names=["G1", "G2", "G3"],
        )
        square_path = joinpath(tempdir(), "sicell_square_dense_test.h5ad")
        Muon.writeh5ad(square_path, square_adata)
        @test_throws ArgumentError read_h5ad(square_path; raw=false)
        square_obj = read_h5ad(square_path; raw=false, dense_orientation=:var_by_obs)
        @test Matrix(square_obj.counts) == Matrix(transpose(square_X))
        @test_throws ArgumentError read_h5ad(dense_path; raw=false, dense_orientation=:invalid)

        rm(path; force=true)
        rm(path32; force=true)
        rm(dense_path; force=true)
        rm(square_path; force=true)
    end

    # This fixture is generated by the Python `anndata` package in CI. It is
    # intentionally generated outside Muon.jl so the H5AD reader is checked
    # against a file produced by the reference Python implementation. Local
    # test runs may omit the fixture unless the generator script is run first.
    python_fixture = joinpath(@__DIR__, "data", "python_anndata_fixture.h5ad")
    if isfile(python_fixture)
        @testset "Python anndata H5AD fixture" begin
            python_obj = read_h5ad(python_fixture; raw=false)
            @test size(python_obj.counts) == (4, 3)
            @test collect(python_obj.meta_data.barcode) == ["py_cell_1", "py_cell_2", "py_cell_3"]
            @test collect(python_obj.meta_data.batch) == ["batch_a", "batch_a", "batch_b"]
            @test collect(python_obj.var_data.gene_name) == ["PY_G1", "PY_G2", "PY_G3", "PY_G4"]
            @test python_obj.reductions["pca"] ≈ Float64[1 0.5; 2 1.5; 3 2.5]
            @test Matrix(python_obj.counts) == Float64[1 0 5; 0 2 0; 3 0 0; 0 4 1]
        end
    else
        @test_skip isfile(python_fixture)
    end
end

# Optional integration with Muon.jl (installed as a test-only dependency).
include("interoperability.jl")
