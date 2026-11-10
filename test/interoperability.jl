using Test
using SiCell
using Muon
using SparseArrays
using DataFrames
using Graphs

@testset "Muon.jl AnnData interoperability" begin
    # AnnData uses observations × variables; SiCell uses genes × cells.
    X = sparse(Float64[
        1 0 3 0;
        0 2 0 4;
        5 0 0 1;
    ])
    obs_names = ["cell_A", "cell_B", "cell_C"]
    var_names = ["GENE1", "GENE2", "GENE3", "GENE4"]
    obs = DataFrame(batch=["b1", "b1", "b2"], donor=["d1", "d2", "d3"])
    var = DataFrame(gene_id=["id1", "id2", "id3", "id4"], biotype=fill("protein_coding", 4))
    umap = reshape(collect(1.0:6.0), 3, 2)
    counts_layer = sparse(Matrix(X) .+ 2.0)

    adata = Muon.AnnData(
        X=X,
        obs=obs,
        obs_names=obs_names,
        var=var,
        var_names=var_names,
        obsm=Dict("X_umap" => umap),
        layers=Dict("counts" => counts_layer),
        uns=Dict("source" => "SiCell test"),
    )

    # AnnData's public Muon.jl API accepts DataFrames in obsm/varm. Use
    # nullable DataFrame columns to verify that missing annotation values are
    # preserved without relying on nullable numeric arrays unsupported by the
    # documented AnnData constructor.
    nullable_obsm = DataFrame(
        dim1=Union{Missing,Float64}[0.1, missing, 0.5],
        dim2=Union{Missing,Float64}[0.2, 0.4, 0.6],
    )
    nullable_varm = DataFrame(
        score=Union{Missing,Float64}[1.0, missing, 3.0, 4.0],
    )
    adata.obsm["qc_annotation"] = nullable_obsm
    adata.varm["gene_annotation"] = nullable_varm
    adata.obs[!, :nullable_score] = Union{Missing,Float64}[0.2, missing, 0.8]

    obj = from_anndata(adata)
    @test obj isa SingleCellObject
    @test size(obj.counts) == (4, 3)
    @test Matrix(obj.counts) == Matrix(transpose(X))
    @test obj.meta_data.batch == obs.batch
    @test obj.var_data.gene_name == var_names
    @test obj.var_data.gene_id == var.gene_id
    @test obj.reductions["umap"] == umap

    exported = to_anndata(obj)
    @test size(exported.X) == size(X)
    @test Matrix(exported.X) == Matrix(X)
    @test collect(exported.obs_names) == obs_names
    @test collect(exported.var_names) == var_names
    @test exported.obs.batch == obs.batch
    @test exported.var.gene_id == var.gene_id
    @test Matrix(exported.obsm["X_umap"]) == umap
    @test Matrix(exported.layers["counts"]) == Matrix(counts_layer)
    @test exported.uns["source"] == "SiCell test"
    @test isequal(exported.obsm["qc_annotation"], nullable_obsm)
    @test isequal(exported.varm["gene_annotation"], nullable_varm)
    @test isequal(exported.obs.nullable_score, adata.obs.nullable_score)

    # If nullable numeric arrays enter the SiCell passthrough payload, convert
    # them to DataFrames because Muon.AnnData's public constructor accepts
    # DataFrames rather than nullable numeric arrays in obsm/varm.
    array_payload_obj = from_anndata(adata)
    nullable_array = Union{Missing,Float64}[0.1 0.2; missing 0.4; 0.5 0.6]
    array_payload_obj.anndata["obsm"]["nullable_array"] = nullable_array
    array_payload_export = to_anndata(array_payload_obj)
    @test array_payload_export.obsm["nullable_array"] isa DataFrame
    @test isequal(Matrix(array_payload_export.obsm["nullable_array"]), nullable_array)

    # SiCell's KNN graph is exported using AnnData's conventional obsp keys.
    graph = SimpleGraph(3)
    add_edge!(graph, 1, 2)
    obj.graphs["neighbors"] = graph
    obj.graphs["neighbor_indices"] = [[1, 2], [2, 1], [3, 2]]
    obj.graphs["neighbor_distances"] = [[0.0, 0.25], [0.0, 0.25], [0.0, 0.5]]
    graph_export = to_anndata(obj)
    @test graph_export.obsp["connectivities"][1, 2] == 1.0
    @test graph_export.obsp["distances"][1, 2] == 0.25

    # A named AnnData layer can be selected as the analysis input.
    from_layer = from_anndata(adata; layer="counts")
    @test Matrix(from_layer.counts) == Matrix(transpose(counts_layer))
    @test Matrix(to_anndata(from_layer).X) == Matrix(counts_layer)
    @test_throws ArgumentError from_anndata(adata; layer="missing")

    # Dense AnnData input must follow the same orientation convention.
    dense_adata = Muon.AnnData(
        X=Matrix(X), obs=obs, obs_names=obs_names, var=var, var_names=var_names,
    )
    @test Matrix(from_anndata(dense_adata).counts) == Matrix(transpose(X))

    # Names and observation-aligned metadata survive SiCell cell subsetting.
    subset = subset_object(obj, [true, false, true])
    subset_export = to_anndata(subset)
    @test collect(subset_export.obs_names) == ["cell_A", "cell_C"]
    @test size(subset_export.X) == (2, 4)
    @test size(subset_export.obsm["X_umap"]) == (2, 2)

    # QC filtering must subset reduction rows (cells), as required by AnnData.
    filtered = from_anndata(adata)
    calculate_qc_metrics!(filtered)
    filter_cells!(filtered; min_genes=0, min_counts=5, max_mito=100.0)
    filtered_export = to_anndata(filtered)
    @test collect(filtered_export.obs_names) == ["cell_B", "cell_C"]
    @test Matrix(filtered_export.obsm["X_umap"]) == umap[[2, 3], :]

    # Normalized export keeps counts in the conventional layers['counts'] slot.
    normalize_data!(obj)
    normalized_export = to_anndata(obj; X=:normalized)
    @test size(normalized_export.X) == size(X)
    @test Matrix(normalized_export.layers["counts"]) == Matrix(X)
end
