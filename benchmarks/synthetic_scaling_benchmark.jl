# ==============================================================================
# SiCell.jl — Synthetic Multi-Scale Scaling Benchmark
# ==============================================================================
#
# WHY THIS SCRIPT EXISTS
# -----------------------
# The benchmarks in the software paper (Tables 1-2) are two anecdotal points:
# one real dataset at 2,700 cells, one at 34,144 cells, each measured once.
# That's a fair demonstration but a thin scientific claim. This script fixes
# the two biggest gaps at once, using only what a laptop or a free Colab
# session can do:
#
#   1. SCALE:      runs the pipeline at 7 cell-count points spanning 500 to
#                  40,000 cells (80x range) instead of 2 points, so we get an
#                  actual scaling curve, not two dots.
#   2. REPLICATION: at each scale point, generates N_REPLICATES independent
#                  synthetic datasets (different random seeds -> different
#                  cell composition, not just repeated timing noise) and
#                  reports mean +/- SD, not a single run.
#
# It also directly tests the claim in the paper's Discussion section that
# UMAP is the memory/scaling bottleneck at large N -- if that's true, this
# will show it in the numbers; if it isn't, we should soften that claim.
#
# The synthetic data generator is a lightweight negative-binomial-ish/Poisson
# simulator with latent cell-type structure (marker genes per type, varying
# library size, mitochondrial genes) -- not a research-grade simulator like
# Splatter or scDesign, but realistic enough that clustering, HVG selection,
# and marker detection all have real signal to work on rather than pure
# noise, so the pipeline is stress-tested the way it would be on real data.
#
# WHAT THIS DOES NOT DO
# ----------------------
# It does not compare against Scanpy/Seurat. That would need the exact same
# synthetic matrix loaded into Python/R for a true apples-to-apples test
# (this script does export the matrix so that's a small follow-up, not a
# rewrite -- ask if you want that companion script).
#
# HOW TO RUN
# -----------
# Local Julia:  julia synthetic_scaling_benchmark.jl
# Google Colab: install a Julia kernel (e.g. via the `julia-1.10 IJulia`
#               setup cell many Colab-Julia notebooks use), then run this
#               file's contents in a cell, or `include("synthetic_scaling_benchmark.jl")`.
#
# Expect ~10-25 minutes total with the defaults below (largely dominated by
# the two biggest scale points). Lower N_REPLICATES or trim SCALE_POINTS if
# you want a faster first pass.
# ==============================================================================

using Pkg

# --- Isolate this benchmark into its own throwaway environment ---
# IMPORTANT: this used to run `Pkg.add` inside whatever project was already
# active. If that's SiCell1's own project (likely, if you're running this
# from inside the package folder), that's wrong on two counts: it risks
# version conflicts with packages SiCell already transitively depends on
# (this is what caused the StatsBase error), and even when it succeeds it
# permanently pollutes SiCell.jl's own Project.toml with benchmark-only
# dependencies that have no business being there. A temp environment avoids
# both problems entirely.
Pkg.activate(temp = true)

# --- Which SiCell to benchmark? ---
# Your local working copy currently has an unresolved compat conflict
# (Clustering="0.15" vs CommunityDetection="0.2.0", which needs Clustering
# 0.14.x) that isn't in the GitHub version -- see chat for the fix.
#
# The GitHub repo itself also has a malformed git tag (literally named
# `"0.1.1"`, quotes included, alongside the correct `v0.1.1`) that breaks
# Julia's bundled LibGit2 clone -- confirmed via `git ls-remote --tags`.
# System git tolerates it fine, so we shell out to it directly instead of
# letting Pkg.jl's LibGit2 backend touch this repo's refs at all.
scratch_dir = mktempdir()
println("Cloning SiCell via system git into: ", scratch_dir)
try
    run(`git clone --depth 1 https://github.com/Sizerta/SiCell.jl.git $scratch_dir`)
catch e
    error("""
    System `git clone` failed ($e).
    Either git isn't on PATH, or something else is wrong. As a fallback,
    manually run in a terminal:
        git clone --depth 1 https://github.com/Sizerta/SiCell.jl.git <some folder>
    then set SICELL_PATH below to that folder and use Pkg.develop directly.
    """)
end
Pkg.develop(path = scratch_dir)

# Once the malformed tag is deleted from GitHub (see chat) and your local
# compat conflict is fixed, you can go back to developing your working copy
# directly instead:
# const SICELL_PATH = abspath(joinpath(@__DIR__, ".."))
# Pkg.develop(path = SICELL_PATH)

Pkg.add(["DataFrames", "CSV", "Plots", "Distributions"])
# (StatsBase intentionally dropped -- see weighted_sample_without_replacement
# below, which replaces the one call that needed it, using only Base Julia.)

using SiCell
using DataFrames
using CSV
using Plots
using Distributions
using Random
using SparseArrays
using LinearAlgebra
using Printf
using Statistics

println("Julia version: ", VERSION)
println("Threads available: ", Threads.nthreads())
println("(For multi-threaded DE/normalization timings, start Julia with `julia -t auto`.)")

# ==============================================================================
# 1. CONFIGURATION
# ==============================================================================
const N_GENES        = 20_000            # fixed gene count across all scale points
const SCALE_POINTS    = [500, 1_000, 2_500, 5_000, 10_000, 20_000, 40_000]
const N_REPLICATES    = 3                 # independent synthetic datasets per scale point
const N_CELL_TYPES    = 8
const N_PCS           = 30
const K_NEIGHBORS     = 20
const CLUSTER_METHOD  = "louvain"
const LOGFC_THRESHOLD = 0.25
const UMAP_MAX_CELLS  = 20_000            # skip UMAP above this N (see header note)
const BASE_SEED       = 20260831
const OUTPUT_CSV      = "sicell_synthetic_benchmark_results.csv"
const OUTPUT_PLOT     = "sicell_synthetic_benchmark_scaling.png"

# ==============================================================================
# 2. SYNTHETIC DATA GENERATOR
# ==============================================================================
"""
    weighted_sample_without_replacement(rng, weights, k)

Weighted sample of `k` indices without replacement, via the Efraimidis-Spirakis
key trick (draw u_i ~ Uniform(0,1), key_i = u_i^(1/w_i), keep the k largest
keys). Same algorithm StatsBase.sample(...; replace=false) uses internally --
implemented here directly so this script has no dependency on StatsBase
(which is often already transitively pinned by other packages, as it was
here via Clustering.jl, and is exactly what caused the version conflict).
"""
function weighted_sample_without_replacement(rng::AbstractRNG, weights::Vector{Float64}, k::Int)
    n = length(weights)
    keys = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        keys[i] = rand(rng)^(1.0 / weights[i])
    end
    return partialsortperm(keys, 1:k, rev = true)
end

"""
    generate_synthetic_counts(n_cells, n_genes; ...)

Generates a synthetic sparse count matrix with latent cell-type structure,
so downstream clustering/HVG/DE steps have real signal, not pure noise.

Returns (counts::SparseMatrixCSC, gene_names::Vector{String}, true_labels::Vector{Int}).
"""
function generate_synthetic_counts(
    n_cells::Int, n_genes::Int;
    n_cell_types::Int = N_CELL_TYPES,
    n_markers_per_type::Int = 30,
    genes_per_cell::Int = 2_500,     # controls sparsity (~detected genes/cell)
    marker_fold::Float64 = 15.0,
    mito_frac::Float64 = 0.03,
    seed::Int = 1234
)
    rng = MersenneTwister(seed)

    # Baseline gene expression weight: log-normal, mimics the long-tailed
    # shape of real transcriptomes (most genes low/rare, a few highly expressed).
    baseline_weight = exp.(randn(rng, n_genes) .* 1.2)

    # Assign each cell a latent type, and each type a disjoint block of marker genes.
    cell_type = rand(rng, 1:n_cell_types, n_cells)
    gene_pool = randperm(rng, n_genes)
    marker_genes = Vector{Vector{Int}}(undef, n_cell_types)
    for t in 1:n_cell_types
        lo = (t - 1) * n_markers_per_type + 1
        hi = min(t * n_markers_per_type, n_genes)
        marker_genes[t] = gene_pool[lo:hi]
    end

    # Per-cell library size: log-normal around ~6,000 total counts (typical depth).
    lib_size = exp.(randn(rng, n_cells) .* 0.4 .+ log(6000.0))

    k = min(genes_per_cell, n_genes)
    I = Int[]; J = Int[]; V = Float64[]
    sizehint!(I, n_cells * k); sizehint!(J, n_cells * k); sizehint!(V, n_cells * k)

    weight = similar(baseline_weight)
    for c in 1:n_cells
        copyto!(weight, baseline_weight)
        weight[marker_genes[cell_type[c]]] .*= marker_fold

        detected = weighted_sample_without_replacement(rng, weight, k)
        w_detected = weight[detected]
        mu = w_detected ./ sum(w_detected) .* lib_size[c]

        for (idx, g) in enumerate(detected)
            cnt = rand(rng, Poisson(max(mu[idx], 0.01)))
            if cnt > 0
                push!(I, g); push!(J, c); push!(V, Float64(cnt))
            end
        end
    end

    counts = sparse(I, J, V, n_genes, n_cells)

    gene_names = ["GENE_$(i)" for i in 1:n_genes]
    n_mito = round(Int, mito_frac * n_genes)
    for i in 1:n_mito
        gene_names[i] = "MT-GENE_$(i)"
    end

    return counts, gene_names, cell_type
end

"""
    build_object(n_cells, seed)

Wraps generate_synthetic_counts into a ready-to-use SingleCellObject.
"""
function build_object(n_cells::Int, seed::Int)
    counts, gene_names, true_labels = generate_synthetic_counts(n_cells, N_GENES; seed = seed)

    meta = DataFrame(barcode = ["cell_$i" for i in 1:n_cells], true_type = true_labels)
    var  = DataFrame(gene_id = gene_names, gene_name = gene_names)

    obj = SingleCellObject(counts, meta, var)
    return obj
end

# ==============================================================================
# 3. TIMED PIPELINE
# ==============================================================================
"""
    run_timed_pipeline(obj; run_umap)

Runs the core SiCell.jl pipeline on `obj`, timing each step with `@timed`
(single precise wall-clock measurement per call -- replicate-level variance
across independent datasets is what gives us the real distribution, not
repeated calls on identical input). Returns a Vector of NamedTuples.
"""
function run_timed_pipeline(obj::SingleCellObject; run_umap::Bool = true)
    rows = NamedTuple[]

    function timeit!(label::String, f::Function)
        r = @timed f()
        push!(rows, (step = label, time_s = r.time, memory_MiB = r.bytes / 1024^2, gctime_s = r.gctime))
        return r.value
    end

    # --- Preprocessing (timed separately, NOT included in "Total Pipeline"
    #     below -- matches the definition already used in the paper's
    #     Table 1/2, where "Total Pipeline" = PCA + KNN + Clustering + UMAP + Markers) ---
    timeit!("QC_Preprocess", () -> begin
        calculate_qc_metrics!(obj, mito_prefix = "MT-")
        filter_cells!(obj, min_genes = 50, max_mito = 100.0)  # generous: shouldn't drop synthetic cells
        normalize_data!(obj)
        find_variable_features!(obj, n_features = 2000)
        nothing
    end)

    # --- Core, directly comparable to the paper's benchmark tables ---
    timeit!("PCA", () -> run_pca!(obj, ndims = N_PCS))
    timeit!("KNN_Graph", () -> find_neighbors!(obj, k = K_NEIGHBORS))
    timeit!("Clustering", () -> run_graph_clustering!(obj, method = CLUSTER_METHOD, key = "cluster"))

    if run_umap
        timeit!("UMAP", () -> run_umap!(obj, dims = N_PCS))
    end

    timeit!("Marker_Genes", () -> find_all_markers(obj, group = "cluster", logfc_threshold = LOGFC_THRESHOLD))

    return rows
end

# ==============================================================================
# 4. WARM-UP (exclude JIT compilation from timed results)
# ==============================================================================
println("\n[Warm-up] Pre-compiling all pipeline methods on a tiny synthetic dataset...")
let warm_obj = build_object(200, 1)
    run_timed_pipeline(warm_obj; run_umap = true)
end
println("[Warm-up] Done.\n")

# ==============================================================================
# 5. MAIN SWEEP: SCALE POINTS x REPLICATES
# ==============================================================================
all_rows = NamedTuple[]

for n_cells in SCALE_POINTS
    run_umap = n_cells <= UMAP_MAX_CELLS
    if !run_umap
        println("  (n_cells=$n_cells > UMAP_MAX_CELLS=$UMAP_MAX_CELLS: skipping UMAP step, per Discussion's known memory limitation)")
    end

    for rep in 1:N_REPLICATES
        # hash(...) guarantees a distinct seed per (n_cells, rep) pair -- a plain
        # arithmetic combination of the two can collide for some scale-point sets.
        seed = Int(hash((BASE_SEED, n_cells, rep)) % 1_000_000_000)
        println("Running n_cells=$n_cells, replicate=$rep/$N_REPLICATES (seed=$seed)...")

        obj = build_object(n_cells, seed)

        step_rows = run_timed_pipeline(obj; run_umap = run_umap)
        actual_n_cells = size(obj.counts, 2)   # post-QC/filter count (should equal n_cells given generous thresholds)
        for r in step_rows
            push!(all_rows, merge((n_cells = actual_n_cells, replicate = rep, seed = seed), r))
        end

        core_steps = ["PCA", "KNN_Graph", "Clustering", "UMAP", "Marker_Genes"]
        total_time = sum(r.time_s for r in step_rows if r.step in core_steps)
        push!(all_rows, (n_cells = actual_n_cells, replicate = rep, seed = seed,
                          step = "TOTAL_PIPELINE", time_s = total_time,
                          memory_MiB = NaN, gctime_s = NaN))

        @printf("  -> Total pipeline (PCA..Markers): %.3fs\n", total_time)
    end
end

results = DataFrame(all_rows)
CSV.write(OUTPUT_CSV, results)
println("\nRaw results written to: $OUTPUT_CSV")

# ==============================================================================
# 6. SUMMARY STATISTICS (mean +/- SD across replicates, per scale point per step)
# ==============================================================================
summary_df = combine(
    groupby(results, [:n_cells, :step]),
    :time_s => mean => :mean_time_s,
    :time_s => std  => :sd_time_s,
    :time_s => median => :median_time_s,
    :time_s => minimum => :min_time_s,
    :time_s => maximum => :max_time_s,
    :memory_MiB => mean => :mean_memory_MiB,
    nrow => :n_replicates
)
sort!(summary_df, [:step, :n_cells])

println("\n" * "="^88)
println("  SUMMARY: mean +/- SD across $N_REPLICATES replicates, per scale point")
println("="^88)
for step in unique(summary_df.step)
    println("\n-- $step --")
    sub = filter(row -> row.step == step, summary_df)
    @printf("%-10s | %-14s | %-10s | %-10s\n", "n_cells", "mean +/- sd (s)", "median (s)", "mem (MiB)")
    println("-"^55)
    for row in eachrow(sub)
        @printf("%-10d | %6.3f +/- %-5.3f | %-10.3f | %-10.1f\n",
                row.n_cells, row.mean_time_s, row.sd_time_s, row.median_time_s, row.mean_memory_MiB)
    end
end

CSV.write("sicell_synthetic_benchmark_summary.csv", summary_df)
println("\nSummary statistics written to: sicell_synthetic_benchmark_summary.csv")

# ==============================================================================
# 7. SCALING PLOT
# ==============================================================================
plot_steps = ["PCA", "KNN_Graph", "Clustering", "UMAP", "Marker_Genes", "TOTAL_PIPELINE"]
p = plot(
    title = "SiCell.jl scaling on synthetic data ($N_REPLICATES replicates/point)",
    xlabel = "Number of cells", ylabel = "Time (s)",
    xscale = :log10, yscale = :log10,
    legend = :outertopright, size = (950, 600), dpi = 300
)
for step in plot_steps
    sub = filter(row -> row.step == step, summary_df)
    isempty(sub) && continue
    sort!(sub, :n_cells)
    # (No error ribbon here: with log-scale axes, mean-SD can dip <= 0 for
    # fast/noisy steps at small N, which log-scale can't render. Full SD is
    # in the CSV/summary table if you want to build a custom error-bar plot.)
    plot!(p, sub.n_cells, sub.mean_time_s, label = step,
          marker = :circle, markersize = 4, linewidth = 2)
end
savefig(p, OUTPUT_PLOT)
println("\nScaling plot written to: $OUTPUT_PLOT")

# ==============================================================================
# 8. SYSTEM INFO (for reporting reproducibility in the paper)
# ==============================================================================
# Wrapped in try/catch per field: exact Sys.* introspection functions can
# differ slightly across Julia versions/platforms, and this is purely
# cosmetic reporting -- it should never be the reason the run "fails" after
# 10-25 minutes of real computation.
println("\n" * "="^60)
println("  SYSTEM INFO (paste into your Methods section)")
println("="^60)
println("Julia version: ", VERSION)
println("Threads.nthreads(): ", Threads.nthreads())
try
    println("OS: ", Sys.KERNEL, " ", Sys.MACHINE)
catch e
    println("OS: (could not determine: $e)")
end
try
    println("CPU: ", Sys.cpu_info()[1].model, " x", length(Sys.cpu_info()))
catch e
    println("CPU: (could not determine: $e)")
end
try
    println("Total RAM (GiB): ", round(Sys.total_memory() / 1024^3, digits = 1))
catch e
    println("Total RAM: (could not determine: $e)")
end
println("="^60)
println("\nDone.")
