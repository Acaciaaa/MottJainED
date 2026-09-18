#!/usr/bin/env julia

import Pkg

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(PROJECT_ROOT; io=devnull)

using Dates
using FuzzifiED
using MottJainED
using Random
using SHA
using TOML

include(joinpath(PROJECT_ROOT, "experimental", "SO3lverED.jl"))
using .SO3lverED

function parse_options(args)
    options = Dict{String,String}()
    for argument in args
        startswith(argument, "--") || throw(ArgumentError(
            "Unknown positional argument '$argument'; use --key=value options",
        ))
        parts = split(argument[3:end], "="; limit=2)
        options[parts[1]] = length(parts) == 2 ? parts[2] : "true"
    end
    return options
end

option(options, name, default) = get(options, name, string(default))
option_bool(options, name, default=false) =
    lowercase(option(options, name, default)) in ("1", "true", "yes", "on")
raw_option_list(options, name, default) =
    filter(!isempty, strip.(split(option(options, name, default), ',')))
option_list(::Type{String}, options, name, default) =
    raw_option_list(options, name, default)
option_list(::Type{T}, options, name, default) where T =
    parse.(T, raw_option_list(options, name, default))

function couplings(options)
    defaults = Couplings()
    return Couplings(; (
        name => parse(Float64, option(options, String(name), getfield(defaults, name)))
        for name in fieldnames(Couplings)
    )...)
end

block_name(representation, ell) = "$(representation)_L$(ell)"

function sha256_file(path)
    return open(path) do io
        bytes2hex(sha256(io))
    end
end

function source_revision()
    project = TOML.parsefile(joinpath(PROJECT_ROOT, "Project.toml"))
    source = project["sources"]["FuzzifiED"]
    return get(source, "rev", "unknown")
end

function score_dict(score)
    return Dict{String,Any}(
        "q" => score.q,
        "factor" => score.factor,
        "delta_s" => score.delta_s,
        "delta_o" => score.delta_o,
        "ground_energy" => score.ground_energy,
        "ground_representation" => String(score.ground_representation),
        "ground_ell" => score.ground_ell,
        "ground_is_singlet_l0" => score.ground_is_singlet_l0,
        "labels" => score.labels,
        "raw_gaps" => score.raw_gaps,
        "target_gaps" => score.target_gaps,
        "scaled_gaps" => score.scaled_gaps,
        "residuals" => score.residuals,
    )
end

options = parse_options(ARGS)
nm = parse(Int, option(options, "nm", 7))
k = parse(Int, option(options, "k", 4))
tol = parse(Float64, option(options, "tol", 1e-8))
ncv = parse(Int, option(options, "ncv", max(2k, k + 10)))
warm_start = option_bool(options, "warm-start", true)
run_solver = option_bool(options, "solve", true)
base = couplings(options)
mus = option_list(Float64, options, "mus", string(base.mu))
ells = option_list(Int, options, "ells", "0,1,2")
representations = Symbol.(lowercase.(option_list(String, options, "representations", "singlet,adjoint")))
all(rep -> rep in (:singlet, :adjoint), representations) || throw(ArgumentError(
    "representations must contain only singlet and/or adjoint",
))
all(ell -> ell >= 0, ells) || throw(ArgumentError("ells must be non-negative"))
k > 0 || throw(ArgumentError("k must be positive"))
isempty(mus) && throw(ArgumentError("mus must not be empty"))
Random.seed!(parse(Int, option(options, "seed", 20260919)))

println("SO(3)lver CFT workload: N=$nm representations=$(join(representations, ',')) " *
        "L=$(join(ells, ',')) k=$k points=$(length(mus))")
println("Julia $(VERSION), threads=$(Threads.nthreads()), FuzzifiED $(Base.pkgversion(FuzzifiED))")
flush(stdout)

workspaces = Dict{Symbol,SO3Workspace}()
hamiltonians = Dict{Tuple{Symbol,Int},SO3Hamiltonian}()
workspace_seconds = Dict{String,Float64}()
operator_seconds = Dict{String,Float64}()
dimensions = Dict{String,Int}()
segment_dimensions = Dict{String,Dict{String,Int}}()
shared_heavy_space = Ref{Any}(nothing)

for representation in representations
    println("building_workspace representation=$representation")
    flush(stdout)
    started = time()
    model = build_so3_model(nm1=nm, representation=representation)
    workspace = build_workspace(
        model; heavy_space=shared_heavy_space[], disp_std=true,
    )
    shared_heavy_space[] = workspace.heavy_space
    elapsed = time() - started
    workspaces[representation] = workspace
    workspace_seconds[String(representation)] = elapsed
    segment_dimensions[String(representation)] = Dict(
        "light" => workspace.light_space.ptr_st[end][end] - 1,
        "heavy" => workspace.heavy_space.ptr_st[end][end] - 1,
    )
    println("workspace_seconds representation=$representation value=$elapsed")
    flush(stdout)

    for ell in ells
        started = time()
        hamiltonian = build_hamiltonian(workspace, ell, base; disp_std=true)
        elapsed = time() - started
        key = (representation, ell)
        name = block_name(representation, ell)
        hamiltonians[key] = hamiltonian
        operator_seconds[name] = elapsed
        dimensions[name] = sector_dimension(hamiltonian)
        println("operator block=$name dimension=$(dimensions[name]) seconds=$elapsed")
        flush(stdout)
    end
end

points = Dict{String,Any}[]
warm_vectors = Dict{Tuple{Symbol,Int},Vector{Float64}}()
if run_solver
    for (point_index, mu) in enumerate(mus)
        point_couplings = MottJainED.with_coupling(base, :mu, mu)
        block_energies = Dict{Tuple{Symbol,Int},Vector{Float64}}()
        block_results = Dict{String,Any}()
        point_started = time()
        for representation in representations, ell in ells
            key = (representation, ell)
            hamiltonian = hamiltonians[key]
            retune_started = time()
            retune!(hamiltonian, point_couplings)
            retune_seconds = time() - retune_started
            initvec = warm_start ? get(warm_vectors, key, nothing) : nothing
            solve_started = time()
            energies, vectors = solve(
                hamiltonian; k, tol, ncv, vectors=true, initvec,
                dense_cutoff=128, disp_std=true,
            )
            solve_seconds = time() - solve_started
            warm_start && !isempty(energies) && (warm_vectors[key] = copy(vectors[:, 1]))
            block_energies[key] = energies
            name = block_name(representation, ell)
            block_results[name] = Dict{String,Any}(
                "dimension" => sector_dimension(hamiltonian),
                "retune_seconds" => retune_seconds,
                "solve_seconds" => solve_seconds,
                "energies" => energies,
            )
            println("point=$point_index mu=$mu block=$name solve_seconds=$solve_seconds")
            flush(stdout)
        end
        point = Dict{String,Any}(
            "index" => point_index,
            "mu" => mu,
            "total_seconds" => time() - point_started,
            "blocks" => block_results,
        )
        if all(key -> haskey(block_energies, key), CFT_BLOCK_KEYS)
            score = score_cft_blocks(block_energies)
            point["score"] = score_dict(score)
            println("point=$point_index q=$(score.q) factor=$(score.factor) " *
                    "delta_s=$(score.delta_s) delta_o=$(score.delta_o)")
            flush(stdout)
        end
        push!(points, point)
    end
end

default_directory = joinpath(PROJECT_ROOT, "output", "so3lver", "workloads")
default_name = "n$(nm)_cft_workload_$(Dates.format(now(), "yyyymmdd_HHMMSS")).toml"
output = normpath(option(options, "output", joinpath(default_directory, default_name)))
mkpath(dirname(output))
fuzzified_root = normpath(joinpath(dirname(pathof(FuzzifiED)), ".."))
result = Dict{String,Any}(
    "completed_at" => string(now()),
    "nm1" => nm,
    "k" => k,
    "tol" => tol,
    "ncv" => ncv,
    "warm_start" => warm_start,
    "solved" => run_solver,
    "mus" => mus,
    "ells" => ells,
    "representations" => String.(representations),
    "workspace_seconds" => workspace_seconds,
    "operator_seconds" => operator_seconds,
    "dimensions" => dimensions,
    "segment_dimensions" => segment_dimensions,
    "shared_heavy_segment" => length(representations) > 1,
    "points" => points,
    "julia_version" => string(VERSION),
    "julia_threads" => Threads.nthreads(),
    "fuzzified_version" => string(Base.pkgversion(FuzzifiED)),
    "fuzzified_expected_revision" => source_revision(),
    "fuzzified_git_revision" => MottJainED.git_revision(fuzzified_root),
    "project_git_revision" => MottJainED.git_revision(PROJECT_ROOT),
    "so3lver_source_sha256" => sha256_file(joinpath(PROJECT_ROOT, "experimental", "SO3lverED.jl")),
    "driver_source_sha256" => sha256_file(@__FILE__),
    "couplings" => Dict(String(name) => getfield(base, name) for name in fieldnames(Couplings)),
)
open(output, "w") do io
    TOML.print(io, result; sorted=true)
end
println("workload_result=$output")
