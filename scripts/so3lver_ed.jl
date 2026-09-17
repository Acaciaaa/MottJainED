#!/usr/bin/env julia

import Pkg

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(PROJECT_ROOT; io=devnull)

using Dates
using FuzzifiED
using MottJainED
using Random
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

function couplings(options)
    defaults = Couplings()
    return Couplings(; (
        name => parse(Float64, option(options, String(name), getfield(defaults, name)))
        for name in fieldnames(Couplings)
    )...)
end

function usage()
    println("""
Usage:
  julia --project=. scripts/so3lver_ed.jl benchmark \\
      --nm=7 --representation=adjoint --ell=2 --k=4 [--solve=false]

Optional Hamiltonian overrides use --Uf=..., --Uf0=..., --U0=..., --Vf=...,
--Vf0=..., --V0=..., --t=..., and --mu=....  Results are written below
output/so3lver/benchmarks unless --output=PATH is supplied.
""")
end

isempty(ARGS) && (usage(); exit(1))
command = lowercase(first(ARGS))
command == "benchmark" || (usage(); throw(ArgumentError("Unknown command '$command'")))
options = parse_options(ARGS[2:end])

nm = parse(Int, option(options, "nm", 7))
representation = Symbol(lowercase(option(options, "representation", "adjoint")))
ell = parse(Int, option(options, "ell", 2))
k = parse(Int, option(options, "k", 4))
tol = parse(Float64, option(options, "tol", 1e-8))
run_solver = option_bool(options, "solve", true)
params = couplings(options)
Random.seed!(parse(Int, option(options, "seed", 20260917)))

println("SO(3)lver benchmark: N=$nm representation=$representation L=$ell k=$k")
println("Julia $(VERSION), threads=$(Threads.nthreads()), FuzzifiED $(Base.pkgversion(FuzzifiED))")
flush(stdout)

started = time()
model = build_so3_model(nm1=nm, representation=representation)
workspace = build_workspace(model; disp_std=true)
workspace_elapsed = time() - started
println("workspace_seconds=$workspace_elapsed")
flush(stdout)

started = time()
hamiltonian = build_hamiltonian(workspace, ell, params; disp_std=true)
operator_elapsed = time() - started
dimension = sector_dimension(hamiltonian)
println("sector_dimension=$dimension operator_seconds=$operator_elapsed")
flush(stdout)

matvec_elapsed = NaN
steady_matvec_elapsed = NaN
solve_elapsed = NaN
energies = Float64[]
if dimension > 0
    vector = randn(dimension)
    started = time()
    hamiltonian.operator * vector
    matvec_elapsed = time() - started
    println("first_matvec_seconds=$matvec_elapsed")
    started = time()
    hamiltonian.operator * vector
    steady_matvec_elapsed = time() - started
    println("steady_matvec_seconds=$steady_matvec_elapsed")
    flush(stdout)
end
if run_solver && dimension > 0
    started = time()
    energies = solve(
        hamiltonian; k, tol, vectors=false, dense_cutoff=128, disp_std=true,
    )
    solve_elapsed = time() - started
    println("solve_seconds=$solve_elapsed energies=$(join(energies, ','))")
    flush(stdout)
end

default_directory = joinpath(PROJECT_ROOT, "output", "so3lver", "benchmarks")
default_name = "n$(nm)_$(representation)_l$(ell)_$(Dates.format(now(), "yyyymmdd_HHMMSS")).toml"
output = normpath(option(options, "output", joinpath(default_directory, default_name)))
mkpath(dirname(output))
result = Dict{String,Any}(
    "completed_at" => string(now()),
    "nm1" => nm,
    "representation" => String(representation),
    "ell" => ell,
    "k" => k,
    "tol" => tol,
    "solved" => run_solver,
    "sector_dimension" => dimension,
    "light_segment_dimension" => workspace.light_space.ptr_st[end][end] - 1,
    "heavy_segment_dimension" => workspace.heavy_space.ptr_st[end][end] - 1,
    "workspace_seconds" => workspace_elapsed,
    "operator_seconds" => operator_elapsed,
    "first_matvec_seconds" => matvec_elapsed,
    "steady_matvec_seconds" => steady_matvec_elapsed,
    "solve_seconds" => solve_elapsed,
    "energies" => energies,
    "julia_version" => string(VERSION),
    "julia_threads" => Threads.nthreads(),
    "fuzzified_version" => string(Base.pkgversion(FuzzifiED)),
    "fuzzified_git_revision" => MottJainED.git_revision(normpath(joinpath(dirname(pathof(FuzzifiED)), ".."))),
    "project_git_revision" => MottJainED.git_revision(PROJECT_ROOT),
    "couplings" => Dict(String(name) => getfield(params, name) for name in fieldnames(Couplings)),
)
open(output, "w") do io
    TOML.print(io, result; sorted=true)
end
println("benchmark_result=$output")
