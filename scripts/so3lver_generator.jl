#!/usr/bin/env julia

import Pkg

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
Base.active_project() == joinpath(PROJECT_ROOT, "Project.toml") ||
    Pkg.activate(PROJECT_ROOT; io=devnull)

using DataFrames
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
        length(parts) == 2 || throw(ArgumentError(
            "Option '$argument' needs a value",
        ))
        options[parts[1]] = parts[2]
    end
    return options
end

function load_couplings(config)
    defaults = Couplings()
    return Couplings(; (
        name => Float64(get(config, String(name), getfield(defaults, name)))
        for name in fieldnames(Couplings)
    )...)
end

function state_label(ell, rank)
    labels = Dict(
        (0, 1) => "G", (0, 2) => "S", (0, 3) => "boxS",
        (1, 1) => "dS", (2, 1) => "T", (2, 2) => "ddS",
    )
    return get(labels, (ell, rank), "L$(ell)_rank$(rank)")
end

options = parse_options(ARGS)
haskey(options, "config") || throw(ArgumentError(
    "Pass an explicit --config=PATH; generator points are never inferred from my_run.toml",
))
config_path = abspath(options["config"])
isfile(config_path) || throw(ArgumentError("configuration not found: $config_path"))
config = TOML.parsefile(config_path)
run = config["run"]
point = load_couplings(config["point"])
MottJainED.validate(point)

nm = Int(run["nm"])
l0_count = Int(get(run, "l0_count", 6))
l2_count = Int(get(run, "l2_count", 6))
tol = Float64(get(run, "tol", 1.0e-9))
ncv = Int(get(run, "ncv", max(18, 2 * max(l0_count, l2_count))))
seed = Int(get(run, "seed", 20260921))
Random.seed!(seed)
output = abspath(get(
    options, "output", joinpath(PROJECT_ROOT, String(run["output"])),
))
mkpath(output)

println("native SO(3)lver generator: N=$nm L0_count=$l0_count L2_count=$l2_count")
println("config=$config_path")
println("output=$output")
println("Julia $(VERSION), threads=$(Threads.nthreads()), FuzzifiED $(Base.pkgversion(FuzzifiED))")
flush(stdout)

started = time()
model = build_so3_model(nm1=nm, representation=:singlet)
workspace = build_workspace(model; disp_std=true)
workspace_seconds = time() - started
println("workspace_seconds=$workspace_seconds")
flush(stdout)

started = time()
result = analyze_scalar_generator(
    workspace, point;
    l0_count, l2_count, eig_tol=tol, ncv, disp_std=true,
)
analysis_seconds = time() - started

MottJainED.atomic_csv(
    joinpath(output, "generator_coefficients.csv"),
    DataFrame(
        name=collect(String.(result.fit.names)),
        coefficient=Float64.(result.fit.coefficients),
    ),
)

overlap_rows = NamedTuple[]
for (ell, energies, overlaps) in (
    (0, result.energies.l0, result.l0_overlap.values),
    (2, result.energies.l2, result.l2_overlap.values),
)
    for rank in eachindex(energies)
        push!(overlap_rows, (
            ell=ell,
            rank=rank,
            label=state_label(ell, rank),
            energy=energies[rank],
            overlap=overlaps[rank],
        ))
    end
end
MottJainED.atomic_csv(
    joinpath(output, "generator_overlaps.csv"), DataFrame(overlap_rows),
)
MottJainED.atomic_csv(
    joinpath(output, "generator_summary.csv"),
    DataFrame([(
        nm1=nm,
        fidelity=result.fit.fidelity,
        numerical_rank=result.fit.numerical_rank,
        candidate_count=length(result.fit.names),
        l0_resolved_overlap=result.l0_overlap.total,
        l2_resolved_overlap=result.l2_overlap.total,
        workspace_seconds=workspace_seconds,
        analysis_seconds=analysis_seconds,
    )]),
)

metadata = Dict{String,Any}(
    "completed_at" => string(now()),
    "config_path" => config_path,
    "nm1" => nm,
    "l0_count" => l0_count,
    "l2_count" => l2_count,
    "tol" => tol,
    "ncv" => ncv,
    "seed" => seed,
    "fidelity" => result.fit.fidelity,
    "numerical_rank" => result.fit.numerical_rank,
    "candidate_names" => collect(String.(result.fit.names)),
    "singular_values" => Float64.(result.fit.singular_values),
    "workspace_seconds" => workspace_seconds,
    "analysis_seconds" => analysis_seconds,
    "julia_version" => string(VERSION),
    "julia_threads" => Threads.nthreads(),
    "fuzzified_version" => string(Base.pkgversion(FuzzifiED)),
    "project_git_revision" => MottJainED.git_revision(PROJECT_ROOT),
    "couplings" => Dict(
        String(name) => getfield(point, name) for name in fieldnames(Couplings)
    ),
)
MottJainED.atomic_toml(joinpath(output, "generator_metadata.toml"), metadata)

println("fit_fidelity=$(result.fit.fidelity) numerical_rank=$(result.fit.numerical_rank)")
println("l0_resolved_overlap=$(result.l0_overlap.total)")
println("l2_resolved_overlap=$(result.l2_overlap.total)")
println("generator_result=$output")
