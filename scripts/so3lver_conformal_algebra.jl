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
        length(parts) == 2 || throw(ArgumentError("Option '$argument' needs a value"))
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

options = parse_options(ARGS)
haskey(options, "config") || throw(ArgumentError(
    "Pass an explicit --config=PATH; algebra points are never inferred from my_run.toml",
))
config_path = abspath(options["config"])
isfile(config_path) || throw(ArgumentError("configuration not found: $config_path"))
config = TOML.parsefile(config_path)
run = config["run"]
algebra = config["algebra"]
point = load_couplings(config["point"])
MottJainED.validate(point)

nm = Int(run["nm"])
heavy_space_mode = Symbol(lowercase(String(run["heavy_space_mode"])))
heavy_space_mode == :laughlin13 || throw(ArgumentError(
    "this projected conformal-algebra pilot requires heavy_space_mode=\"laughlin13\"",
))
tol = Float64(get(run, "tol", 1.0e-8))
ncv = Int(get(run, "ncv", 18))
seed = Int(get(run, "seed", 20260924))
Random.seed!(seed)
output = abspath(get(
    options, "output", joinpath(PROJECT_ROOT, String(run["output"])),
))
mkpath(output)

block_counts = Dict{Tuple{Symbol,Int},Int}(
    (:singlet, 0) => Int(get(algebra, "singlet_l0_count", 5)),
    (:singlet, 2) => Int(get(algebra, "singlet_l2_count", 4)),
    (:adjoint, 0) => Int(get(algebra, "adjoint_l0_count", 4)),
    (:adjoint, 1) => Int(get(algebra, "adjoint_l1_count", 4)),
)
factor_bounds = (
    Float64(algebra["factor_lower"]),
    Float64(algebra["factor_upper"]),
)
fit_primaries = Symbol.(String.(algebra["fit_primaries"]))
vacuum_weight = Float64(get(algebra, "vacuum_weight", 1.0))
dilatation_weight = Float64(get(algebra, "dilatation_weight", 1.0))
descendant_weight = Float64(get(algebra, "descendant_weight", 1.0))
descendant_state_count = Int(get(algebra, "descendant_state_count", 6))
gram_rtol = Float64(get(algebra, "gram_rtol", 1.0e-10))

println("projected native SO(3)lver conformal-algebra pilot: N=$nm")
println("config=$config_path")
println("output=$output")
println("fit_primaries=$(join(String.(fit_primaries), ','))")
println("factor_bounds=$factor_bounds")
println("Julia $(VERSION), threads=$(Threads.nthreads()), FuzzifiED $(Base.pkgversion(FuzzifiED))")
flush(stdout)

started = time()
singlet_model = build_so3_model(nm1=nm, representation=:singlet)
singlet_workspace = build_workspace(
    singlet_model; heavy_space_mode, disp_std=true,
)
adjoint_model = build_so3_model(nm1=nm, representation=:adjoint)
adjoint_workspace = build_workspace(
    adjoint_model;
    heavy_space=singlet_workspace.heavy_space,
    heavy_space_mode,
    disp_std=true,
)
workspace_seconds = time() - started
println("workspace_seconds=$workspace_seconds")
flush(stdout)

started = time()
result = analyze_so3_conformal_algebra(
    singlet_workspace, adjoint_workspace, point;
    block_counts,
    fit_primary_labels=fit_primaries,
    factor_bounds,
    vacuum_weight,
    dilatation_weight,
    descendant_weight,
    descendant_state_count,
    eig_tol=tol,
    ncv,
    gram_rtol,
    disp_std=true,
)
analysis_seconds = time() - started

MottJainED.atomic_csv(
    joinpath(output, "algebra_generator_coefficients.csv"),
    DataFrame(
        name=collect(String.(result.generator_candidate_names)),
        coefficient=Float64.(result.fit.coefficients),
    ),
)
MottJainED.atomic_csv(
    joinpath(output, "algebra_primary_residuals.csv"),
    DataFrame(result.primary_rows),
)
MottJainED.atomic_csv(
    joinpath(output, "algebra_channel_residuals.csv"),
    DataFrame(result.channel_rows),
)

k2_rows = NamedTuple[]
k2_vector_rows = NamedTuple[]
for ((representation, ell), sector) in sort!(collect(result.k2); by=first)
    for mode in eachindex(sector.eigenvalues)
        push!(k2_rows, (
            representation=representation,
            ell=ell,
            mode=mode,
            k2=sector.eigenvalues[mode],
            lambda_norm2=sector.lambda_norm2[mode],
            k_fraction=sector.k_fraction[mode],
            energy_expectation=sector.energy_expectation[mode],
        ))
        for rank in eachindex(sector.source_energies)
            push!(k2_vector_rows, (
                representation=representation,
                ell=ell,
                mode=mode,
                energy_rank=rank,
                source_energy=sector.source_energies[rank],
                coefficient=sector.eigenvectors[rank, mode],
            ))
        end
    end
end
MottJainED.atomic_csv(joinpath(output, "algebra_k2_modes.csv"), DataFrame(k2_rows))
MottJainED.atomic_csv(
    joinpath(output, "algebra_k2_eigenvectors.csv"), DataFrame(k2_vector_rows),
)

maximum_k_fraction = maximum(getproperty.(result.primary_rows, :k_fraction))
maximum_p_dilatation_fraction = maximum(
    getproperty.(result.primary_rows, :p_dilatation_fraction),
)
maximum_p_low_energy_leakage_fraction = maximum(
    getproperty.(result.primary_rows, :p_low_energy_leakage_fraction),
)
scalar_commutator_fractional_residuals = [
    row.kp_commutator_fractional_residual for row in result.primary_rows
    if isfinite(row.kp_commutator_fractional_residual)
]
scalar_commutator_rms_fraction = isempty(scalar_commutator_fractional_residuals) ?
    NaN : sqrt(sum(abs2, scalar_commutator_fractional_residuals) /
               length(scalar_commutator_fractional_residuals))
summary = DataFrame([(
    nm1=nm,
    heavy_space_mode=String(result.heavy_space_mode),
    algebra_loss=result.fit.value,
    factor=result.fit.factor,
    factor_at_boundary=result.fit.factor_at_boundary,
    normalization_rank=result.fit.normalization_rank,
    vacuum_norm2=result.vacuum_norm2,
    maximum_primary_k_fraction=maximum_k_fraction,
    maximum_primary_p_dilatation_fraction=maximum_p_dilatation_fraction,
    maximum_primary_p_low_energy_leakage_fraction=
        maximum_p_low_energy_leakage_fraction,
    scalar_commutator_rms_fraction=scalar_commutator_rms_fraction,
    commutator_normalization_valid=result.fit.commutator_normalization_valid,
    workspace_seconds=workspace_seconds,
    analysis_seconds=analysis_seconds,
)])
MottJainED.atomic_csv(joinpath(output, "algebra_summary.csv"), summary)

metadata = Dict{String,Any}(
    "completed_at" => string(now()),
    "config_path" => config_path,
    "nm1" => nm,
    "heavy_space_mode" => String(result.heavy_space_mode),
    "fit_primary_labels" => String.(result.fit.fit_primary_labels),
    "factor" => result.fit.factor,
    "factor_bounds" => collect(result.fit.factor_bounds),
    "factor_at_boundary" => result.fit.factor_at_boundary,
    "factor_scan" => result.fit.factor_scan,
    "factor_scan_values" => result.fit.factor_scan_values,
    "algebra_loss" => result.fit.value,
    "vacuum_norm2" => result.vacuum_norm2,
    "vacuum_weight" => result.fit.vacuum_weight,
    "dilatation_weight" => result.fit.dilatation_weight,
    "descendant_weight" => result.fit.descendant_weight,
    "descendant_state_count" => descendant_state_count,
    "commutator_normalization_valid" => result.fit.commutator_normalization_valid,
    "commutator_normalization_scale" => result.fit.commutator_normalization_scale,
    "scalar_commutator_residuals" => result.fit.scalar_commutator_residuals,
    "normalization_rank" => result.fit.normalization_rank,
    "generalized_eigenvalues" => result.fit.spectrum,
    "normalization_eigenvalues" => result.fit.normalization_eigenvalues,
    "dimensions" => Dict(
        "$(key[1])_L$(key[2])" => value for (key, value) in result.dimensions
    ),
    "block_counts" => Dict(
        "$(key[1])_L$(key[2])" => value for (key, value) in block_counts
    ),
    "tol" => tol,
    "ncv" => ncv,
    "seed" => seed,
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
MottJainED.atomic_toml(joinpath(output, "algebra_metadata.toml"), metadata)

println("algebra_loss=$(result.fit.value)")
println("factor=$(result.fit.factor) factor_at_boundary=$(result.fit.factor_at_boundary)")
println("vacuum_norm2=$(result.vacuum_norm2)")
for row in result.primary_rows
    println(
        "primary=$(row.label) K2/Lambda2=$(row.k_fraction) " *
        "P_dilatation=$(row.p_dilatation_fraction) " *
        "P_leakage=$(row.p_low_energy_leakage_fraction) " *
        "KP_residual=$(row.kp_commutator_fractional_residual)",
    )
end
println("algebra_result=$output")
