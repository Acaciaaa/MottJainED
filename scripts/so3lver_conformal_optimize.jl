#!/usr/bin/env julia

import Pkg

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(PROJECT_ROOT; io=devnull)

using DataFrames
using Dates
using FuzzifiED
using LinearAlgebra
using MottJainED
using Optim
using Random
using SHA
using TOML

include(joinpath(PROJECT_ROOT, "experimental", "SO3lverED.jl"))
using .SO3lverED
include(joinpath(PROJECT_ROOT, "experimental", "SO3ParameterSearch.jl"))
using .SO3ParameterSearch

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

sha256_file(path) = open(path) do io
    bytes2hex(sha256(io))
end

function dict_couplings(values)
    defaults = Couplings()
    return Couplings(; (
        name => Float64(get(values, String(name), getfield(defaults, name)))
        for name in fieldnames(Couplings)
    )...)
end

options = parse_options(ARGS)
config_path = abspath(get(
    options, "config",
    joinpath(
        PROJECT_ROOT, "config", "so3lver",
        "n6_projected_conformal_optimization.toml",
    ),
))
isfile(config_path) || throw(ArgumentError("configuration not found: $config_path"))
config = TOML.parsefile(config_path)
run_config = config["run"]
algebra_config = config["algebra"]
search_config = config["search"]

nm = Int(run_config["nm"])
allow_test_size = lowercase(get(options, "allow-test-size", "false")) == "true"
(nm == 6 || allow_test_size) || error(
    "production conformal optimization is restricted to N=6; " *
    "pass --allow-test-size=true only for a smoke test",
)
tol = Float64(run_config["tol"])
ncv = Int(run_config["ncv"])
seed = Int(run_config["seed"])
heavy_space_mode = Symbol(lowercase(String(run_config["heavy_space_mode"])))
heavy_space_mode == :laughlin13 || error(
    "conformal optimization requires heavy_space_mode=\"laughlin13\"",
)
base = dict_couplings(config["base"])
iszero(base.V0) || error("V0 must be fixed to zero after Laughlin-1/3 projection")

fit_labels = Symbol.(algebra_config["fit_primaries"])
training_labels = Symbol.(algebra_config["training_primaries"])
holdout_labels = Symbol.(algebra_config["holdout_primaries"])
isempty(intersect(Set(training_labels), Set(holdout_labels))) || error(
    "training and holdout primary labels must be disjoint",
)
Set(fit_labels) == Set(training_labels) || error(
    "this pilot requires generator fit primaries to equal the outer training primaries",
)
factor_bounds = (
    Float64(algebra_config["factor_lower"]),
    Float64(algebra_config["factor_upper"]),
)
block_counts = Dict{Tuple{Symbol,Int},Int}(
    (:singlet, 0) => Int(algebra_config["singlet_l0_count"]),
    (:singlet, 2) => Int(algebra_config["singlet_l2_count"]),
    (:adjoint, 0) => Int(algebra_config["adjoint_l0_count"]),
    (:adjoint, 1) => Int(algebra_config["adjoint_l1_count"]),
)
term_weights = Dict(
    term => Float64(get(algebra_config["weights"], String(term), 0.0))
    for term in CONFORMAL_OBJECTIVE_TERMS
)
holdout_weights = copy(term_weights)
holdout_weights[:vacuum] = 0.0
worst_weight = Float64(algebra_config["worst_weight"])

parameters = Symbol.(search_config["parameters"])
length(unique(parameters)) == length(parameters) || error(
    "search parameters must not contain duplicates",
)
all(parameter -> parameter in fieldnames(Couplings), parameters) || error(
    "search contains an unknown Hamiltonian parameter",
)
:mu in parameters || error("mu (the finite-size muc) must be a search parameter")
:U0 in parameters && error("U0 is linked to Uf and cannot be an independent axis")
:V0 in parameters && error("V0 is redundant and fixed to zero in the projected model")
u0_over_uf = Float64(search_config["u0_over_uf"])
isapprox(base.U0, u0_over_uf * base.Uf; atol=1e-12, rtol=0) || error(
    "base point must obey U0 = u0_over_uf * Uf",
)
center = Float64[getfield(base, parameter) for parameter in parameters]
half_widths = Float64[
    search_config["initial_half_widths"][String(parameter)]
    for parameter in parameters
]
lower = Float64[
    search_config["bounds"][String(parameter)][1] for parameter in parameters
]
upper = Float64[
    search_config["bounds"][String(parameter)][2] for parameter in parameters
]
all(lower .< upper) || error("all search bounds must have positive width")
all((lower .<= center) .& (center .<= upper)) || error(
    "base point lies outside the search bounds",
)
all(half_widths .> 0) || error("initial half-widths must be positive")
minimum_identity_overlap = Float64(search_config["minimum_identity_overlap"])
0 <= minimum_identity_overlap <= 1 || error(
    "minimum_identity_overlap must lie in [0,1]",
)
penalty = Float64(search_config["penalty"])
initial_samples = parse(Int, get(
    options, "initial-samples", string(Int(search_config["initial_samples"])),
))
local_iterations = parse(Int, get(
    options, "local-iterations", string(Int(search_config["local_iterations"])),
))
simplex_step = Float64(search_config["simplex_step"])
initial_samples >= 0 || error("initial_samples must be non-negative")
local_iterations >= 0 || error("local_iterations must be non-negative")

output = abspath(get(
    options, "output", joinpath(PROJECT_ROOT, String(run_config["output"])),
))
validate_only = lowercase(get(options, "validate-only", "false")) == "true"
signature = MottJainED.stable_id(
    "projected-conformal-hamiltonian-optimization-v1",
    nm, tol, ncv, seed, String(heavy_space_mode),
    parameters, center, lower, upper, half_widths,
    fit_labels, training_labels, holdout_labels,
    factor_bounds,
    [term_weights[term] for term in CONFORMAL_OBJECTIVE_TERMS],
    worst_weight,
    sha256_file(config_path),
    sha256_file(joinpath(PROJECT_ROOT, "experimental", "SO3lverED.jl")),
    sha256_file(@__FILE__),
)

println("N=$nm projected conformal-algebra Hamiltonian optimization")
println("parameters=$(join(parameters, ',')); mu is optimized directly as muc")
println("training=$(join(training_labels, ',')) holdout=$(join(holdout_labels, ','))")
println("fixed V0=0; linked U0=$(u0_over_uf)*Uf")
println("config=$config_path")
println("output=$output")
println("Julia $(VERSION), threads=$(Threads.nthreads()), FuzzifiED $(Base.pkgversion(FuzzifiED))")
flush(stdout)
validate_only && exit(0)

mkpath(output)
trace_path = joinpath(output, "algebra_optimization_evaluations.csv")
constraint_path = joinpath(output, "algebra_optimization_constraints.csv")
identity_path = joinpath(output, "algebra_optimization_identity.csv")
for path in (trace_path, constraint_path, identity_path)
    isfile(path) && error(
        "output already contains $(basename(path)); choose a new directory",
    )
end
cp(config_path, joinpath(output, "optimization_config.toml"); force=true)
Random.seed!(seed)

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
problem = build_so3_conformal_problem(
    singlet_workspace, adjoint_workspace, base; disp_std=true,
)
println("problem_build_seconds=$(time() - started)")
flush(stdout)

function physical_couplings(values)
    couplings = base
    for (parameter, value) in zip(parameters, values)
        couplings = MottJainED.with_coupling(couplings, parameter, Float64(value))
    end
    couplings = MottJainED.with_coupling(couplings, :V0, 0.0)
    couplings = MottJainED.with_coupling(
        couplings, :U0, u0_over_uf * couplings.Uf,
    )
    MottJainED.validate(couplings)
    return couplings
end

function point_key(values)
    return MottJainED.stable_id(signature, round.(Float64.(values); digits=11))
end

reference_states = Dict{Tuple{Symbol,Int,Int},Vector{Float64}}()
cache = Dict{String,NamedTuple}()
evaluation = Ref(0)

function identity_assessment(result)
    rows = NamedTuple[]
    for spec in DEFAULT_CONFORMAL_PRIMARY_SPECS
        key = (spec.representation, spec.ell, spec.rank)
        source = view(result.states[(spec.representation, spec.ell)], :, spec.rank)
        overlap = isempty(reference_states) ? 1.0 :
            abs2(dot(reference_states[key], source))
        push!(rows, (
            label=spec.label,
            representation=spec.representation,
            ell=spec.ell,
            rank=spec.rank,
            overlap=overlap,
        ))
    end
    return (minimum=minimum(getproperty.(rows, :overlap)), rows=rows)
end

function record_constraints!(id, source, score, family)
    for row in score.rows
        MottJainED.append_csv(constraint_path, (
            signature=signature,
            evaluation=evaluation[],
            point_id=id,
            source=String(source),
            family=String(family),
            label=String(row.label),
            term=String(row.term),
            value=row.value,
            weight=row.weight,
            weighted_value=row.weighted_value,
        ))
    end
end

function evaluate_point(values; source="search", force=false)
    physical = Float64.(values)
    id = point_key(physical)
    !force && haskey(cache, id) && return cache[id]
    evaluation[] += 1
    point_started = time()
    couplings = physical_couplings(physical)
    try
        result = analyze_so3_conformal_algebra(
            singlet_workspace, adjoint_workspace, couplings;
            block_counts,
            fit_primary_labels=fit_labels,
            factor_bounds,
            vacuum_weight=Float64(algebra_config["vacuum_weight"]),
            dilatation_weight=Float64(algebra_config["dilatation_weight"]),
            descendant_weight=Float64(algebra_config["descendant_weight"]),
            descendant_state_count=Int(algebra_config["descendant_state_count"]),
            eig_tol=tol,
            ncv,
            prepared_problem=problem,
            disp_std=false,
        )
        training = score_so3_conformal_algebra(
            result; labels=training_labels, term_weights,
            worst_weight,
        )
        holdout = score_so3_conformal_algebra(
            result; labels=holdout_labels, term_weights=holdout_weights,
            worst_weight=0.0,
        )
        identity = identity_assessment(result)
        overlap = identity.minimum
        valid = result.fit.commutator_normalization_valid &&
                !result.fit.factor_at_boundary &&
                overlap >= minimum_identity_overlap
        reason = !result.fit.commutator_normalization_valid ?
            "commutator_normalization" :
            result.fit.factor_at_boundary ? "factor_boundary" :
            overlap < minimum_identity_overlap ? "state_identity" : "ok"
        objective = valid ? training.objective : penalty + training.objective
        summary = (
            id=id,
            values=physical,
            couplings=couplings,
            objective=objective,
            training=training,
            holdout=holdout,
            valid=valid,
            reason=reason,
            minimum_identity_overlap=overlap,
            factor=result.fit.factor,
            algebra_fit_loss=result.fit.value,
            identity=identity,
            result=result,
        )
        record_constraints!(id, source, training, :training)
        record_constraints!(id, source, holdout, :holdout)
        for row in identity.rows
            MottJainED.append_csv(identity_path, (
                signature=signature,
                evaluation=evaluation[],
                point_id=id,
                source=String(source),
                label=String(row.label),
                representation=String(row.representation),
                ell=row.ell,
                rank=row.rank,
                overlap=row.overlap,
                passed=row.overlap >= minimum_identity_overlap,
            ))
        end
        MottJainED.append_csv(trace_path, (
            signature=signature,
            evaluation=evaluation[],
            point_id=id,
            source=String(source),
            timestamp=string(now()),
            valid=valid,
            reason=reason,
            objective=objective,
            training_mean=training.mean,
            training_worst=training.worst,
            holdout_mean=holdout.mean,
            holdout_worst=holdout.worst,
            factor=result.fit.factor,
            factor_at_boundary=result.fit.factor_at_boundary,
            algebra_fit_loss=result.fit.value,
            minimum_identity_overlap=overlap,
            seconds=time() - point_started,
            Uf=couplings.Uf,
            U0=couplings.U0,
            Uf0=couplings.Uf0,
            Vf=couplings.Vf,
            Vf0=couplings.Vf0,
            V0=couplings.V0,
            t=couplings.t,
            muc=couplings.mu,
        ))
        compact = merge(summary, (result=nothing,))
        cache[id] = compact
        println(
            "evaluation=$(evaluation[]) source=$source objective=$objective " *
            "holdout=$(holdout.mean) valid=$valid muc=$(couplings.mu)",
        )
        flush(stdout)
        return force ? summary : compact
    catch err
        reason = replace(sprint(showerror, err), '\n' => ' ')
        summary = (
            id=id,
            values=physical,
            couplings=couplings,
            objective=penalty,
            training=nothing,
            holdout=nothing,
            valid=false,
            reason=reason,
            minimum_identity_overlap=0.0,
            factor=NaN,
            algebra_fit_loss=Inf,
            identity=nothing,
            result=nothing,
        )
        MottJainED.append_csv(trace_path, (
            signature=signature,
            evaluation=evaluation[],
            point_id=id,
            source=String(source),
            timestamp=string(now()),
            valid=false,
            reason=reason,
            objective=penalty,
            training_mean=Inf,
            training_worst=Inf,
            holdout_mean=Inf,
            holdout_worst=Inf,
            factor=NaN,
            factor_at_boundary=true,
            algebra_fit_loss=Inf,
            minimum_identity_overlap=0.0,
            seconds=time() - point_started,
            Uf=couplings.Uf,
            U0=couplings.U0,
            Uf0=couplings.Uf0,
            Vf=couplings.Vf,
            Vf0=couplings.Vf0,
            V0=couplings.V0,
            t=couplings.t,
            muc=couplings.mu,
        ))
        cache[id] = summary
        @error "conformal Hamiltonian point failed" source couplings exception=(err, catch_backtrace())
        return summary
    end
end

anchor = evaluate_point(center; source="anchor", force=true)
anchor.valid || error("the conformal optimization anchor is invalid: $(anchor.reason)")
for spec in DEFAULT_CONFORMAL_PRIMARY_SPECS
    key = (spec.representation, spec.ell, spec.rank)
    reference_states[key] = copy(view(
        anchor.result.states[(spec.representation, spec.ell)], :, spec.rank,
    ))
end
cache[anchor.id] = merge(anchor, (result=nothing,))

rng = MersenneTwister(seed)
sample_lower = max.(lower, center .- half_widths)
sample_upper = min.(upper, center .+ half_widths)
samples = initial_samples == 0 ? zeros(Float64, 0, length(parameters)) :
    latin_hypercube_points(initial_samples, sample_lower, sample_upper, rng)
for row in axes(samples, 1)
    evaluate_point(vec(samples[row, :]); source="latin_hypercube_$row")
end

valid_initial = filter(value -> value.valid, collect(values(cache)))
isempty(valid_initial) && error("initial conformal search found no valid point")
start = first(sort!(valid_initial; by=value -> value.objective))
local_result = nothing
if local_iterations > 0
    start_search = (start.values .- center) ./ half_widths
    function local_objective(search_values)
        physical = center .+ Float64.(search_values) .* half_widths
        if all((lower .<= physical) .& (physical .<= upper))
            return evaluate_point(physical; source="nelder_mead").objective
        end
        below = max.((lower .- physical) ./ half_widths, 0.0)
        above = max.((physical .- upper) ./ half_widths, 0.0)
        return penalty + sum(abs2, below .+ above)
    end
    local_result = optimize(
        local_objective,
        start_search,
        NelderMead(initial_simplex=Optim.AffineSimplexer(a=simplex_step, b=0.0)),
        Optim.Options(
            iterations=local_iterations,
            f_reltol=1.0e-4,
            show_trace=false,
            store_trace=false,
        ),
    )
    local_values = center .+ Optim.minimizer(local_result) .* half_widths
    all((lower .<= local_values) .& (local_values .<= upper)) &&
        evaluate_point(local_values; source="nelder_mead_minimum")
end

valid_results = sort!(
    filter(value -> value.valid && isfinite(value.objective), collect(values(cache)));
    by=value -> value.objective,
)
isempty(valid_results) && error("conformal optimization produced no valid point")
best_cached = first(valid_results)
best = evaluate_point(best_cached.values; source="best_recheck", force=true)
best.valid || error("best conformal point failed recheck: $(best.reason)")

primary_rows = DataFrame(best.result.primary_rows)
mixed_rows = DataFrame(best.result.mixed_commutator_rows)
component_rows = DataFrame(best.result.mixed_commutator_pair_rows)
MottJainED.atomic_csv(joinpath(output, "best_primary_residuals.csv"), primary_rows)
MottJainED.atomic_csv(joinpath(output, "best_commutator_residuals.csv"), mixed_rows)
MottJainED.atomic_csv(joinpath(output, "best_commutator_components.csv"), component_rows)
MottJainED.atomic_csv(
    joinpath(output, "best_training_constraints.csv"), DataFrame(best.training.rows),
)
MottJainED.atomic_csv(
    joinpath(output, "best_holdout_constraints.csv"), DataFrame(best.holdout.rows),
)
MottJainED.atomic_csv(
    joinpath(output, "best_identity_overlaps.csv"), DataFrame(best.identity.rows),
)

best_dict = Dict{String,Any}(
    "completed_at" => string(now()),
    "signature" => signature,
    "nm1" => nm,
    "heavy_space_mode" => String(heavy_space_mode),
    "parameters" => String.(parameters),
    "training_primaries" => String.(training_labels),
    "holdout_primaries" => String.(holdout_labels),
    "objective" => best.objective,
    "training_mean" => best.training.mean,
    "training_worst" => best.training.worst,
    "holdout_mean" => best.holdout.mean,
    "holdout_worst" => best.holdout.worst,
    "factor" => best.factor,
    "algebra_fit_loss" => best.algebra_fit_loss,
    "minimum_identity_overlap" => best.minimum_identity_overlap,
    "anchor_objective" => anchor.objective,
    "anchor_training_mean" => anchor.training.mean,
    "anchor_holdout_mean" => anchor.holdout.mean,
    "evaluation_count" => evaluation[],
    "initial_samples" => initial_samples,
    "local_iterations" => local_iterations,
    "optimizer_converged" => isnothing(local_result) ? false : Optim.converged(local_result),
    "couplings" => Dict(
        String(name) => getfield(best.couplings, name) for name in fieldnames(Couplings)
    ),
    "term_weights" => Dict(String(key) => value for (key, value) in term_weights),
    "project_git_revision" => MottJainED.git_revision(PROJECT_ROOT),
    "fuzzified_version" => string(Base.pkgversion(FuzzifiED)),
)
MottJainED.atomic_toml(joinpath(output, "best.toml"), best_dict)

println("anchor_objective=$(anchor.objective) anchor_holdout=$(anchor.holdout.mean)")
println("best_objective=$(best.objective) best_holdout=$(best.holdout.mean)")
println(
    "best_parameters=Uf=$(best.couplings.Uf),U0=$(best.couplings.U0)," *
    "Uf0=$(best.couplings.Uf0),Vf0=$(best.couplings.Vf0)," *
    "V0=$(best.couplings.V0),muc=$(best.couplings.mu)",
)
println("optimization_result=$(joinpath(output, "best.toml"))")
