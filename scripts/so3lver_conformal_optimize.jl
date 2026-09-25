#!/usr/bin/env julia

import Pkg

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(PROJECT_ROOT; io=devnull)

using CSV
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
    joinpath(PROJECT_ROOT, "config", "so3lver",
             "n6_projected_conformal_optimization.toml"),
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
    "generator fit primaries must equal the outer training primaries",
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
parameters == [:Uf, :Uf0, :Vf0, :mu] || error(
    "production conformal search parameters must be ordered Uf,Uf0,Vf0,mu",
)
u0_over_uf = Float64(search_config["u0_over_uf"])
isapprox(base.U0, u0_over_uf * base.Uf; atol=1e-12, rtol=0) || error(
    "base point must obey U0 = u0_over_uf * Uf",
)
center = Float64[getfield(base, parameter) for parameter in parameters]
half_widths = Float64[
    search_config["initial_half_widths"][String(parameter)] for parameter in parameters
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
local_starts = parse(Int, get(
    options, "local-starts", string(Int(search_config["local_starts"])),
))
local_iterations = parse(Int, get(
    options, "local-iterations", string(Int(search_config["local_iterations"])),
))
simplex_step = Float64(search_config["simplex_step"])
maximum_rounds = Int(search_config["maximum_rounds"])
expansion_factor = Float64(search_config["expansion_factor"])
boundary_fraction = Float64(search_config["boundary_fraction"])
minimum_start_distance = Float64(search_config["minimum_start_distance"])
top_candidate_count = Int(search_config["top_candidate_count"])
cold_recheck_tolerance = Float64(search_config["cold_recheck_tolerance"])
holdout_maximum_anchor_ratio = Float64(search_config["holdout_maximum_anchor_ratio"])
consensus_objective_tolerance = Float64(
    search_config["consensus_objective_tolerance"],
)
consensus_parameter_tolerances = Float64[
    search_config["consensus_parameter_tolerances"][String(parameter)]
    for parameter in parameters
]
robustness_steps = Float64[
    search_config["robustness_steps"][String(parameter)] for parameter in parameters
]
robustness_objective_tolerance = Float64(
    search_config["robustness_objective_tolerance"],
)
hard_boundary_fraction = Float64(search_config["hard_boundary_fraction"])

initial_samples >= 1 || error("initial_samples must be positive")
local_starts >= 1 || error("local_starts must be positive")
local_iterations >= 1 || error("local_iterations must be positive")
maximum_rounds >= 1 || error("maximum_rounds must be positive")
expansion_factor > 1 || error("expansion_factor must exceed one")
0 <= boundary_fraction < 0.5 || error("boundary_fraction must lie in [0,0.5)")
minimum_start_distance >= 0 || error("minimum_start_distance must be nonnegative")
all(consensus_parameter_tolerances .> 0) || error(
    "consensus parameter tolerances must be positive",
)
all(robustness_steps .> 0) || error("robustness steps must be positive")
cold_recheck_tolerance >= 0 || error("cold recheck tolerance must be nonnegative")
holdout_maximum_anchor_ratio > 0 || error(
    "holdout maximum anchor ratio must be positive",
)

output = abspath(get(
    options, "output", joinpath(PROJECT_ROOT, String(run_config["output"])),
))
validate_only = lowercase(get(options, "validate-only", "false")) == "true"
so3_path = joinpath(PROJECT_ROOT, "experimental", "SO3lverED.jl")
signature = MottJainED.stable_id(
    "projected-conformal-hamiltonian-optimization-v3-cg-norm",
    nm, tol, ncv, seed, String(heavy_space_mode), parameters, center,
    lower, upper, half_widths, initial_samples, local_starts, local_iterations,
    simplex_step, maximum_rounds, expansion_factor, boundary_fraction,
    minimum_start_distance, fit_labels, training_labels, holdout_labels,
    factor_bounds, [term_weights[term] for term in CONFORMAL_OBJECTIVE_TERMS],
    worst_weight, minimum_identity_overlap, cold_recheck_tolerance,
    holdout_maximum_anchor_ratio, consensus_parameter_tolerances,
    consensus_objective_tolerance, robustness_steps,
    robustness_objective_tolerance, hard_boundary_fraction,
    sha256_file(config_path), sha256_file(so3_path), sha256_file(@__FILE__),
)

trace_path = joinpath(output, "algebra_optimization_evaluations.csv")
constraint_path = joinpath(output, "algebra_optimization_constraints.csv")
identity_path = joinpath(output, "algebra_optimization_identity.csv")

function validate_trace_signature(path)
    isfile(path) || return nothing
    previous = CSV.read(path, DataFrame)
    "signature" in names(previous) || error(
        "$(basename(path)) has no signature; choose a new output directory",
    )
    all(String(value) == signature for value in previous.signature) || error(
        "$(basename(path)) belongs to another search; choose a new output directory",
    )
    return previous
end

previous_trace = validate_trace_signature(trace_path)
validate_trace_signature(constraint_path)
validate_trace_signature(identity_path)
resuming = !isnothing(previous_trace)
evaluation = Ref(resuming ? maximum(Int.(previous_trace.evaluation)) : 0)

println("N=$nm projected conformal-algebra Hamiltonian optimization")
println("parameters=$(join(parameters, ',')); mu is optimized directly as muc")
println("training=$(join(training_labels, ',')) holdout=$(join(holdout_labels, ','))")
println("fixed V0=0; linked U0=$(u0_over_uf)*Uf")
println(
    "seed=$seed initial_samples=$initial_samples local_starts=$local_starts " *
    "maximum_rounds=$maximum_rounds resume=$resuming",
)
println("config=$config_path")
println("output=$output")
println("Julia $(VERSION), threads=$(Threads.nthreads()), FuzzifiED $(Base.pkgversion(FuzzifiED))")
flush(stdout)
if validate_only
    println("VALIDATION_OK: configuration, dependencies, and source signatures loaded")
    exit(0)
end

mkpath(output)
cp(config_path, joinpath(output, "optimization_config.toml"); force=true)

started = time()
singlet_model = build_so3_model(nm1=nm, representation=:singlet)
singlet_workspace = build_workspace(singlet_model; heavy_space_mode, disp_std=true)
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

point_key(values) = MottJainED.stable_id(
    signature, round.(Float64.(values); digits=11),
)

reference_states = Dict{Tuple{Symbol,Int,Int},Vector{Float64}}()
cache = Dict{String,NamedTuple}()

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

function evaluate_point(values; source="search", force=false, record=true)
    physical = Float64.(values)
    id = point_key(physical)
    !force && haskey(cache, id) && return cache[id]
    record && (evaluation[] += 1)
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
            result; labels=training_labels, term_weights, worst_weight,
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
            id=id, source=String(source), values=physical, couplings=couplings,
            objective=objective,
            training=training, holdout=holdout, valid=valid, reason=reason,
            minimum_identity_overlap=overlap, factor=result.fit.factor,
            algebra_fit_loss=result.fit.value, identity=identity, result=result,
        )
        if record
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
        end
        compact = merge(summary, (result=nothing,))
        cache[id] = compact
        if record
            println(
                "evaluation=$(evaluation[]) source=$source objective=$objective " *
                "holdout=$(holdout.mean) valid=$valid muc=$(couplings.mu)",
            )
            flush(stdout)
        end
        return force ? summary : compact
    catch err
        reason = replace(sprint(showerror, err), '\n' => ' ')
        summary = (
            id=id, source=String(source), values=physical, couplings=couplings,
            objective=penalty,
            training=nothing, holdout=nothing, valid=false, reason=reason,
            minimum_identity_overlap=0.0, factor=NaN, algebra_fit_loss=Inf,
            identity=nothing, result=nothing,
        )
        if record
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
            @error "conformal Hamiltonian point failed" source couplings exception=(err, catch_backtrace())
        end
        cache[id] = summary
        return summary
    end
end

# Always solve the anchor in the current process so identity tracking uses fresh
# eigenvectors.  A resumed job does not append that administrative solve.
anchor_full = evaluate_point(center; source="anchor", force=true, record=!resuming)
anchor_full.valid || error("the conformal optimization anchor is invalid: $(anchor_full.reason)")
for spec in DEFAULT_CONFORMAL_PRIMARY_SPECS
    key = (spec.representation, spec.ell, spec.rank)
    reference_states[key] = copy(view(
        anchor_full.result.states[(spec.representation, spec.ell)], :, spec.rank,
    ))
end
anchor = merge(anchor_full, (result=nothing,))
cache[anchor.id] = anchor

# Final audits are excluded from the resume cache so they cannot silently alter
# the deterministic search trajectory.
if resuming
    excluded_final_audits = Ref(0)
    for row in eachrow(previous_trace)
        source = String(row.source)
        if source in ("best_recheck", "best_cold_recheck") ||
                startswith(source, "consensus_") ||
                startswith(source, "robustness_")
            excluded_final_audits[] += 1
            continue
        end
        values = Float64[row.Uf, row.Uf0, row.Vf0, row.muc]
        id = point_key(values)
        cache[id] = (
            id=id,
            source=source,
            values=values,
            couplings=physical_couplings(values),
            objective=Float64(row.objective),
            training=(mean=Float64(row.training_mean), worst=Float64(row.training_worst)),
            holdout=(mean=Float64(row.holdout_mean), worst=Float64(row.holdout_worst)),
            valid=Bool(row.valid),
            reason=String(row.reason),
            minimum_identity_overlap=Float64(row.minimum_identity_overlap),
            factor=Float64(row.factor),
            algebra_fit_loss=Float64(row.algebra_fit_loss),
            identity=nothing,
            result=nothing,
        )
    end
    cache[anchor.id] = anchor
    println(
        "resumed_search_points=$(length(cache)) " *
        "excluded_final_audits=$(excluded_final_audits[])",
    )
    flush(stdout)
end
anchor_full = nothing
GC.gc()

function diverse_starts(candidates, count, region_lower, region_upper)
    valid = sort(
        filter(row -> row.valid && isfinite(row.objective), candidates);
        by=row -> row.objective,
    )
    selected = Any[]
    width = region_upper .- region_lower
    for candidate in valid
        far_enough = all(
            norm((candidate.values .- existing.values) ./ width) >=
            minimum_start_distance for existing in selected
        )
        (isempty(selected) || far_enough) && push!(selected, candidate)
        length(selected) >= count && break
    end
    return selected
end

function run_search()
    rng = MersenneTwister(seed)
    region_center = copy(center)
    region_half_widths = copy(half_widths)
    local_run_rows = Dict{String,Any}[]
    local_endpoints = Dict{Tuple{Int,Int},Any}()
    rounds_completed = 0

    for round_index in 1:maximum_rounds
        rounds_completed = round_index
        region_lower = max.(lower, region_center .- region_half_widths)
        region_upper = min.(upper, region_center .+ region_half_widths)
        println(
            "round=$round_index lower=$(join(region_lower, ',')) " *
            "upper=$(join(region_upper,','))",
        )
        flush(stdout)

        candidates = Any[evaluate_point(
            region_center; source="round_$(round_index)_center",
        )]
        samples = latin_hypercube_points(
            initial_samples, region_lower, region_upper, rng,
        )
        for row_index in axes(samples, 1)
            push!(candidates, evaluate_point(
                vec(samples[row_index, :]);
                source="round_$(round_index)_lhs_$(row_index)",
            ))
        end
        current_local_prefix = "round_$(round_index)_local_"
        for candidate in values(cache)
            # On a resumed round, do not promote that same round's previous
            # simplex endpoint into a new starting point.  The original
            # pre-local candidate set is reconstructed exactly, while points
            # from earlier trust-region rounds remain eligible.
            startswith(candidate.source, current_local_prefix) && continue
            all((region_lower .<= candidate.values) .&
                (candidate.values .<= region_upper)) && push!(candidates, candidate)
        end

        starts = diverse_starts(candidates, local_starts, region_lower, region_upper)
        isempty(starts) && error("round $round_index produced no valid local start")
        println("round=$round_index selected_local_starts=$(length(starts))")
        flush(stdout)
        for (start_index, start) in enumerate(starts)
            search_start = (start.values .- region_center) ./ region_half_widths
            function local_objective(search_values)
                physical = region_center .+
                           Float64.(search_values) .* region_half_widths
                if all((region_lower .<= physical) .&
                       (physical .<= region_upper))
                    return evaluate_point(
                        physical;
                        source="round_$(round_index)_local_$(start_index)",
                    ).objective
                end
                below = max.((region_lower .- physical) ./ region_half_widths, 0.0)
                above = max.((physical .- region_upper) ./ region_half_widths, 0.0)
                return penalty + sum(abs2, below .+ above)
            end
            local_started = time()
            result = optimize(
                local_objective,
                search_start,
                NelderMead(
                    initial_simplex=Optim.AffineSimplexer(a=simplex_step, b=0.0),
                ),
                Optim.Options(
                    iterations=local_iterations,
                    f_reltol=1.0e-4,
                    g_tol=1.0e-6,
                    show_trace=false,
                    store_trace=false,
                ),
            )
            minimizer = clamp.(
                region_center .+ Optim.minimizer(result) .* region_half_widths,
                region_lower,
                region_upper,
            )
            endpoint = evaluate_point(
                minimizer;
                source="round_$(round_index)_local_$(start_index)_minimum",
            )
            local_endpoints[(round_index, start_index)] = endpoint
            push!(local_run_rows, Dict{String,Any}(
                "round" => round_index,
                "start_index" => start_index,
                "start" => start.values,
                "start_objective" => start.objective,
                "minimizer" => minimizer,
                "objective" => endpoint.objective,
                "valid" => endpoint.valid,
                "converged" => Optim.converged(result),
                "iterations" => Optim.iterations(result),
                "seconds" => time() - local_started,
            ))
        end

        valid_results = sort(
            filter(row -> row.valid && isfinite(row.objective), collect(values(cache)));
            by=row -> row.objective,
        )
        isempty(valid_results) && error("search produced no identity-safe point")
        round_best = first(valid_results)
        near_boundary = any(
            min(round_best.values[index] - region_lower[index],
                region_upper[index] - round_best.values[index]) <=
            boundary_fraction * (region_upper[index] - region_lower[index])
            for index in eachindex(parameters)
        )
        if round_index < maximum_rounds && near_boundary
            region_center = copy(round_best.values)
            region_half_widths = min.(
                region_half_widths .* expansion_factor,
                (upper .- lower) ./ 2,
            )
            println("expanding_trust_region=true new_center=$(join(region_center,','))")
            flush(stdout)
        else
            break
        end
    end

    valid_results = sort(
        filter(row -> row.valid && isfinite(row.objective), collect(values(cache)));
        by=row -> row.objective,
    )
    isempty(valid_results) && error("conformal optimization produced no valid point")
    best_recorded = first(valid_results)

    # The decisive result is recomputed after clearing every eigensolver warm
    # vector; only coupling-independent spaces and exact operator matrices remain.
    empty!(problem.warm_vectors)
    best = evaluate_point(
        best_recorded.values; source="best_cold_recheck", force=true,
    )
    cold_recheck_delta = abs(best.objective - best_recorded.objective)
    cold_recheck_passed = best.valid &&
        cold_recheck_delta <= cold_recheck_tolerance

    final_round_runs = filter(
        run -> Int(run["round"]) == rounds_completed, local_run_rows,
    )
    consensus_rows = NamedTuple[]
    for run in local_run_rows
        round_index = Int(run["round"])
        start_index = Int(run["start_index"])
        considered = round_index == rounds_completed
        endpoint = local_endpoints[(round_index, start_index)]
        parameter_deltas = abs.(endpoint.values .- best.values)
        agrees = considered && endpoint.valid &&
            all(parameter_deltas .<= consensus_parameter_tolerances) &&
            endpoint.objective <= best.objective + consensus_objective_tolerance
        push!(consensus_rows, (
            round=round_index,
            start_index=start_index,
            considered=considered,
            valid=endpoint.valid,
            agrees_with_best=agrees,
            objective=endpoint.objective,
            delta_objective=endpoint.objective - best.objective,
            Uf=endpoint.values[1],
            Uf0=endpoint.values[2],
            Vf0=endpoint.values[3],
            muc=endpoint.values[4],
            delta_Uf=parameter_deltas[1],
            delta_Uf0=parameter_deltas[2],
            delta_Vf0=parameter_deltas[3],
            delta_muc=parameter_deltas[4],
        ))
    end
    MottJainED.atomic_csv(
        joinpath(output, "multistart_convergence.csv"), DataFrame(consensus_rows),
    )
    final_consensus_rows = filter(row -> row.considered, consensus_rows)
    multistart_consensus_fraction = isempty(final_consensus_rows) ? 0.0 :
        count(row -> row.agrees_with_best, final_consensus_rows) /
        length(final_consensus_rows)
    multistart_all_agree = !isempty(final_consensus_rows) &&
        all(row -> row.agrees_with_best, final_consensus_rows)
    multistart_complete = length(final_consensus_rows) == local_starts

    neighbor_rows = NamedTuple[]
    local_minimum_confirmed = true
    for (index, parameter) in enumerate(parameters)
        for direction in (-1, 1)
            neighbor_values = copy(best.values)
            neighbor_values[index] += direction * robustness_steps[index]
            if !(lower[index] <= neighbor_values[index] <= upper[index])
                local_minimum_confirmed = false
                push!(neighbor_rows, (
                    parameter=String(parameter), direction=direction,
                    value=neighbor_values[index], valid=false, objective=Inf,
                    delta_objective=Inf, holdout_mean=Inf,
                ))
                continue
            end
            neighbor = evaluate_point(
                neighbor_values;
                source="robustness_$(parameter)_$(direction > 0 ? "plus" : "minus")",
                force=true,
            )
            acceptable = neighbor.valid &&
                neighbor.objective >= best.objective - robustness_objective_tolerance
            local_minimum_confirmed &= acceptable
            push!(neighbor_rows, (
                parameter=String(parameter), direction=direction,
                value=neighbor_values[index], valid=neighbor.valid,
                objective=neighbor.objective,
                delta_objective=neighbor.objective - best.objective,
                holdout_mean=isnothing(neighbor.holdout) ? Inf : neighbor.holdout.mean,
            ))
        end
    end
    MottJainED.atomic_csv(
        joinpath(output, "robustness_neighbors.csv"), DataFrame(neighbor_rows),
    )

    top_rows = NamedTuple[]
    for (rank, candidate) in enumerate(first(
        valid_results, min(top_candidate_count, length(valid_results)),
    ))
        push!(top_rows, (
            rank=rank,
            point_id=candidate.id,
            objective=candidate.objective,
            training_mean=candidate.training.mean,
            training_worst=candidate.training.worst,
            holdout_mean=candidate.holdout.mean,
            minimum_identity_overlap=candidate.minimum_identity_overlap,
            factor=candidate.factor,
            Uf=candidate.values[1],
            U0=u0_over_uf * candidate.values[1],
            Uf0=candidate.values[2],
            Vf0=candidate.values[3],
            V0=0.0,
            muc=candidate.values[4],
        ))
    end
    MottJainED.atomic_csv(joinpath(output, "top_candidates.csv"), DataFrame(top_rows))

    MottJainED.atomic_csv(
        joinpath(output, "best_primary_residuals.csv"),
        DataFrame(best.result.primary_rows),
    )
    MottJainED.atomic_csv(
        joinpath(output, "best_commutator_residuals.csv"),
        DataFrame(best.result.mixed_commutator_rows),
    )
    MottJainED.atomic_csv(
        joinpath(output, "best_commutator_components.csv"),
        DataFrame(best.result.mixed_commutator_pair_rows),
    )
    MottJainED.atomic_csv(
        joinpath(output, "best_training_constraints.csv"), DataFrame(best.training.rows),
    )
    MottJainED.atomic_csv(
        joinpath(output, "best_holdout_constraints.csv"), DataFrame(best.holdout.rows),
    )
    MottJainED.atomic_csv(
        joinpath(output, "best_identity_overlaps.csv"), DataFrame(best.identity.rows),
    )

    # Legacy energy relations are computed only after algebraic selection.  They
    # never enter sampling, local optimization, ranking, or acceptance.
    spectrum_five = score_cft_blocks(best.result.energies)
    spectrum_six = score_cft_blocks(
        best.result.energies; terms=CFT_STABLE_SIX_TERMS,
    )
    spectrum_rows = NamedTuple[]
    for (family, score) in (("five", spectrum_five), ("stable_six", spectrum_six))
        for index in eachindex(score.labels)
            push!(spectrum_rows, (
                family=family,
                term=String(score.terms[index]),
                label=score.labels[index],
                scaled_gap=score.scaled_gaps[index],
                target=score.target_gaps[index],
                residual=score.residuals[index],
            ))
        end
    end
    MottJainED.atomic_csv(
        joinpath(output, "best_spectrum_diagnostics.csv"), DataFrame(spectrum_rows),
    )

    hard_clearances = min.(
        (best.values .- lower) ./ (upper .- lower),
        (upper .- best.values) ./ (upper .- lower),
    )
    hard_boundary_passed = minimum(hard_clearances) > hard_boundary_fraction
    objective_improved = best.objective < anchor.objective
    holdout_guard_passed = best.holdout.mean <=
        holdout_maximum_anchor_ratio * anchor.holdout.mean
    accepted = best.valid && objective_improved && cold_recheck_passed &&
        holdout_guard_passed && hard_boundary_passed &&
        local_minimum_confirmed && multistart_complete && multistart_all_agree

    best_dict = Dict{String,Any}(
        "completed_at" => string(now()),
        "accepted" => accepted,
        "signature" => signature,
        "nm1" => nm,
        "seed" => seed,
        "heavy_space_mode" => String(heavy_space_mode),
        "parameters" => String.(parameters),
        "training_primaries" => String.(training_labels),
        "holdout_primaries" => String.(holdout_labels),
        "objective" => best.objective,
        "objective_improved" => objective_improved,
        "training_mean" => best.training.mean,
        "training_worst" => best.training.worst,
        "holdout_mean" => best.holdout.mean,
        "holdout_worst" => best.holdout.worst,
        "holdout_guard_passed" => holdout_guard_passed,
        "holdout_maximum_anchor_ratio" => holdout_maximum_anchor_ratio,
        "factor" => best.factor,
        "algebra_fit_loss" => best.algebra_fit_loss,
        "minimum_identity_overlap" => best.minimum_identity_overlap,
        "anchor_objective" => anchor.objective,
        "anchor_training_mean" => anchor.training.mean,
        "anchor_holdout_mean" => anchor.holdout.mean,
        "cold_recheck_delta" => cold_recheck_delta,
        "cold_recheck_tolerance" => cold_recheck_tolerance,
        "cold_recheck_passed" => cold_recheck_passed,
        "hard_boundary_clearances" => hard_clearances,
        "hard_boundary_passed" => hard_boundary_passed,
        "local_minimum_confirmed" => local_minimum_confirmed,
        "multistart_all_agree" => multistart_all_agree,
        "multistart_complete" => multistart_complete,
        "multistart_requested_starts" => local_starts,
        "multistart_final_start_count" => length(final_round_runs),
        "multistart_consensus_fraction" => multistart_consensus_fraction,
        "multistart_final_round" => rounds_completed,
        "consensus_parameter_tolerances" => consensus_parameter_tolerances,
        "consensus_objective_tolerance" => consensus_objective_tolerance,
        "spectrum_five_q_diagnostic_only" => spectrum_five.q,
        "spectrum_stable_six_q_diagnostic_only" => spectrum_six.q,
        "spectrum_delta_s_diagnostic_only" => spectrum_six.delta_s,
        "spectrum_delta_o_diagnostic_only" => spectrum_six.delta_o,
        "evaluation_count" => evaluation[],
        "cached_search_point_count" => length(cache),
        "initial_samples_per_round" => initial_samples,
        "maximum_rounds" => maximum_rounds,
        "rounds_completed" => rounds_completed,
        "local_iterations" => local_iterations,
        "local_runs" => local_run_rows,
        "couplings" => Dict(
            String(name) => getfield(best.couplings, name)
            for name in fieldnames(Couplings)
        ),
        "term_weights" => Dict(String(key) => value for (key, value) in term_weights),
        "config_source_sha256" => sha256_file(config_path),
        "so3lver_source_sha256" => sha256_file(so3_path),
        "driver_source_sha256" => sha256_file(@__FILE__),
        "project_git_revision" => MottJainED.git_revision(PROJECT_ROOT),
        "fuzzified_version" => string(Base.pkgversion(FuzzifiED)),
    )
    MottJainED.atomic_toml(joinpath(output, "best.toml"), best_dict)

    println("search_accepted=$accepted")
    println("objective_improved=$objective_improved")
    println("cold_recheck_passed=$cold_recheck_passed delta=$cold_recheck_delta")
    println("holdout_guard_passed=$holdout_guard_passed")
    println("hard_boundary_passed=$hard_boundary_passed")
    println("local_minimum_confirmed=$local_minimum_confirmed")
    println(
        "multistart_all_agree=$multistart_all_agree " *
        "multistart_complete=$multistart_complete " *
        "consensus_fraction=$multistart_consensus_fraction",
    )
    println("anchor_objective=$(anchor.objective) anchor_holdout=$(anchor.holdout.mean)")
    println("best_objective=$(best.objective) best_holdout=$(best.holdout.mean)")
    println(
        "diagnostic_only_spectrum_q5=$(spectrum_five.q) q6=$(spectrum_six.q)",
    )
    println(
        "best_parameters=Uf=$(best.couplings.Uf),U0=$(best.couplings.U0)," *
        "Uf0=$(best.couplings.Uf0),Vf0=$(best.couplings.Vf0)," *
        "V0=$(best.couplings.V0),muc=$(best.couplings.mu)",
    )
    println("optimization_result=$(joinpath(output, "best.toml"))")
    accepted || println(
        "SEARCH_NOT_ACCEPTED: inspect best.toml, multistart_convergence.csv, " *
        "and robustness_neighbors.csv",
    )
    return best_dict
end

run_search()
