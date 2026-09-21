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

# The first production N=6 run wrote 407 reusable point evaluations before a
# missing LinearAlgebra import stopped the diverse-start selection.  Keep the
# original driver hash in the resume signature because importing `norm` only
# restores the already intended algorithm.  Any future scientific change to
# this driver must deliberately update this compatibility token.
const NESTED_SEARCH_SIGNATURE_DRIVER_SHA256 =
    "59f0f080ed6c085b696a4fde033b760eb37a5b820879e40b9c2f063fa34824bb"

function dict_couplings(values)
    defaults = Couplings()
    return Couplings(; (
        name => Float64(get(values, String(name), getfield(defaults, name)))
        for name in fieldnames(Couplings)
    )...)
end

function with_values(base, parameters, values)
    couplings = base
    for (parameter, value) in zip(parameters, values)
        couplings = MottJainED.with_coupling(couplings, parameter, value)
    end
    return couplings
end

function interpolate_couplings(left, right, fraction)
    return Couplings(; (
        name => (1 - fraction) * getfield(left, name) +
                fraction * getfield(right, name)
        for name in fieldnames(Couplings)
    )...)
end

function finite_maximum(values)
    isempty(values) && return Inf
    return maximum(abs, values)
end

options = parse_options(ARGS)
config_path = abspath(get(
    options, "config",
    joinpath(PROJECT_ROOT, "config", "so3lver", "n6_nested_six_term_search.toml"),
))
isfile(config_path) || throw(ArgumentError("configuration not found: $config_path"))
config = TOML.parsefile(config_path)
run_config = config["run"]
score_config = config["score"]
mu_config = config["mu"]
outer_config = config["outer"]
continuation_config = config["continuation"]
robustness_config = config["robustness"]

nm = Int(run_config["nm"])
allow_test_size = lowercase(get(options, "allow-test-size", "false")) == "true"
(nm == 6 || allow_test_size) || error(
    "this driver is intentionally restricted to N=6; " *
    "--allow-test-size=true is only for a small end-to-end smoke test",
)
k = Int(run_config["k"])
k >= 3 || error("the stable six-term score needs at least three levels per block")
tol = Float64(run_config["tol"])
ncv = Int(run_config["ncv"])
seed = Int(run_config["seed"])
base = dict_couplings(config["base"])

terms = Symbol.(score_config["terms"])
terms == collect(CFT_STABLE_SIX_TERMS) || error(
    "the nested search must use exactly the audited stable six relations",
)
worst_weight = Float64(score_config["worst_residual_weight"])
minimum_overlap = Float64(score_config["minimum_overlap"])
penalty = Float64(score_config["penalty"])
maximum_q5_ratio = Float64(score_config["maximum_q5_ratio"])

mu_lower = Float64(mu_config["lower"])
mu_upper = Float64(mu_config["upper"])
mu_lower < mu_upper || error("mu bounds must have positive width")
mu_grid_count = Int(mu_config["grid_count"])
mu_grid_count >= 3 || error("mu grid needs at least three points")
mu_refine_basins = Int(mu_config["refine_basins"])
mu_refine_iterations = Int(mu_config["refine_iterations"])
mu_absolute_tolerance = Float64(mu_config["absolute_tolerance"])
mu_boundary_tolerance = Float64(mu_config["boundary_tolerance"])

outer_parameters = Symbol.(outer_config["parameters"])
outer_parameters == [:Uf0, :Vf0, :V0] || error(
    "outer parameters must be ordered as Uf0,Vf0,V0",
)
initial_samples = Int(outer_config["initial_samples"])
local_starts = Int(outer_config["local_starts"])
local_iterations = Int(outer_config["local_iterations"])
simplex_step = Float64(outer_config["simplex_step"])
maximum_rounds = Int(outer_config["maximum_rounds"])
expansion_factor = Float64(outer_config["expansion_factor"])
boundary_fraction = Float64(outer_config["boundary_fraction"])
minimum_start_distance = Float64(outer_config["minimum_start_distance"])
top_candidate_count = Int(outer_config["top_candidate_count"])

initial_width_config = outer_config["initial_half_widths"]
hard_bounds_config = outer_config["hard_bounds"]
initial_half_widths = Float64[
    initial_width_config[String(parameter)] for parameter in outer_parameters
]
hard_lower = Float64[
    hard_bounds_config[String(parameter)][1] for parameter in outer_parameters
]
hard_upper = Float64[
    hard_bounds_config[String(parameter)][2] for parameter in outer_parameters
]
all(initial_half_widths .> 0) || error("initial half-widths must be positive")
all(hard_lower .< hard_upper) || error("hard bounds must have positive width")
anchor_values = Float64[getfield(base, parameter) for parameter in outer_parameters]
all((hard_lower .<= anchor_values) .& (anchor_values .<= hard_upper)) || error(
    "anchor outer parameters must lie inside hard bounds",
)

continuation_steps_config = continuation_config["steps"]
continuation_steps = Dict(
    Symbol(name) => Float64(value)
    for (name, value) in continuation_steps_config
)
all(value > 0 for value in values(continuation_steps)) || error(
    "continuation steps must be positive",
)
maximum_continuation_steps = Int(continuation_config["maximum_steps"])

robustness_steps_config = robustness_config["steps"]
robustness_steps = Float64[
    robustness_steps_config[String(parameter)] for parameter in outer_parameters
]
robustness_tolerance = Float64(robustness_config["objective_tolerance"])

output = abspath(get(
    options, "output", joinpath(PROJECT_ROOT, String(run_config["output"])),
))
validate_only = lowercase(get(options, "validate-only", "false")) == "true"

so3_path = joinpath(PROJECT_ROOT, "experimental", "SO3lverED.jl")
search_path = joinpath(PROJECT_ROOT, "experimental", "SO3ParameterSearch.jl")
signature = MottJainED.stable_id(
    "so3lver-n6-nested-six-term-v1", nm, k, tol, ncv, seed,
    String.(terms), String.(outer_parameters), anchor_values,
    hard_lower, hard_upper, initial_half_widths,
    mu_lower, mu_upper, mu_grid_count, mu_refine_basins,
    sha256_file(config_path), sha256_file(so3_path), sha256_file(search_path),
    NESTED_SEARCH_SIGNATURE_DRIVER_SHA256,
)

point_path = joinpath(output, "point_evaluations.csv")
residual_path = joinpath(output, "point_residuals.csv")
tracking_path = joinpath(output, "state_tracking.csv")
outer_path = joinpath(output, "outer_evaluations.csv")
evaluation = Ref(0)
if isfile(point_path)
    previous = CSV.read(point_path, DataFrame)
    "signature" in names(previous) || error(
        "existing point trace has no signature; choose a new output directory",
    )
    all(String(value) == signature for value in previous.signature) || error(
        "existing point trace belongs to another search; choose a new output directory",
    )
    evaluation[] = nrow(previous)
end
if isfile(outer_path)
    previous = CSV.read(outer_path, DataFrame)
    "signature" in names(previous) || error(
        "existing outer trace has no signature; choose a new output directory",
    )
    all(String(value) == signature for value in previous.signature) || error(
        "existing outer trace belongs to another search; choose a new output directory",
    )
end

println("N=$nm nested six-term SO(3)lver parameter search")
println("outer_parameters=$(join(outer_parameters, ',')); mu is fully reprofiled")
println("output=$output")
println("Julia $(VERSION), threads=$(Threads.nthreads()), FuzzifiED $(Base.pkgversion(FuzzifiED))")
flush(stdout)
if validate_only
    println("VALIDATION_OK: configuration, dependencies, and source signatures loaded")
    exit(0)
end
mkpath(output)
cp(config_path, joinpath(output, "search_config.toml"); force=true)

problem = build_cft_problem(nm, base; disp_std=true)
anchor = solve_cft_blocks!(
    problem.hamiltonians, base; k, tol, ncv,
    warm_vectors=Dict{Tuple{Symbol,Int},Vector{Float64}}(), disp_std=true,
)
anchor_six = score_cft_blocks(anchor.energies; terms)
anchor_five = score_cft_blocks(anchor.energies; terms=CFT_SCORE_TERMS)
anchor_six.ground_is_singlet_l0 || error("anchor is not on the singlet L=0 ground branch")

anchor_warm() = Dict(
    key => copy(anchor.vectors[key][:, 1]) for key in CFT_BLOCK_KEYS
)

function characterize(solved; reference_vectors=anchor.vectors)
    six = score_cft_blocks(solved.energies; terms)
    five = score_cft_blocks(solved.energies; terms=CFT_SCORE_TERMS)
    tracking = track_reference_states(
        reference_vectors, solved.vectors;
        specifications=STABLE_SIX_TRACKED_STATE_SPECS,
        minimum_overlap, require_same_rank=true,
    )
    min_overlap = minimum(row.expected_overlap for row in tracking.rows)
    valid = tracking.passed && six.ground_is_singlet_l0
    reason = valid ? "ok" : (
        !tracking.passed ? "state_identity_gate" : "ground_branch_gate"
    )
    raw_objective = six.q + worst_weight * finite_maximum(six.residuals)
    return (
        solved=solved, six=six, five=five, tracking=tracking,
        min_overlap=min_overlap, valid=valid, reason=reason,
        raw_objective=raw_objective,
    )
end

function required_continuation_steps(target)
    counts = Int[]
    for (parameter, step) in continuation_steps
        delta = abs(getfield(target, parameter) - getfield(base, parameter))
        push!(counts, ceil(Int, delta / step))
    end
    return max(1, maximum(counts))
end

function continuation_assessment(target)
    count = required_continuation_steps(target)
    count <= maximum_continuation_steps || return (
        assessment=nothing, minimum_overlap=0.0,
        reason="continuation_requires_$(count)_steps",
    )
    previous_vectors = anchor.vectors
    warm = anchor_warm()
    path_minimum = 1.0
    final_assessment = nothing
    for index in 1:count
        couplings = interpolate_couplings(base, target, index / count)
        solved = solve_cft_blocks!(
            problem.hamiltonians, couplings; k, tol, ncv,
            warm_vectors=warm, disp_std=false,
        )
        assessment = characterize(solved; reference_vectors=previous_vectors)
        path_minimum = min(path_minimum, assessment.min_overlap)
        if !assessment.valid
            return (
                assessment=assessment, minimum_overlap=path_minimum,
                reason="continuation_step_$(index)_$(assessment.reason)",
            )
        end
        previous_vectors = solved.vectors
        final_assessment = assessment
    end
    return (
        assessment=final_assessment, minimum_overlap=path_minimum,
        reason="ok",
    )
end

function record_point!(source, outer_id, couplings, assessment, objective,
                       valid, reason, identity_mode, min_overlap, seconds)
    evaluation[] += 1
    six = isnothing(assessment) ? nothing : assessment.six
    five = isnothing(assessment) ? nothing : assessment.five
    tracking = isnothing(assessment) ? nothing : assessment.tracking
    MottJainED.append_csv(point_path, (
        signature=signature,
        evaluation=evaluation[],
        outer_id=String(outer_id),
        source=String(source),
        timestamp=string(now()),
        valid=valid,
        reason=String(reason),
        identity_mode=String(identity_mode),
        objective=Float64(objective),
        q6=isnothing(six) ? Inf : six.q,
        q5=isnothing(five) ? Inf : five.q,
        worst_residual=isnothing(six) ? Inf : finite_maximum(six.residuals),
        factor=isnothing(six) ? NaN : six.factor,
        delta_s=isnothing(six) ? NaN : six.delta_s,
        delta_o=isnothing(six) ? NaN : six.delta_o,
        ground_is_singlet_l0=isnothing(six) ? false : six.ground_is_singlet_l0,
        minimum_expected_overlap=Float64(min_overlap),
        seconds=Float64(seconds),
        Uf0=couplings.Uf0,
        Vf0=couplings.Vf0,
        V0=couplings.V0,
        mu=couplings.mu,
    ))
    if !isnothing(six)
        for index in eachindex(six.terms)
            MottJainED.append_csv(residual_path, (
                signature=signature,
                evaluation=evaluation[],
                outer_id=String(outer_id),
                source=String(source),
                term=String(six.terms[index]),
                label=six.labels[index],
                raw_gap=six.raw_gaps[index],
                target=six.target_gaps[index],
                scaled_gap=six.scaled_gaps[index],
                residual=six.residuals[index],
            ))
        end
    end
    if !isnothing(tracking)
        for row in tracking.rows
            MottJainED.append_csv(tracking_path, (
                signature=signature,
                evaluation=evaluation[],
                outer_id=String(outer_id),
                source=String(source),
                identity_mode=String(identity_mode),
                label=String(row.label),
                representation=String(row.representation),
                ell=row.ell,
                expected_rank=row.expected_rank,
                best_rank=row.best_rank,
                expected_overlap=row.expected_overlap,
                best_overlap=row.best_overlap,
                same_rank=row.same_rank,
                passed=row.passed,
            ))
        end
    end
    return nothing
end

function evaluate_target(couplings, warm; source, outer_id)
    started = time()
    assessment = nothing
    identity_mode = "direct"
    minimum_path_overlap = 0.0
    reason = "unknown"
    try
        solved = solve_cft_blocks!(
            problem.hamiltonians, couplings; k, tol, ncv,
            warm_vectors=warm, disp_std=false,
        )
        assessment = characterize(solved)
        minimum_path_overlap = assessment.min_overlap
        reason = assessment.reason
        if !assessment.valid && assessment.reason == "state_identity_gate"
            continued = continuation_assessment(couplings)
            identity_mode = "continuation"
            minimum_path_overlap = continued.minimum_overlap
            reason = continued.reason
            if continued.reason == "ok" && !isnothing(continued.assessment)
                assessment = continued.assessment
                reason = assessment.reason
            end
        end
        valid = assessment.valid && reason == "ok"
        raw_objective = assessment.raw_objective
        objective = valid ? raw_objective : penalty + raw_objective
        record_point!(
            source, outer_id, couplings, assessment, objective, valid, reason,
            identity_mode, minimum_path_overlap, time() - started,
        )
        println("point=$(evaluation[]) source=$source outer=$outer_id " *
                "objective=$objective valid=$valid mu=$(couplings.mu)")
        flush(stdout)
        # The search only consumes the score and tracking summaries after the
        # point has been recorded.  Do not retain the much larger eigensystem
        # inside every mu/profile checkpoint during the local search.
        compact_assessment = merge(assessment, (solved=nothing,))
        return (
            objective=objective, valid=valid, reason=reason,
            assessment=compact_assessment, couplings=couplings,
            identity_mode=identity_mode,
            minimum_overlap=minimum_path_overlap,
        )
    catch err
        reason = replace(sprint(showerror, err), '\n' => ' ')
        record_point!(
            source, outer_id, couplings, assessment, penalty, false, reason,
            identity_mode, minimum_path_overlap, time() - started,
        )
        @error "N=$nm point evaluation failed" outer_id source couplings exception=(err, catch_backtrace())
        return (
            objective=penalty, valid=false, reason=reason,
            assessment=nothing, couplings=couplings,
            identity_mode=identity_mode, minimum_overlap=minimum_path_overlap,
        )
    end
end

function outer_id(values)
    rounded = round.(Float64.(values); digits=11)
    return MottJainED.stable_id(signature, rounded)
end

profiles = Dict{String,Any}()
if isfile(outer_path)
    previous = CSV.read(outer_path, DataFrame)
    "signature" in names(previous) || error(
        "existing outer trace has no signature; choose a new output directory",
    )
    all(String(value) == signature for value in previous.signature) || error(
        "existing outer trace belongs to another search; choose a new output directory",
    )
    excluded_final_audits = Ref(0)
    for row in eachrow(previous)
        source = String(row.source)
        if source == "best_recheck" || startswith(source, "robustness_")
            # These rows are downstream audits of the already selected search
            # minimum.  If a job dies while writing the final files, feeding a
            # slightly better robustness neighbor back into the optimizer on
            # resume would silently change the search trajectory.
            excluded_final_audits[] += 1
            continue
        end
        values = Float64[row.Uf0, row.Vf0, row.V0]
        profiles[String(row.outer_id)] = (
            id=String(row.outer_id), values=values, valid=Bool(row.valid),
            reason=String(row.reason), best_mu=Float64(row.best_mu),
            objective=Float64(row.objective), q6=Float64(row.q6),
            q5=Float64(row.q5), worst_residual=Float64(row.worst_residual),
            minimum_overlap=Float64(row.minimum_expected_overlap),
            mu_at_boundary=Bool(row.mu_at_boundary), best_result=nothing,
        )
    end
    println("resumed_complete_outer_profiles=$(length(profiles)) " *
            "excluded_final_audit_rows=$(excluded_final_audits[])")
end

function record_outer!(summary, source, seconds)
    MottJainED.append_csv(outer_path, (
        signature=signature,
        outer_id=summary.id,
        source=String(source),
        timestamp=string(now()),
        valid=summary.valid,
        reason=summary.reason,
        objective=summary.objective,
        q6=summary.q6,
        q5=summary.q5,
        worst_residual=summary.worst_residual,
        minimum_expected_overlap=summary.minimum_overlap,
        best_mu=summary.best_mu,
        mu_at_boundary=summary.mu_at_boundary,
        seconds=seconds,
        Uf0=summary.values[1],
        Vf0=summary.values[2],
        V0=summary.values[3],
    ))
    return nothing
end

function profile_mu(values; source="outer", force=false)
    physical = Float64.(values)
    id = outer_id(physical)
    !force && haskey(profiles, id) && return profiles[id]
    started = time()
    fixed = with_values(base, outer_parameters, physical)
    warm = anchor_warm()
    mu_results = Dict{String,Any}()
    function evaluate_mu(mu, mu_source)
        bounded_mu = Float64(mu)
        key = string(round(bounded_mu; digits=12))
        haskey(mu_results, key) && return mu_results[key]
        couplings = MottJainED.with_coupling(fixed, :mu, bounded_mu)
        result = evaluate_target(
            couplings, warm; source="$(source)_$(mu_source)", outer_id=id,
        )
        mu_results[key] = result
        return result
    end

    grid = collect(range(mu_lower, mu_upper; length=mu_grid_count))
    grid_results = [evaluate_mu(mu, "grid") for mu in grid]
    brackets = mu_refinement_brackets(
        grid, getproperty.(grid_results, :objective),
        getproperty.(grid_results, :valid); maximum_count=mu_refine_basins,
    )
    for (basin_index, bracket) in enumerate(brackets)
        objective(mu) = evaluate_mu(mu, "refine_$(basin_index)").objective
        try
            result = optimize(
                objective, bracket.lower, bracket.upper, Brent();
                abs_tol=mu_absolute_tolerance, iterations=mu_refine_iterations,
            )
            evaluate_mu(Optim.minimizer(result), "refine_$(basin_index)_minimum")
        catch err
            @warn "mu basin refinement failed" id basin_index exception=(err, catch_backtrace())
        end
    end

    valid_results = filter(result -> result.valid && isfinite(result.objective),
                           collect(Base.values(mu_results)))
    if isempty(valid_results)
        summary = (
            id=id, values=physical, valid=false, reason="no_valid_mu",
            best_mu=NaN, objective=penalty, q6=Inf, q5=Inf,
            worst_residual=Inf, minimum_overlap=0.0,
            mu_at_boundary=true, best_result=nothing,
        )
    else
        best = valid_results[argmin(getproperty.(valid_results, :objective))]
        six = best.assessment.six
        summary = (
            id=id, values=physical, valid=true, reason="ok",
            best_mu=best.couplings.mu, objective=best.objective,
            q6=six.q, q5=best.assessment.five.q,
            worst_residual=finite_maximum(six.residuals),
            minimum_overlap=best.minimum_overlap,
            mu_at_boundary=min(best.couplings.mu - mu_lower,
                               mu_upper - best.couplings.mu) <= mu_boundary_tolerance,
            best_result=best,
        )
    end
    record_outer!(summary, source, time() - started)
    profiles[id] = summary
    println("outer=$id source=$source objective=$(summary.objective) " *
            "mu=$(summary.best_mu) valid=$(summary.valid)")
    flush(stdout)
    return summary
end

function diverse_starts(candidates, count, lower, upper)
    valid = sort(
        filter(row -> row.valid && isfinite(row.objective), candidates);
        by=row -> row.objective,
    )
    selected = Any[]
    width = upper .- lower
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

function run_nested_search()
    rng = MersenneTwister(seed)
    region_center = copy(anchor_values)
    region_half_widths = copy(initial_half_widths)
    local_run_rows = Dict{String,Any}[]

for round_index in 1:maximum_rounds
    region_lower = max.(hard_lower, region_center .- region_half_widths)
    region_upper = min.(hard_upper, region_center .+ region_half_widths)
    println("round=$round_index lower=$(join(region_lower, ',')) " *
            "upper=$(join(region_upper, ','))")
    flush(stdout)

    candidates = Any[]
    push!(candidates, profile_mu(region_center; source="round_$(round_index)_center"))
    samples = latin_hypercube_points(
        initial_samples, region_lower, region_upper, rng,
    )
    for row_index in axes(samples, 1)
        push!(candidates, profile_mu(
            vec(samples[row_index, :]);
            source="round_$(round_index)_lhs_$(row_index)",
        ))
    end
    for profile in values(profiles)
        all((region_lower .<= profile.values) .&
            (profile.values .<= region_upper)) && push!(candidates, profile)
    end

    starts = diverse_starts(candidates, local_starts, region_lower, region_upper)
    isempty(starts) && error("round $round_index produced no valid nested-mu start")
    for (start_index, start) in enumerate(starts)
        search_start = (start.values .- region_center) ./ region_half_widths
        function local_objective(search_values)
            physical = region_center .+ Float64.(search_values) .* region_half_widths
            if all((region_lower .<= physical) .& (physical .<= region_upper))
                return profile_mu(
                    physical;
                    source="round_$(round_index)_local_$(start_index)",
                ).objective
            end
            below = max.((region_lower .- physical) ./ region_half_widths, 0)
            above = max.((physical .- region_upper) ./ region_half_widths, 0)
            return penalty + sum(abs2, below .+ above)
        end
        method = NelderMead(
            initial_simplex=Optim.AffineSimplexer(a=simplex_step, b=0.0),
        )
        started = time()
        result = optimize(
            local_objective, search_start, method,
            Optim.Options(
                iterations=local_iterations, f_reltol=1.0e-4,
                g_tol=1.0e-6, show_trace=false, store_trace=false,
            ),
        )
        minimizer = region_center .+ Optim.minimizer(result) .* region_half_widths
        push!(local_run_rows, Dict{String,Any}(
            "round" => round_index,
            "start_index" => start_index,
            "start" => start.values,
            "minimizer" => minimizer,
            "minimum" => Optim.minimum(result),
            "converged" => Optim.converged(result),
            "iterations" => Optim.iterations(result),
            "seconds" => time() - started,
        ))
    end

    valid_profiles = sort(
        filter(row -> row.valid && isfinite(row.objective), collect(values(profiles)));
        by=row -> row.objective,
    )
    isempty(valid_profiles) && error("search produced no identity-safe profile")
    round_best = first(valid_profiles)
    near_boundary = any(
        min(round_best.values[index] - region_lower[index],
            region_upper[index] - round_best.values[index]) <=
        boundary_fraction * (region_upper[index] - region_lower[index])
        for index in eachindex(outer_parameters)
    )
    if round_index < maximum_rounds && near_boundary
        region_center = copy(round_best.values)
        region_half_widths = min.(
            region_half_widths .* expansion_factor,
            (hard_upper .- hard_lower) ./ 2,
        )
        println("expanding_trust_region=true new_center=$(join(region_center, ','))")
        flush(stdout)
    else
        break
    end
end

valid_profiles = sort(
    filter(row -> row.valid && isfinite(row.objective), collect(values(profiles)));
    by=row -> row.objective,
)
isempty(valid_profiles) && error("nested search produced no valid candidate")
best_recorded = first(valid_profiles)

# Re-run the complete mu profile, rather than merely reusing the optimizer's
# last warm state, before applying the acceptance guards.
best = profile_mu(best_recorded.values; source="best_recheck", force=true)
best.valid || error("best candidate failed its complete mu-profile recheck")

neighbor_rows = NamedTuple[]
local_minimum_confirmed = true
for (index, parameter) in enumerate(outer_parameters)
    for direction in (-1, 1)
        neighbor_values = copy(best.values)
        neighbor_values[index] += direction * robustness_steps[index]
        if !(hard_lower[index] <= neighbor_values[index] <= hard_upper[index])
            local_minimum_confirmed = false
            push!(neighbor_rows, (
                parameter=String(parameter), direction=direction,
                value=neighbor_values[index], valid=false, objective=Inf,
                delta_objective=Inf, best_mu=NaN,
            ))
            continue
        end
        neighbor = profile_mu(
            neighbor_values;
            source="robustness_$(parameter)_$(direction > 0 ? "plus" : "minus")",
            force=true,
        )
        acceptable = neighbor.valid &&
                     neighbor.objective >= best.objective - robustness_tolerance
        local_minimum_confirmed &= acceptable
        push!(neighbor_rows, (
            parameter=String(parameter), direction=direction,
            value=neighbor_values[index], valid=neighbor.valid,
            objective=neighbor.objective,
            delta_objective=neighbor.objective - best.objective,
            best_mu=neighbor.best_mu,
        ))
    end
end
MottJainED.atomic_csv(joinpath(output, "robustness_neighbors.csv"), DataFrame(neighbor_rows))

all_valid_profiles = sort(
    filter(row -> row.valid && isfinite(row.objective), collect(values(profiles)));
    by=row -> row.objective,
)
top_rows = NamedTuple[]
for (rank, candidate) in enumerate(first(
    all_valid_profiles, min(top_candidate_count, length(all_valid_profiles)),
))
    push!(top_rows, (
        rank=rank, outer_id=candidate.id, objective=candidate.objective,
        q6=candidate.q6, q5=candidate.q5,
        worst_residual=candidate.worst_residual,
        minimum_expected_overlap=candidate.minimum_overlap,
        mu_at_boundary=candidate.mu_at_boundary,
        Uf0=candidate.values[1], Vf0=candidate.values[2],
        V0=candidate.values[3], mu=candidate.best_mu,
    ))
end
MottJainED.atomic_csv(joinpath(output, "top_candidates.csv"), DataFrame(top_rows))

q5_guard_passed = best.q5 <= maximum_q5_ratio * anchor_five.q
accepted = best.valid && !best.mu_at_boundary && q5_guard_passed &&
           local_minimum_confirmed
best_result = best.best_result
six = best_result.assessment.six

best_dict = Dict{String,Any}(
    "completed_at" => string(now()),
    "accepted" => accepted,
    "signature" => signature,
    "nm1" => nm,
    "objective" => best.objective,
    "q6" => best.q6,
    "q5" => best.q5,
    "anchor_q6" => anchor_six.q,
    "anchor_q5" => anchor_five.q,
    "q5_guard_passed" => q5_guard_passed,
    "mu_at_boundary" => best.mu_at_boundary,
    "local_minimum_confirmed" => local_minimum_confirmed,
    "factor" => six.factor,
    "delta_s" => six.delta_s,
    "delta_o" => six.delta_o,
    "worst_residual" => best.worst_residual,
    "minimum_expected_overlap" => best.minimum_overlap,
    "training_terms" => String.(terms),
    "labels" => six.labels,
    "residuals" => six.residuals,
    "parameters" => Dict(
        "Uf" => base.Uf, "Uf0" => best.values[1], "U0" => base.U0,
        "Vf" => base.Vf, "Vf0" => best.values[2], "V0" => best.values[3],
        "t" => base.t, "mu" => best.best_mu,
    ),
    "outer_profile_count" => length(profiles),
    "point_evaluation_count" => evaluation[],
    "local_runs" => local_run_rows,
    "config_source_sha256" => sha256_file(config_path),
    "so3lver_source_sha256" => sha256_file(so3_path),
    "search_source_sha256" => sha256_file(search_path),
    "driver_source_sha256" => sha256_file(@__FILE__),
    "driver_signature_sha256" => NESTED_SEARCH_SIGNATURE_DRIVER_SHA256,
    "project_git_revision" => MottJainED.git_revision(PROJECT_ROOT),
    "fuzzified_version" => string(Base.pkgversion(FuzzifiED)),
)
MottJainED.atomic_toml(joinpath(output, "best.toml"), best_dict)

println("search_accepted=$accepted")
println("q5_guard_passed=$q5_guard_passed")
println("local_minimum_confirmed=$local_minimum_confirmed")
println("best_objective=$(best.objective) q6=$(best.q6) q5=$(best.q5)")
println("best_parameters=Uf0=$(best.values[1]),Vf0=$(best.values[2])," *
        "V0=$(best.values[3]),mu=$(best.best_mu)")
println("search_result=$(joinpath(output, "best.toml"))")
accepted || println("SEARCH_NOT_ACCEPTED: inspect best.toml and robustness_neighbors.csv")
    return best_dict
end

run_nested_search()
