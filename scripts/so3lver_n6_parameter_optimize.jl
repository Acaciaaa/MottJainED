#!/usr/bin/env julia

import Pkg

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(PROJECT_ROOT; io=devnull)

using CSV
using DataFrames
using Dates
using FuzzifiED
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

function sha256_file(path)
    return open(path) do io
        bytes2hex(sha256(io))
    end
end

function dict_couplings(values)
    defaults = Couplings()
    return Couplings(; (
        name => Float64(get(values, String(name), getfield(defaults, name)))
        for name in fieldnames(Couplings)
    )...)
end

function couplings_with_parameters(base, parameters, values)
    result = base
    for (parameter, value) in zip(parameters, values)
        result = MottJainED.with_coupling(result, parameter, value)
    end
    return result
end

function verify_audit!(summary, config_path, so3_path, search_path, audit_driver_path)
    Bool(summary["passed"]) || error(
        "reference audit did not certify ddS/boxS; seven-term optimization is forbidden",
    )
    checks = (
        ("config_source_sha256", config_path),
        ("so3lver_source_sha256", so3_path),
        ("search_source_sha256", search_path),
        ("driver_source_sha256", audit_driver_path),
    )
    for (field, path) in checks
        expected = String(summary[field])
        actual = sha256_file(path)
        actual == expected || error(
            "audit identity mismatch for $path: expected $expected, got $actual; rerun the audit",
        )
    end
    return nothing
end

options = parse_options(ARGS)
config_path = abspath(get(
    options, "config",
    joinpath(PROJECT_ROOT, "config", "so3lver", "n6_parameter_reference_audit.toml"),
))
isfile(config_path) || throw(ArgumentError("configuration not found: $config_path"))
config = TOML.parsefile(config_path)
run_config = config["run"]
audit_config = config["audit"]
optimization_config = config["optimization"]
base = dict_couplings(config["base"])
parameters = Symbol.(audit_config["free_parameters"])
terms = Symbol.(audit_config["score_terms"])
parameters == [:Uf0, :Vf0, :V0, :mu] || error(
    "this first audited search requires parameter order Uf0,Vf0,V0,mu",
)
bounds = Dict(
    Symbol(name) => (Float64(value[1]), Float64(value[2]))
    for (name, value) in audit_config["bounds"]
)
lower = Float64[bounds[parameter][1] for parameter in parameters]
upper = Float64[bounds[parameter][2] for parameter in parameters]
all(lower .< upper) || error("all optimization bounds must have positive width")
center = Float64[getfield(base, parameter) for parameter in parameters]

so3_path = joinpath(PROJECT_ROOT, "experimental", "SO3lverED.jl")
search_path = joinpath(PROJECT_ROOT, "experimental", "SO3ParameterSearch.jl")
audit_driver_path = joinpath(PROJECT_ROOT, "scripts", "so3lver_n6_parameter_audit.jl")
audit_path = abspath(joinpath(PROJECT_ROOT, String(optimization_config["audit_summary"])))
isfile(audit_path) || error("audit summary not found: $audit_path")
audit_summary = TOML.parsefile(audit_path)
verify_audit!(
    audit_summary, config_path, so3_path, search_path, audit_driver_path,
)

nm = Int(run_config["nm"])
k = Int(run_config["k"])
tol = Float64(run_config["tol"])
ncv = Int(run_config["ncv"])
max_iterations = Int(optimization_config["max_iterations_per_start"])
simplex_step = Float64(optimization_config["simplex_step"])
scale_config = optimization_config["parameter_scales"]
parameter_scales = Float64[
    scale_config[String(parameter)] for parameter in parameters
]
all(parameter_scales .> 0) || error(
    "all optimization parameter scales must be positive",
)
worst_weight = Float64(optimization_config["worst_residual_weight"])
penalty = Float64(optimization_config["penalty"])
minimum_overlap = Float64(optimization_config["minimum_overlap"])
wide_mu_count = Int(optimization_config["wide_mu_count"])
wide_tolerance = Float64(optimization_config["wide_disagreement_tolerance"])
starts = [Float64.(values) for values in optimization_config["starts"]]
all(length(values) == length(parameters) for values in starts) || error(
    "every optimization start must contain Uf0,Vf0,V0,mu",
)
all(values -> all((lower .<= values) .& (values .<= upper)), starts) || error(
    "every optimization start must lie inside the configured bounds",
)
output = abspath(get(
    options, "output", joinpath(PROJECT_ROOT, String(optimization_config["output"])),
))
mkpath(output)
Random.seed!(Int(run_config["seed"]))

signature = MottJainED.stable_id(
    "so3lver-n6-seven-term-search-v1", nm, k, tol, ncv,
    String.(parameters), String.(terms), lower, upper, starts,
    max_iterations, simplex_step, parameter_scales, worst_weight, minimum_overlap,
    sha256_file(config_path), sha256_file(so3_path), sha256_file(search_path),
)
trace_path = joinpath(output, "evaluations.csv")
residual_path = joinpath(output, "evaluation_residuals.csv")
evaluation = Ref(0)
if isfile(trace_path)
    previous = CSV.read(trace_path, DataFrame)
    "signature" in names(previous) || error(
        "existing trace has no signature; choose a new output directory",
    )
    all(String(value) == signature for value in previous.signature) || error(
        "existing trace belongs to another search; choose a new output directory",
    )
    evaluation[] = nrow(previous)
    valid_rows = filter(row -> Bool(row.valid) && isfinite(row.objective), previous)
    if nrow(valid_rows) > 0
        best = valid_rows[argmin(valid_rows.objective), :]
        resumed = Float64[best[parameter] for parameter in parameters]
        pushfirst!(starts, resumed)
        @info "resuming from best recorded point" objective=best.objective parameters=resumed
    end
end

println("N=$nm seven-term SO(3)lver multi-start optimization")
println("certified_audit=$audit_path")
println("output=$output")
println("starts=$(length(starts)) max_iterations_per_start=$max_iterations")
println("Julia $(VERSION), threads=$(Threads.nthreads()), FuzzifiED $(Base.pkgversion(FuzzifiED))")
flush(stdout)

problem = build_cft_problem(nm, base; disp_std=true)
anchor = solve_cft_blocks!(
    problem.hamiltonians, base; k, tol, ncv,
    warm_vectors=Dict{Tuple{Symbol,Int},Vector{Float64}}(), disp_std=true,
)
anchor_score = score_cft_blocks(anchor.energies; terms=terms)
warm_vectors = Dict(
    key => copy(anchor.vectors[key][:, 1]) for key in CFT_BLOCK_KEYS
)

function record_evaluation!(source, values, couplings, solved, training, five, tracking,
                            objective, valid, reason, seconds)
    minimum_expected_overlap = isnothing(tracking) ? NaN :
        minimum(row.expected_overlap for row in tracking.rows)
    worst_residual = isnothing(training) ? Inf : maximum(abs, training.residuals)
    parameter_fields = (;
        (parameter => Float64(value) for (parameter, value) in zip(parameters, values))...,
    )
    MottJainED.append_csv(trace_path, (
        signature=signature,
        evaluation=evaluation[],
        source=String(source),
        timestamp=string(now()),
        valid=valid,
        reason=reason,
        objective=Float64(objective),
        q_training=isnothing(training) ? Inf : training.q,
        q5=isnothing(five) ? Inf : five.q,
        worst_residual=worst_residual,
        factor=isnothing(training) ? NaN : training.factor,
        delta_s=isnothing(training) ? NaN : training.delta_s,
        delta_o=isnothing(training) ? NaN : training.delta_o,
        ground_is_singlet_l0=isnothing(training) ? false : training.ground_is_singlet_l0,
        identity_passed=isnothing(tracking) ? false : tracking.passed,
        minimum_expected_overlap=minimum_expected_overlap,
        seconds=seconds,
        parameter_fields...,
    ))
    if !isnothing(training)
        for (index, term) in enumerate(training.terms)
            MottJainED.append_csv(residual_path, (
                signature=signature,
                evaluation=evaluation[],
                source=String(source),
                term=String(term),
                label=training.labels[index],
                raw_gap=training.raw_gaps[index],
                target=training.target_gaps[index],
                scaled_gap=training.scaled_gaps[index],
                residual=training.residuals[index],
            ))
        end
    end
    return nothing
end

function evaluate_physical(values; source="optimizer")
    evaluation[] += 1
    if !all((lower .<= values) .& (values .<= upper))
        record_evaluation!(
            source, values, base, nothing, nothing, nothing, nothing,
            penalty, false, "outside_bounds", 0.0,
        )
        return (objective=penalty, valid=false, score=nothing, five=nothing,
                tracking=nothing, couplings=base)
    end
    couplings = couplings_with_parameters(base, parameters, values)
    started = time()
    try
        solved = solve_cft_blocks!(
            problem.hamiltonians, couplings; k, tol, ncv,
            warm_vectors, disp_std=true,
        )
        training = score_cft_blocks(solved.energies; terms=terms)
        five = score_cft_blocks(solved.energies; terms=CFT_SCORE_TERMS)
        tracking = track_reference_states(
            anchor.vectors, solved.vectors;
            minimum_overlap, require_same_rank=true,
        )
        valid = tracking.passed && training.ground_is_singlet_l0
        reason = valid ? "ok" : (
            !tracking.passed ? "state_identity_gate" : "ground_branch_gate"
        )
        raw_objective = training.q + worst_weight * maximum(abs, training.residuals)
        objective = valid ? raw_objective : penalty + raw_objective
        seconds = time() - started
        record_evaluation!(
            source, values, couplings, solved, training, five, tracking,
            objective, valid, reason, seconds,
        )
        println("evaluation=$(evaluation[]) source=$source objective=$objective " *
                "q=$(training.q) valid=$valid values=$(join(values, ','))")
        flush(stdout)
        return (objective=objective, valid=valid, score=training, five=five,
                tracking=tracking, couplings=couplings)
    catch err
        seconds = time() - started
        reason = replace(sprint(showerror, err), '\n' => ' ')
        record_evaluation!(
            source, values, couplings, nothing, nothing, nothing, nothing,
            penalty, false, reason, seconds,
        )
        @error "parameter evaluation failed" values exception=(err, catch_backtrace())
        return (objective=penalty, valid=false, score=nothing, five=nothing,
                tracking=nothing, couplings=couplings)
    end
end

run_summaries = Dict{String,Any}[]
for (start_index, physical_start) in enumerate(starts)
    search_start = physical_to_parameter_search(
        physical_start, center, parameter_scales,
    )
    objective(search_values) = begin
        physical = parameter_search_to_physical(
            search_values, center, parameter_scales,
        )
        if all((lower .<= physical) .& (physical .<= upper))
            evaluate_physical(physical; source="start_$start_index").objective
        else
            below = max.((lower .- physical) ./ parameter_scales, 0)
            above = max.((physical .- upper) ./ parameter_scales, 0)
            penalty + sum(abs2, below .+ above)
        end
    end
    method = NelderMead(
        initial_simplex=Optim.AffineSimplexer(a=simplex_step, b=0.0),
    )
    started = time()
    result = optimize(
        objective, search_start, method,
        Optim.Options(
            iterations=max_iterations,
            f_reltol=1.0e-5,
            g_tol=1.0e-7,
            show_trace=false,
            store_trace=false,
        ),
    )
    minimizer = parameter_search_to_physical(
        Optim.minimizer(result), center, parameter_scales,
    )
    push!(run_summaries, Dict{String,Any}(
        "start_index" => start_index,
        "start" => physical_start,
        "minimizer" => minimizer,
        "minimum" => Optim.minimum(result),
        "converged" => Optim.converged(result),
        "iterations" => Optim.iterations(result),
        "seconds" => time() - started,
    ))
    println("start=$start_index converged=$(Optim.converged(result)) " *
            "minimum=$(Optim.minimum(result)) minimizer=$(join(minimizer, ','))")
    flush(stdout)
end

evaluations = CSV.read(trace_path, DataFrame)
valid_rows = filter(row -> Bool(row.valid) && isfinite(row.objective), evaluations)
nrow(valid_rows) > 0 || error("optimization produced no identity-safe valid point")
best_row = valid_rows[argmin(valid_rows.objective), :]
best_values = Float64[best_row[parameter] for parameter in parameters]
best = evaluate_physical(best_values; source="best_recheck")
best.valid || error("best recorded point did not pass deterministic recheck")

mu_index = findfirst(==(:mu), parameters)
mu_grid = collect(range(lower[mu_index], upper[mu_index]; length=wide_mu_count))
guard_rows = NamedTuple[]
for mu in mu_grid
    values = copy(best_values)
    values[mu_index] = mu
    result = evaluate_physical(values; source="wide_mu_guard")
    push!(guard_rows, (
        mu=mu,
        valid=result.valid,
        objective=result.objective,
        q=isnothing(result.score) ? Inf : result.score.q,
        minimum_expected_overlap=isnothing(result.tracking) ? NaN :
            minimum(row.expected_overlap for row in result.tracking.rows),
    ))
end
MottJainED.atomic_csv(joinpath(output, "wide_mu_guard.csv"), DataFrame(guard_rows))
valid_guard = filter(row -> row.valid && isfinite(row.objective), guard_rows)
guard_best = isempty(valid_guard) ? nothing : valid_guard[argmin(getproperty.(valid_guard, :objective))]
wide_disagreement = !isnothing(guard_best) &&
                    guard_best.objective + wide_tolerance < best.objective
accepted = best.valid && !wide_disagreement

best_dict = Dict{String,Any}(
    "completed_at" => string(now()),
    "accepted" => accepted,
    "wide_disagreement" => wide_disagreement,
    "signature" => signature,
    "nm1" => nm,
    "training_terms" => String.(terms),
    "objective" => best.objective,
    "q_training" => best.score.q,
    "q5" => best.five.q,
    "factor" => best.score.factor,
    "delta_s" => best.score.delta_s,
    "delta_o" => best.score.delta_o,
    "worst_residual" => maximum(abs, best.score.residuals),
    "residuals" => best.score.residuals,
    "labels" => best.score.labels,
    "minimum_expected_overlap" => minimum(
        row.expected_overlap for row in best.tracking.rows
    ),
    "evaluations" => evaluation[],
    "run_summaries" => run_summaries,
    "wide_guard_best_mu" => isnothing(guard_best) ? NaN : guard_best.mu,
    "wide_guard_best_objective" => isnothing(guard_best) ? Inf : guard_best.objective,
    "parameters" => Dict(
        String(parameter) => value for (parameter, value) in zip(parameters, best_values)
    ),
    "fixed_parameters" => Dict(
        "Uf" => base.Uf, "U0" => base.U0, "Vf" => base.Vf, "t" => base.t,
    ),
    "audit_summary" => audit_path,
    "config_source_sha256" => sha256_file(config_path),
    "so3lver_source_sha256" => sha256_file(so3_path),
    "search_source_sha256" => sha256_file(search_path),
    "driver_source_sha256" => sha256_file(@__FILE__),
    "fuzzified_expected_revision" => audit_summary["fuzzified_expected_revision"],
)
MottJainED.atomic_toml(joinpath(output, "best.toml"), best_dict)

println("optimization_accepted=$accepted")
println("wide_disagreement=$wide_disagreement")
println("best_objective=$(best.objective) q=$(best.score.q)")
println("best_parameters=$(join(best_values, ','))")
println("optimization_result=$(joinpath(output, "best.toml"))")
accepted || println("SEARCH_NOT_ACCEPTED: inspect the wide-mu guard and evaluation trace")
