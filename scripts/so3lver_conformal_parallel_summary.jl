#!/usr/bin/env julia

import Pkg

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(PROJECT_ROOT; io=devnull)

using DataFrames
using Dates
using MottJainED
using TOML

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

options = parse_options(ARGS)
haskey(options, "root") || error("--root is required")
root = abspath(options["root"])
config_path = abspath(get(
    options, "config",
    joinpath(
        PROJECT_ROOT, "config", "so3lver",
        "n6_projected_conformal_optimization.toml",
    ),
))
worker_count = parse(Int, get(options, "workers", "6"))
worker_count >= 2 || error("parallel summary needs at least two workers")
isdir(root) || error("parallel result root not found: $root")
isfile(config_path) || error("configuration not found: $config_path")

config = TOML.parsefile(config_path)
run = config["run"]
search = config["search"]
base = config["base"]
parameter_names = ["Uf", "Uf0", "Vf0", "mu"]
u0_over_uf = Float64(search["u0_over_uf"])
parameter_tolerances = Float64[
    search["consensus_parameter_tolerances"][name]
    for name in parameter_names
]
objective_tolerance = Float64(search["consensus_objective_tolerance"])

worker_results = Dict{String,Any}[]
for worker in 1:worker_count
    best_path = joinpath(root, "worker-$worker", "best.toml")
    isfile(best_path) || error(
        "worker $worker did not produce best.toml: $best_path",
    )
    best = TOML.parsefile(best_path)
    Int(get(best, "worker_index", 0)) == worker || error(
        "worker $worker output has the wrong worker_index",
    )
    Int(get(best, "parallel_workers", 0)) == worker_count || error(
        "worker $worker output was not produced for $worker_count workers",
    )
    String(get(best, "heavy_space_mode", "")) == "laughlin13" || error(
        "worker $worker did not use the Laughlin-1/3 projected space",
    )
    Int(get(best, "nm1", 0)) == Int(run["nm"]) || error(
        "worker $worker used the wrong system size",
    )
    couplings = best["couplings"]
    values = Float64[couplings[name] for name in parameter_names]
    objective = Float64(best["objective"])
    isfinite(objective) || error("worker $worker has a non-finite objective")
    push!(worker_results, Dict{String,Any}(
        "worker" => worker,
        "accepted" => Bool(best["accepted"]),
        "objective" => objective,
        "values" => values,
        "holdout_mean" => Float64(best["holdout_mean"]),
        "minimum_identity_overlap" => Float64(best["minimum_identity_overlap"]),
        "cold_recheck_passed" => Bool(best["cold_recheck_passed"]),
        "holdout_guard_passed" => Bool(best["holdout_guard_passed"]),
        "hard_boundary_passed" => Bool(best["hard_boundary_passed"]),
        "local_minimum_confirmed" => Bool(best["local_minimum_confirmed"]),
        "signature" => String(best["signature"]),
        "project_git_revision" => String(best["project_git_revision"]),
        "config_source_sha256" => String(best["config_source_sha256"]),
        "so3lver_source_sha256" => String(best["so3lver_source_sha256"]),
        "driver_source_sha256" => String(best["driver_source_sha256"]),
        "fuzzified_version" => String(best["fuzzified_version"]),
    ))
end

for field in (
    "project_git_revision", "config_source_sha256", "so3lver_source_sha256",
    "driver_source_sha256", "fuzzified_version",
)
    values = unique(String(result[field]) for result in worker_results)
    length(values) == 1 || error(
        "parallel workers disagree on $field: $(join(values, ", "))",
    )
end

best_result = worker_results[argmin(
    Float64[result["objective"] for result in worker_results],
)]
best_values = Float64.(best_result["values"])
best_objective = Float64(best_result["objective"])

rows = NamedTuple[]
summary_workers = Dict{String,Any}[]
for result in worker_results
    values = Float64.(result["values"])
    deltas = abs.(values .- best_values)
    delta_objective = Float64(result["objective"]) - best_objective
    agrees = all(deltas .<= parameter_tolerances) &&
             delta_objective <= objective_tolerance
    push!(rows, (
        worker=Int(result["worker"]),
        accepted=Bool(result["accepted"]),
        agrees_with_best=agrees,
        objective=Float64(result["objective"]),
        delta_objective=delta_objective,
        holdout_mean=Float64(result["holdout_mean"]),
        Uf=values[1],
        Uf0=values[2],
        Vf0=values[3],
        muc=values[4],
        delta_Uf=deltas[1],
        delta_Uf0=deltas[2],
        delta_Vf0=deltas[3],
        delta_muc=deltas[4],
        minimum_identity_overlap=Float64(result["minimum_identity_overlap"]),
        cold_recheck_passed=Bool(result["cold_recheck_passed"]),
        holdout_guard_passed=Bool(result["holdout_guard_passed"]),
        hard_boundary_passed=Bool(result["hard_boundary_passed"]),
        local_minimum_confirmed=Bool(result["local_minimum_confirmed"]),
    ))
    push!(summary_workers, Dict{String,Any}(
        "worker" => Int(result["worker"]),
        "accepted" => Bool(result["accepted"]),
        "agrees_with_best" => agrees,
        "objective" => Float64(result["objective"]),
        "delta_objective" => delta_objective,
        "Uf" => values[1],
        "Uf0" => values[2],
        "Vf0" => values[3],
        "muc" => values[4],
        "delta_Uf" => deltas[1],
        "delta_Uf0" => deltas[2],
        "delta_Vf0" => deltas[3],
        "delta_muc" => deltas[4],
    ))
end

all_workers_agree = all(row -> row.agrees_with_best, rows)
all_workers_accepted = all(row -> row.accepted, rows)
accepted = all_workers_agree && all_workers_accepted
MottJainED.atomic_csv(
    joinpath(root, "parallel_convergence.csv"), DataFrame(rows),
)

summary = Dict{String,Any}(
    "completed_at" => string(now()),
    "accepted" => accepted,
    "worker_count" => worker_count,
    "all_workers_agree" => all_workers_agree,
    "all_workers_accepted" => all_workers_accepted,
    "best_worker" => Int(best_result["worker"]),
    "best_objective" => best_objective,
    "best_holdout_mean" => Float64(best_result["holdout_mean"]),
    "best_parameters" => Dict(
        "Uf" => best_values[1],
        "U0" => u0_over_uf * best_values[1],
        "Uf0" => best_values[2],
        "Vf" => Float64(base["Vf"]),
        "Vf0" => best_values[3],
        "V0" => Float64(base["V0"]),
        "t" => Float64(base["t"]),
        "muc" => best_values[4],
    ),
    "consensus_parameter_tolerances" => parameter_tolerances,
    "consensus_objective_tolerance" => objective_tolerance,
    "workers" => summary_workers,
    "project_git_revisions" => sort(unique(
        String(result["project_git_revision"]) for result in worker_results
    )),
)
MottJainED.atomic_toml(joinpath(root, "parallel_summary.toml"), summary)

println("parallel_search_accepted=$accepted")
println("all_workers_agree=$all_workers_agree " *
        "all_workers_accepted=$all_workers_accepted")
println("best_worker=$(best_result["worker"]) " *
        "best_objective=$best_objective")
println(
    "best_parameters=Uf=$(best_values[1])," *
    "U0=$(u0_over_uf * best_values[1])," *
    "Uf0=$(best_values[2]),Vf0=$(best_values[3])," *
    "V0=$(base["V0"])," *
    "muc=$(best_values[4])",
)
println("parallel_summary=$(joinpath(root, "parallel_summary.toml"))")
accepted || println(
    "PARALLEL_SEARCH_NOT_ACCEPTED: inspect parallel_convergence.csv and worker outputs",
)
