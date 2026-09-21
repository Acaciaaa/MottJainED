#!/usr/bin/env julia

import Pkg

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(PROJECT_ROOT; io=devnull)

using CSV
using DataFrames
using Dates
using FuzzifiED
using MottJainED
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

function require_columns(frame, required, label)
    missing = setdiff(String.(required), names(frame))
    isempty(missing) || error("$label is missing columns: $(join(missing, ','))")
end

options = parse_options(ARGS)
config_path = abspath(get(
    options, "config",
    joinpath(PROJECT_ROOT, "config", "so3lver", "n6_nested_six_term_search.toml"),
))
isfile(config_path) || error("configuration not found: $config_path")
config = TOML.parsefile(config_path)
run_config = config["run"]
score_config = config["score"]
outer_config = config["outer"]
robustness_config = config["robustness"]
nm = Int(run_config["nm"])
k = Int(run_config["k"])
tol = Float64(run_config["tol"])
ncv = Int(run_config["ncv"])
base = dict_couplings(config["base"])
terms = Symbol.(score_config["terms"])
terms == collect(CFT_STABLE_SIX_TERMS) || error(
    "finalization requires the audited stable-six score",
)
maximum_q5_ratio = Float64(score_config["maximum_q5_ratio"])
outer_parameters = Symbol.(outer_config["parameters"])
top_candidate_count = Int(outer_config["top_candidate_count"])
robustness_tolerance = Float64(robustness_config["objective_tolerance"])

output = abspath(get(
    options, "output", joinpath(PROJECT_ROOT, String(run_config["output"])),
))
point_path = joinpath(output, "point_evaluations.csv")
residual_path = joinpath(output, "point_residuals.csv")
outer_path = joinpath(output, "outer_evaluations.csv")
robustness_path = joinpath(output, "robustness_neighbors.csv")
for path in (point_path, residual_path, outer_path, robustness_path)
    isfile(path) || error("required finalization checkpoint not found: $path")
end

points = CSV.read(point_path, DataFrame)
residuals = CSV.read(residual_path, DataFrame)
outer = CSV.read(outer_path, DataFrame)
robustness = CSV.read(robustness_path, DataFrame)
require_columns(points, (
    :signature, :evaluation, :outer_id, :source, :valid, :objective,
    :q6, :q5, :worst_residual, :factor, :delta_s, :delta_o,
    :minimum_expected_overlap, :Uf0, :Vf0, :V0, :mu,
), "point trace")
require_columns(residuals, (:evaluation, :term, :label, :residual), "residual trace")
require_columns(outer, (
    :signature, :outer_id, :source, :valid, :objective, :q6, :q5,
    :worst_residual, :minimum_expected_overlap, :best_mu,
    :mu_at_boundary, :Uf0, :Vf0, :V0,
), "outer trace")
require_columns(robustness, (
    :parameter, :direction, :valid, :objective, :delta_objective, :best_mu,
), "robustness trace")

signatures = unique(String.(outer.signature))
length(signatures) == 1 || error("outer trace contains multiple signatures")
all(String(value) == only(signatures) for value in points.signature) || error(
    "point and outer traces have different signatures",
)
signature = only(signatures)

best_rows = [row for row in eachrow(outer) if String(row.source) == "best_recheck"]
isempty(best_rows) && error("no completed best_recheck row is available")
best_outer = last(best_rows)
best_outer_id = String(best_outer.outer_id)
best_points = [
    row for row in eachrow(points)
    if String(row.outer_id) == best_outer_id &&
       startswith(String(row.source), "best_recheck") && Bool(row.valid)
]
isempty(best_points) && error("best_recheck has no valid point evaluations")
best_point = best_points[argmin(Float64[row.objective for row in best_points])]
isapprox(Float64(best_point.objective), Float64(best_outer.objective); rtol=1e-8, atol=1e-10) ||
    error("best point and outer objectives do not match")
isapprox(Float64(best_point.mu), Float64(best_outer.best_mu); rtol=1e-8, atol=1e-10) ||
    error("best point and outer chemical potentials do not match")

evaluation = Int(best_point.evaluation)
labels = String[]
best_residuals = Float64[]
for term in terms
    rows = [
        row for row in eachrow(residuals)
        if Int(row.evaluation) == evaluation && String(row.term) == String(term)
    ]
    length(rows) == 1 || error(
        "evaluation $evaluation has $(length(rows)) residual rows for term $term",
    )
    push!(labels, String(only(rows).label))
    push!(best_residuals, Float64(only(rows).residual))
end

expected_neighbors = Set(
    (String(parameter), direction)
    for parameter in outer_parameters for direction in (-1, 1)
)
actual_neighbors = Set(
    (String(row.parameter), Int(row.direction)) for row in eachrow(robustness)
)
nrow(robustness) == length(expected_neighbors) && actual_neighbors == expected_neighbors ||
    error("robustness checkpoint is incomplete or contains duplicate neighbors")
local_minimum_confirmed = all(eachrow(robustness)) do row
    Bool(row.valid) &&
        Float64(row.objective) >= Float64(best_outer.objective) - robustness_tolerance
end

# Recompute only the fixed anchor needed by the q5 acceptance guard.  No outer
# profile, local search, best recheck, or robustness eigensolve is repeated.
println("recomputing_anchor_only=true N=$nm")
problem = build_cft_problem(nm, base; disp_std=true)
anchor = solve_cft_blocks!(
    problem.hamiltonians, base; k, tol, ncv,
    warm_vectors=Dict{Tuple{Symbol,Int},Vector{Float64}}(), disp_std=true,
)
anchor_six = score_cft_blocks(anchor.energies; terms)
anchor_five = score_cft_blocks(anchor.energies; terms=CFT_SCORE_TERMS)
q5_guard_passed = Float64(best_outer.q5) <= maximum_q5_ratio * anchor_five.q
accepted = Bool(best_outer.valid) && !Bool(best_outer.mu_at_boundary) &&
           q5_guard_passed && local_minimum_confirmed

search_profiles = Dict{String,NamedTuple}()
for row in eachrow(outer)
    source = String(row.source)
    (source == "best_recheck" || startswith(source, "robustness_")) && continue
    search_profiles[String(row.outer_id)] = (
        id=String(row.outer_id), objective=Float64(row.objective),
        q6=Float64(row.q6), q5=Float64(row.q5),
        worst_residual=Float64(row.worst_residual),
        minimum_overlap=Float64(row.minimum_expected_overlap),
        mu_at_boundary=Bool(row.mu_at_boundary),
        values=Float64[row.Uf0, row.Vf0, row.V0],
        best_mu=Float64(row.best_mu), valid=Bool(row.valid),
    )
end
valid_profiles = sort(
    filter(row -> row.valid && isfinite(row.objective), collect(values(search_profiles)));
    by=row -> row.objective,
)
top_rows = NamedTuple[]
for (rank, candidate) in enumerate(first(
    valid_profiles, min(top_candidate_count, length(valid_profiles)),
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

best_dict = Dict{String,Any}(
    "completed_at" => string(now()),
    "completion_mode" => "recovered_from_complete_final_audit_checkpoints",
    "source_job_id" => get(options, "job-id", "unknown"),
    "accepted" => accepted,
    "signature" => signature,
    "nm1" => nm,
    "objective" => Float64(best_outer.objective),
    "q6" => Float64(best_outer.q6),
    "q5" => Float64(best_outer.q5),
    "anchor_q6" => anchor_six.q,
    "anchor_q5" => anchor_five.q,
    "q5_guard_passed" => q5_guard_passed,
    "mu_at_boundary" => Bool(best_outer.mu_at_boundary),
    "local_minimum_confirmed" => local_minimum_confirmed,
    "factor" => Float64(best_point.factor),
    "delta_s" => Float64(best_point.delta_s),
    "delta_o" => Float64(best_point.delta_o),
    "worst_residual" => Float64(best_outer.worst_residual),
    "minimum_expected_overlap" => Float64(best_outer.minimum_expected_overlap),
    "training_terms" => String.(terms),
    "labels" => labels,
    "residuals" => best_residuals,
    "parameters" => Dict(
        "Uf" => base.Uf, "Uf0" => Float64(best_outer.Uf0), "U0" => base.U0,
        "Vf" => base.Vf, "Vf0" => Float64(best_outer.Vf0),
        "V0" => Float64(best_outer.V0), "t" => base.t,
        "mu" => Float64(best_outer.best_mu),
    ),
    "outer_profile_count" => length(search_profiles),
    "point_evaluation_count" => nrow(points),
    "local_run_trace" => "outer_evaluations.csv",
    "local_runs_recovered" => false,
    "config_source_sha256" => sha256_file(config_path),
    "so3lver_source_sha256" => sha256_file(joinpath(PROJECT_ROOT, "experimental", "SO3lverED.jl")),
    "search_source_sha256" => sha256_file(joinpath(PROJECT_ROOT, "experimental", "SO3ParameterSearch.jl")),
    "finalizer_source_sha256" => sha256_file(@__FILE__),
    "project_git_revision" => MottJainED.git_revision(PROJECT_ROOT),
    "fuzzified_version" => string(Base.pkgversion(FuzzifiED)),
)
MottJainED.atomic_toml(joinpath(output, "best.toml"), best_dict)

println("recovered_finalization=true")
println("search_accepted=$accepted")
println("q5_guard_passed=$q5_guard_passed")
println("local_minimum_confirmed=$local_minimum_confirmed")
println("best_objective=$(best_outer.objective) q6=$(best_outer.q6) q5=$(best_outer.q5)")
println("best_parameters=Uf0=$(best_outer.Uf0),Vf0=$(best_outer.Vf0)," *
        "V0=$(best_outer.V0),mu=$(best_outer.best_mu)")
println("search_result=$(joinpath(output, "best.toml"))")
