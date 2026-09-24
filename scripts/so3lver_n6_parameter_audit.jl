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

function source_revision()
    project = TOML.parsefile(joinpath(PROJECT_ROOT, "Project.toml"))
    source = project["sources"]["FuzzifiED"]
    return get(source, "rev", "unknown")
end

function dict_couplings(values)
    defaults = Couplings()
    return Couplings(; (
        name => Float64(get(values, String(name), getfield(defaults, name)))
        for name in fieldnames(Couplings)
    )...)
end

function set_parameter(couplings::Couplings, parameter::Symbol, value::Real)
    parameter in fieldnames(Couplings) || throw(ArgumentError(
        "unknown Hamiltonian parameter $parameter",
    ))
    return MottJainED.with_coupling(couplings, parameter, value)
end

function score_rows(point_name, parameter, direction, value, score, family)
    return [(
        point=point_name,
        parameter=String(parameter),
        direction=direction,
        value=Float64(value),
        score_family=family,
        term=String(term),
        label=score.labels[index],
        raw_gap=score.raw_gaps[index],
        target=score.target_gaps[index],
        scaled_gap=score.scaled_gaps[index],
        residual=score.residuals[index],
        q=score.q,
        factor=score.factor,
        delta_s=score.delta_s,
        delta_o=score.delta_o,
        ground_is_singlet_l0=score.ground_is_singlet_l0,
    ) for (index, term) in enumerate(score.terms)]
end

function tracking_rows(point_name, parameter, direction, value, tracking)
    return [(
        point=point_name,
        parameter=String(parameter),
        direction=direction,
        value=Float64(value),
        label=String(row.label),
        representation=String(row.representation),
        ell=row.ell,
        expected_rank=row.expected_rank,
        best_rank=row.best_rank,
        expected_overlap=row.expected_overlap,
        best_overlap=row.best_overlap,
        same_rank=row.same_rank,
        passed=row.passed,
    ) for row in tracking.rows]
end

function generator_dict(result)
    return Dict{String,Any}(
        "passed" => result.passed,
        "fit_fidelity" => result.fit_fidelity,
        "numerical_rank" => result.numerical_rank,
        "candidate_names" => String.(result.candidate_names),
        "coefficients" => result.coefficients,
        "singular_values" => result.singular_values,
        "l0_overlaps_by_rank" => result.l0_overlaps,
        "l2_overlaps_by_rank" => result.l2_overlaps,
        "l0_expected_subspace_overlap" => result.l0_expected_subspace_overlap,
        "l2_expected_subspace_overlap" => result.l2_expected_subspace_overlap,
        "l0_expected_subspace_fraction" => result.l0_expected_subspace_fraction,
        "l2_expected_subspace_fraction" => result.l2_expected_subspace_fraction,
        "l0_resolved_overlap" => result.l0_resolved_overlap,
        "l2_resolved_overlap" => result.l2_resolved_overlap,
        "l0_unresolved_overlap_upper_bound" => result.l0_unresolved_overlap_upper_bound,
        "l2_unresolved_overlap_upper_bound" => result.l2_unresolved_overlap_upper_bound,
        "boxs_competitor_overlap_upper_bound" => result.boxs_competitor_overlap_upper_bound,
        "dds_competitor_overlap_upper_bound" => result.dds_competitor_overlap_upper_bound,
        "boxs_leading_margin" => result.boxs_leading_margin,
        "dds_leading_margin" => result.dds_leading_margin,
        "boxs_is_leading_non_s" => result.boxs_is_leading_non_s,
        "dds_is_leading_non_t" => result.dds_is_leading_non_t,
        "low_level_count" => result.low_level_count,
        "minimum_l0_level_count" => result.minimum_l0_level_count,
        "minimum_l2_level_count" => result.minimum_l2_level_count,
        "l0_level_count" => result.l0_level_count,
        "l2_level_count" => result.l2_level_count,
        "l0_available_level_count" => result.l0_available_level_count,
        "l2_available_level_count" => result.l2_available_level_count,
    )
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
generator_config = config["generator"]
nm = Int(run_config["nm"])
k = Int(run_config["k"])
tol = Float64(run_config["tol"])
ncv = Int(run_config["ncv"])
base = dict_couplings(config["base"])
parameters = Symbol.(audit_config["free_parameters"])
terms = Symbol.(audit_config["score_terms"])
steps = Dict(Symbol(name) => Float64(value) for (name, value) in audit_config["steps"])
bounds = Dict(
    Symbol(name) => (Float64(value[1]), Float64(value[2]))
    for (name, value) in audit_config["bounds"]
)
minimum_overlap = Float64(audit_config["minimum_overlap"])
require_same_rank = Bool(audit_config["require_same_rank"])
require_full_rank = Bool(audit_config["require_full_jacobian_rank"])
maximum_condition_number = Float64(audit_config["maximum_jacobian_condition_number"])
maximum_condition_number > 1 || error(
    "audit.maximum_jacobian_condition_number must be greater than one",
)
all(haskey(steps, parameter) for parameter in parameters) || error(
    "every free parameter needs an audit step",
)
all(haskey(bounds, parameter) for parameter in parameters) || error(
    "every free parameter needs an audit bound",
)
for parameter in parameters
    lower, upper = bounds[parameter]
    center = getfield(base, parameter)
    lower <= center - steps[parameter] < center + steps[parameter] <= upper || error(
        "audit probes for $parameter lie outside the configured bounds",
    )
end

output = abspath(get(options, "output", joinpath(PROJECT_ROOT, String(run_config["output"]))))
mkpath(output)
Random.seed!(Int(run_config["seed"]))

println("N=$nm SO(3)lver parameter-reference audit")
println("config=$config_path")
println("output=$output")
println("Julia $(VERSION), threads=$(Threads.nthreads()), FuzzifiED $(Base.pkgversion(FuzzifiED))")
flush(stdout)

generator_enabled = Bool(generator_config["enabled"])
generator_error = Ref("")
generator_result = if generator_enabled
    println("generator_anchor_audit=started")
    flush(stdout)
    started = time()
    try
        result = audit_scalar_generator(
            nm, base;
            k=Int(generator_config["k"]),
            eig_tol=tol,
            low_level_count=Int(generator_config["low_level_count"]),
            minimum_l0_level_count=Int(generator_config["minimum_l0_level_count"]),
            minimum_l2_level_count=Int(generator_config["minimum_l2_level_count"]),
            minimum_fit_fidelity=Float64(generator_config["minimum_fit_fidelity"]),
            minimum_expected_subspace_overlap=Float64(
                generator_config["minimum_expected_subspace_overlap"],
            ),
            minimum_expected_subspace_fraction=Float64(
                generator_config["minimum_expected_subspace_fraction"],
            ),
        )
        println("generator_anchor_audit=finished passed=$(result.passed) seconds=$(time()-started)")
        flush(stdout)
        result
    catch err
        generator_error[] = replace(sprint(showerror, err), '\n' => ' ')
        @error "generator anchor audit failed" exception=(err, catch_backtrace())
        println("generator_anchor_audit=failed seconds=$(time()-started)")
        flush(stdout)
        nothing
    end
else
    nothing
end
GC.gc(true)

println("so3_workspace=started")
flush(stdout)
problem = build_cft_problem(nm, base; disp_std=true)
println("so3_workspace=finished")
flush(stdout)

anchor_started = time()
anchor = solve_cft_blocks!(
    problem.hamiltonians, base; k, tol, ncv,
    warm_vectors=Dict{Tuple{Symbol,Int},Vector{Float64}}(), disp_std=true,
)
anchor_seconds = time() - anchor_started
anchor_five = score_cft_blocks(anchor.energies; terms=CFT_SCORE_TERMS)
anchor_training = score_cft_blocks(anchor.energies; terms=terms)
println("anchor q5=$(anchor_five.q) q_training=$(anchor_training.q) seconds=$anchor_seconds")
flush(stdout)

score_table = NamedTuple[]
tracking_table = NamedTuple[]
append!(score_table, score_rows("anchor", :none, 0, NaN, anchor_five, "fixed_five"))
append!(score_table, score_rows("anchor", :none, 0, NaN, anchor_training, "training"))

plus_scores = Dict{Symbol,Any}()
minus_scores = Dict{Symbol,Any}()
probe_summaries = Dict{String,Any}[]
all_tracking_passed = Ref(true)
all_ground_passed = Ref(anchor_training.ground_is_singlet_l0)

for parameter in parameters, direction in (-1, 1)
    value = getfield(base, parameter) + direction * steps[parameter]
    point_couplings = set_parameter(base, parameter, value)
    point_name = "$(parameter)_$(direction < 0 ? "minus" : "plus")"
    warm = Dict(
        key => copy(anchor.vectors[key][:, 1]) for key in CFT_BLOCK_KEYS
    )
    started = time()
    solved = solve_cft_blocks!(
        problem.hamiltonians, point_couplings; k, tol, ncv,
        warm_vectors=warm, disp_std=true,
    )
    elapsed = time() - started
    five = score_cft_blocks(solved.energies; terms=CFT_SCORE_TERMS)
    training = score_cft_blocks(solved.energies; terms=terms)
    tracking = track_reference_states(
        anchor.vectors, solved.vectors;
        minimum_overlap, require_same_rank,
    )
    direction < 0 ? (minus_scores[parameter] = training) :
                    (plus_scores[parameter] = training)
    append!(score_table, score_rows(
        point_name, parameter, direction, value, five, "fixed_five",
    ))
    append!(score_table, score_rows(
        point_name, parameter, direction, value, training, "training",
    ))
    append!(tracking_table, tracking_rows(
        point_name, parameter, direction, value, tracking,
    ))
    all_tracking_passed[] &= tracking.passed
    all_ground_passed[] &= training.ground_is_singlet_l0
    minimum_point_overlap = minimum(row.expected_overlap for row in tracking.rows)
    push!(probe_summaries, Dict{String,Any}(
        "point" => point_name,
        "parameter" => String(parameter),
        "direction" => direction,
        "value" => value,
        "seconds" => elapsed,
        "q5" => five.q,
        "q_training" => training.q,
        "factor_training" => training.factor,
        "ground_is_singlet_l0" => training.ground_is_singlet_l0,
        "tracking_passed" => tracking.passed,
        "minimum_expected_overlap" => minimum_point_overlap,
    ))
    println("probe=$point_name value=$value q=$(training.q) " *
            "tracking=$(tracking.passed) min_overlap=$minimum_point_overlap seconds=$elapsed")
    flush(stdout)
end

jacobian = normalized_residual_jacobian(plus_scores, minus_scores, parameters)
jacobian_rank_passed = !require_full_rank ||
                       jacobian.numerical_rank == length(parameters)
jacobian_condition_passed = isfinite(jacobian.condition_number) &&
                            jacobian.condition_number <= maximum_condition_number
jacobian_passed = jacobian_rank_passed && jacobian_condition_passed
generator_passed = !generator_enabled ||
                   (!isnothing(generator_result) && generator_result.passed)
passed = generator_passed && all_tracking_passed[] && all_ground_passed[] && jacobian_passed

MottJainED.atomic_csv(joinpath(output, "scores.csv"), DataFrame(score_table))
MottJainED.atomic_csv(joinpath(output, "state_tracking.csv"), DataFrame(tracking_table))
jacobian_rows = NamedTuple[]
for (row_index, term) in enumerate(anchor_training.terms),
    (column_index, parameter) in enumerate(parameters)
    push!(jacobian_rows, (
        term=String(term), parameter=String(parameter),
        normalized_derivative=jacobian.matrix[row_index, column_index],
    ))
end
MottJainED.atomic_csv(joinpath(output, "normalized_jacobian.csv"), DataFrame(jacobian_rows))

summary = Dict{String,Any}(
    "completed_at" => string(now()),
    "passed" => passed,
    "generator_passed" => generator_passed,
    "generator_enabled" => generator_enabled,
    "generator_error" => generator_error[],
    "state_tracking_passed" => all_tracking_passed[],
    "ground_branch_passed" => all_ground_passed[],
    "jacobian_passed" => jacobian_passed,
    "jacobian_rank_passed" => jacobian_rank_passed,
    "jacobian_condition_passed" => jacobian_condition_passed,
    "nm1" => nm,
    "k" => k,
    "tol" => tol,
    "ncv" => ncv,
    "training_terms" => String.(terms),
    "free_parameters" => String.(parameters),
    "minimum_overlap" => minimum_overlap,
    "require_same_rank" => require_same_rank,
    "anchor_seconds" => anchor_seconds,
    "anchor_q5" => anchor_five.q,
    "anchor_q_training" => anchor_training.q,
    "anchor_factor_training" => anchor_training.factor,
    "workspace_seconds" => problem.workspace_seconds,
    "operator_seconds" => problem.operator_seconds,
    "probe_summaries" => probe_summaries,
    "jacobian_singular_values" => collect(jacobian.singular_values),
    "jacobian_numerical_rank" => jacobian.numerical_rank,
    "jacobian_condition_number" => jacobian.condition_number,
    "maximum_jacobian_condition_number" => maximum_condition_number,
    "julia_version" => string(VERSION),
    "julia_threads" => Threads.nthreads(),
    "fuzzified_version" => string(Base.pkgversion(FuzzifiED)),
    "fuzzified_expected_revision" => source_revision(),
    "project_git_revision" => MottJainED.git_revision(PROJECT_ROOT),
    "so3lver_source_sha256" => sha256_file(joinpath(PROJECT_ROOT, "experimental", "SO3lverED.jl")),
    "search_source_sha256" => sha256_file(joinpath(PROJECT_ROOT, "experimental", "SO3ParameterSearch.jl")),
    "driver_source_sha256" => sha256_file(@__FILE__),
    "config_source_sha256" => sha256_file(config_path),
    "base" => Dict(String(name) => getfield(base, name) for name in fieldnames(Couplings)),
)
isnothing(generator_result) || (summary["generator"] = generator_dict(generator_result))
MottJainED.atomic_toml(joinpath(output, "audit_summary.toml"), summary)

println("audit_passed=$passed")
println("generator_passed=$generator_passed")
println("state_tracking_passed=$(all_tracking_passed[])")
println("ground_branch_passed=$(all_ground_passed[])")
println("jacobian_rank=$(jacobian.numerical_rank)/$(length(parameters))")
println("jacobian_condition_number=$(jacobian.condition_number)")
println("jacobian_passed=$jacobian_passed")
println("audit_result=$(joinpath(output, "audit_summary.toml"))")
passed || println("AUDIT_NOT_CERTIFIED: do not enable the seven-term optimizer")
