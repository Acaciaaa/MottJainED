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

finite_maximum(values) = isempty(values) ? Inf : maximum(abs, values)

function main()
    options = parse_options(ARGS)
    config_path = abspath(get(
        options, "config",
        joinpath(PROJECT_ROOT, "config", "so3lver", "n6_linked_uf_v0_extension.toml"),
    ))
    isfile(config_path) || throw(ArgumentError("configuration not found: $config_path"))
    config = TOML.parsefile(config_path)
    run_config = config["run"]
    score_config = config["score"]
    mu_config = config["mu"]
    scan_config = config["scan"]
    continuation_config = config["continuation"]
    robustness_config = config["robustness"]

    nm = Int(run_config["nm"])
    allow_test_size = lowercase(get(options, "allow-test-size", "false")) == "true"
    (nm == 6 || allow_test_size) || error(
        "this driver is intentionally restricted to N=6; " *
        "--allow-test-size=true is only for an end-to-end smoke test",
    )
    k = Int(run_config["k"])
    k >= 3 || error("the stable six-term score needs at least three levels per block")
    tol = Float64(run_config["tol"])
    ncv = Int(run_config["ncv"])
    ncv > k || error("ncv must exceed k")
    base = dict_couplings(config["base"])

    terms = Symbol.(score_config["terms"])
    terms == collect(CFT_STABLE_SIX_TERMS) || error(
        "the linked scan must use exactly the audited stable six relations",
    )
    worst_weight = Float64(score_config["worst_residual_weight"])
    minimum_overlap = Float64(score_config["minimum_overlap"])
    penalty = Float64(score_config["penalty"])
    maximum_q5_ratio = Float64(score_config["maximum_q5_ratio"])
    0 <= minimum_overlap <= 1 || error("minimum_overlap must lie in [0,1]")
    penalty > 0 || error("penalty must be positive")

    mu_lower = Float64(mu_config["lower"])
    mu_upper = Float64(mu_config["upper"])
    mu_lower < mu_upper || error("mu bounds must have positive width")
    mu_grid_count = Int(mu_config["grid_count"])
    mu_grid_count >= 3 || error("mu grid needs at least three points")
    mu_refine_basins = Int(mu_config["refine_basins"])
    mu_refine_iterations = Int(mu_config["refine_iterations"])
    mu_absolute_tolerance = Float64(mu_config["absolute_tolerance"])
    mu_boundary_tolerance = Float64(mu_config["boundary_tolerance"])

    u0_over_uf = Float64(scan_config["u0_over_uf"])
    isfinite(u0_over_uf) && u0_over_uf > 0 || error(
        "scan.u0_over_uf must be finite and positive",
    )
    isapprox(base.U0, u0_over_uf * base.Uf; atol=1e-12, rtol=0) || error(
        "the reference point must obey U0 = scan.u0_over_uf * Uf",
    )
    uf_values = Float64.(scan_config["uf_values"])
    uf0_values = Float64.(scan_config["uf0_values"])
    vf0_values = Float64.(scan_config["vf0_values"])
    v0_values = Float64.(scan_config["v0_values"])
    linked_grid = linked_uf_grid_values(
        uf_values, uf0_values, vf0_values, v0_values; u0_over_uf,
    )
    expected_profile_count = Int(scan_config["expected_profile_count"])
    top_candidate_count = Int(scan_config["top_candidate_count"])
    top_candidate_count > 0 || error("top_candidate_count must be positive")

    parameter_names = (:Uf, :Uf0, :Vf0, :V0)
    allowed_config = scan_config["allowed_bounds"]
    allowed_lower = Float64[
        allowed_config[String(parameter)][1] for parameter in parameter_names
    ]
    allowed_upper = Float64[
        allowed_config[String(parameter)][2] for parameter in parameter_names
    ]
    all(allowed_lower .< allowed_upper) || error("allowed bounds must have positive width")

    continuation_steps = Dict(
        Symbol(name) => Float64(value)
        for (name, value) in continuation_config["steps"]
    )
    Set(keys(continuation_steps)) == Set((:Uf, :Uf0, :Vf0, :V0, :mu)) || error(
        "continuation steps must contain exactly Uf,Uf0,Vf0,V0,mu",
    )
    all(value > 0 for value in values(continuation_steps)) || error(
        "continuation steps must be positive",
    )
    maximum_continuation_steps = Int(continuation_config["maximum_steps"])
    maximum_continuation_steps > 0 || error("maximum continuation steps must be positive")

    robustness_steps = Dict(
        Symbol(name) => Float64(value)
        for (name, value) in robustness_config["steps"]
    )
    Set(keys(robustness_steps)) == Set(parameter_names) || error(
        "robustness steps must contain exactly Uf,Uf0,Vf0,V0",
    )
    all(value > 0 for value in values(robustness_steps)) || error(
        "robustness steps must be positive",
    )
    robustness_tolerance = Float64(robustness_config["objective_tolerance"])
    robustness_tolerance >= 0 || error("robustness tolerance must be nonnegative")

    function profile_couplings(values; mu=base.mu)
        length(values) == length(parameter_names) || throw(DimensionMismatch(
            "a linked scan profile must contain Uf,Uf0,Vf0,V0",
        ))
        uf, uf0, vf0, v0 = Float64.(values)
        return Couplings(
            Uf=uf,
            U0=u0_over_uf * uf,
            Uf0=uf0,
            Vf=base.Vf,
            Vf0=vf0,
            V0=v0,
            t=base.t,
            mu=Float64(mu),
        )
    end

    scan_points = NamedTuple[]
    seen_points = Set{NTuple{4,Float64}}()
    function add_scan_point!(source, values)
        physical = Float64.(values)
        all((allowed_lower .<= physical) .& (physical .<= allowed_upper)) || error(
            "scan point $source lies outside allowed bounds: $physical",
        )
        key = Tuple(round.(physical; digits=12))
        key in seen_points && error("duplicate scan point $source at $physical")
        push!(seen_points, key)
        push!(scan_points, (source=String(source), values=physical))
        return nothing
    end

    for extra in get(scan_config, "extra_profiles", Any[])
        uf = Float64(extra["Uf"])
        add_scan_point!(String(extra["name"]), Float64[
            uf, extra["Uf0"], extra["Vf0"], extra["V0"],
        ])
    end
    for (index, point) in enumerate(linked_grid)
        add_scan_point!("grid_$(lpad(index, 3, '0'))", Float64[
            point.Uf, point.Uf0, point.Vf0, point.V0,
        ])
    end
    length(scan_points) == expected_profile_count || error(
        "expected $expected_profile_count unique profiles, built $(length(scan_points))",
    )

    output = abspath(get(
        options, "output", joinpath(PROJECT_ROOT, String(run_config["output"])),
    ))
    validate_only = lowercase(get(options, "validate-only", "false")) == "true"
    so3_path = joinpath(PROJECT_ROOT, "experimental", "SO3lverED.jl")
    search_path = joinpath(PROJECT_ROOT, "experimental", "SO3ParameterSearch.jl")
    driver_path = abspath(@__FILE__)
    signature = MottJainED.stable_id(
        "so3lver-n6-linked-uf-v0-extension-v1", nm, k, tol, ncv,
        String.(terms), u0_over_uf, scan_points, allowed_lower, allowed_upper,
        mu_lower, mu_upper, mu_grid_count, mu_refine_basins,
        sha256_file(config_path), sha256_file(so3_path), sha256_file(search_path),
        sha256_file(driver_path),
    )

    println("N=$nm deterministic linked-Uf/U0 SO(3)lver extension")
    println("profiles=$(length(scan_points)) outer_parameters=Uf(linked U0),Uf0,Vf0,V0")
    println("mu is fully reprofiled on [$mu_lower,$mu_upper]")
    println("output=$output")
    println("Julia $(VERSION), threads=$(Threads.nthreads()), FuzzifiED $(Base.pkgversion(FuzzifiED))")
    flush(stdout)
    if validate_only
        println("VALIDATION_OK: deterministic grid, bounds, linkage, and signatures loaded")
        return nothing
    end

    point_path = joinpath(output, "point_evaluations.csv")
    residual_path = joinpath(output, "point_residuals.csv")
    tracking_path = joinpath(output, "state_tracking.csv")
    profile_path = joinpath(output, "profile_evaluations.csv")
    robustness_path = joinpath(output, "robustness_neighbors.csv")
    best_path = joinpath(output, "best.toml")

    if isfile(best_path)
        previous_best = TOML.parsefile(best_path)
        String(get(previous_best, "signature", "")) == signature || error(
            "existing best.toml belongs to another scan; choose a new output directory",
        )
        println("ALREADY_COMPLETE: existing best.toml matches this scan")
        return nothing
    end

    function validate_existing_trace(path)
        isfile(path) || return DataFrame()
        previous = CSV.read(path, DataFrame)
        "signature" in names(previous) || error(
            "existing trace $(basename(path)) has no signature",
        )
        all(String(value) == signature for value in previous.signature) || error(
            "existing trace $(basename(path)) belongs to another scan",
        )
        return previous
    end

    previous_points = validate_existing_trace(point_path)
    previous_profiles = validate_existing_trace(profile_path)
    validate_existing_trace(residual_path)
    validate_existing_trace(tracking_path)
    evaluation = Ref(nrow(previous_points))

    mkpath(output)
    cp(config_path, joinpath(output, "search_config.toml"); force=true)
    manifest_rows = [(
        source=point.source,
        Uf=point.values[1],
        U0=u0_over_uf * point.values[1],
        Uf0=point.values[2],
        Vf0=point.values[3],
        V0=point.values[4],
    ) for point in scan_points]
    MottJainED.atomic_csv(joinpath(output, "scan_manifest.csv"), DataFrame(manifest_rows))

    problem = build_cft_problem(nm, base; disp_std=true)
    anchor = solve_cft_blocks!(
        problem.hamiltonians, base; k, tol, ncv,
        warm_vectors=Dict{Tuple{Symbol,Int},Vector{Float64}}(), disp_std=true,
    )
    anchor_six = score_cft_blocks(anchor.energies; terms)
    anchor_five = score_cft_blocks(anchor.energies; terms=CFT_SCORE_TERMS)
    anchor_six.ground_is_singlet_l0 || error(
        "anchor is not on the singlet L=0 ground branch",
    )
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
            assessment=nothing,
            minimum_overlap=0.0,
            reason="continuation_requires_$(count)_steps",
        )
        previous_vectors = anchor.vectors
        warm = anchor_warm()
        path_minimum = 1.0
        final_assessment = nothing
        for index in 1:count
            fraction = index / count
            couplings = Couplings(; (
                name => (1 - fraction) * getfield(base, name) +
                        fraction * getfield(target, name)
                for name in fieldnames(Couplings)
            )...)
            solved = solve_cft_blocks!(
                problem.hamiltonians, couplings; k, tol, ncv,
                warm_vectors=warm, disp_std=false,
            )
            assessment = characterize(solved; reference_vectors=previous_vectors)
            path_minimum = min(path_minimum, assessment.min_overlap)
            if !assessment.valid
                return (
                    assessment=assessment,
                    minimum_overlap=path_minimum,
                    reason="continuation_step_$(index)_$(assessment.reason)",
                )
            end
            previous_vectors = solved.vectors
            final_assessment = assessment
        end
        return (
            assessment=final_assessment,
            minimum_overlap=path_minimum,
            reason="ok",
        )
    end

    function record_point!(source, profile_id, couplings, assessment, objective,
                           valid, reason, identity_mode, min_overlap, seconds)
        evaluation[] += 1
        six = isnothing(assessment) ? nothing : assessment.six
        five = isnothing(assessment) ? nothing : assessment.five
        tracking = isnothing(assessment) ? nothing : assessment.tracking
        MottJainED.append_csv(point_path, (
            signature=signature,
            evaluation=evaluation[],
            profile_id=String(profile_id),
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
            Uf=couplings.Uf,
            Uf0=couplings.Uf0,
            U0=couplings.U0,
            Vf=couplings.Vf,
            Vf0=couplings.Vf0,
            V0=couplings.V0,
            t=couplings.t,
            mu=couplings.mu,
        ))
        if !isnothing(six)
            for index in eachindex(six.terms)
                MottJainED.append_csv(residual_path, (
                    signature=signature,
                    evaluation=evaluation[],
                    profile_id=String(profile_id),
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
                    profile_id=String(profile_id),
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

    function evaluate_target(couplings, warm; source, profile_id)
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
            objective = valid ? assessment.raw_objective : penalty + assessment.raw_objective
            record_point!(
                source, profile_id, couplings, assessment, objective, valid, reason,
                identity_mode, minimum_path_overlap, time() - started,
            )
            println("point=$(evaluation[]) source=$source profile=$profile_id " *
                    "objective=$objective valid=$valid mu=$(couplings.mu)")
            flush(stdout)
            compact_assessment = merge(assessment, (solved=nothing,))
            return (
                objective=objective,
                valid=valid,
                reason=reason,
                assessment=compact_assessment,
                couplings=couplings,
                identity_mode=identity_mode,
                minimum_overlap=minimum_path_overlap,
            )
        catch err
            reason = replace(sprint(showerror, err), '\n' => ' ')
            record_point!(
                source, profile_id, couplings, assessment, penalty, false, reason,
                identity_mode, minimum_path_overlap, time() - started,
            )
            @error "N=$nm point evaluation failed" profile_id source couplings exception=(err, catch_backtrace())
            return (
                objective=penalty,
                valid=false,
                reason=reason,
                assessment=nothing,
                couplings=couplings,
                identity_mode=identity_mode,
                minimum_overlap=minimum_path_overlap,
            )
        end
    end

    profile_id(values) = MottJainED.stable_id(
        signature, round.(Float64.(values); digits=11),
    )

    search_profiles = Dict{String,Any}()
    if nrow(previous_profiles) > 0
        for row in eachrow(previous_profiles)
            source = String(row.source)
            (source == "best_recheck" || startswith(source, "robustness_")) && continue
            values = Float64[row.Uf, row.Uf0, row.Vf0, row.V0]
            search_profiles[String(row.profile_id)] = (
                id=String(row.profile_id),
                values=values,
                valid=Bool(row.valid),
                reason=String(row.reason),
                best_mu=Float64(row.best_mu),
                objective=Float64(row.objective),
                q6=Float64(row.q6),
                q5=Float64(row.q5),
                worst_residual=Float64(row.worst_residual),
                minimum_overlap=Float64(row.minimum_expected_overlap),
                mu_at_boundary=Bool(row.mu_at_boundary),
                best_result=nothing,
                source=source,
            )
        end
        println("resumed_complete_profiles=$(length(search_profiles))")
    end

    function record_profile!(summary, source, seconds)
        couplings = profile_couplings(summary.values; mu=summary.best_mu)
        MottJainED.append_csv(profile_path, (
            signature=signature,
            profile_id=summary.id,
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
            Uf=couplings.Uf,
            Uf0=couplings.Uf0,
            U0=couplings.U0,
            Vf=couplings.Vf,
            Vf0=couplings.Vf0,
            V0=couplings.V0,
            t=couplings.t,
        ))
        return nothing
    end

    function profile_mu(values; source="profile", force=false, store_search=true)
        physical = Float64.(values)
        id = profile_id(physical)
        !force && haskey(search_profiles, id) && return search_profiles[id]
        started = time()
        fixed = profile_couplings(physical)
        warm = anchor_warm()
        mu_results = Dict{String,Any}()
        function evaluate_mu(mu, mu_source)
            bounded_mu = Float64(mu)
            key = string(round(bounded_mu; digits=12))
            haskey(mu_results, key) && return mu_results[key]
            couplings = Couplings(; (
                name => (name == :mu ? bounded_mu : getfield(fixed, name))
                for name in fieldnames(Couplings)
            )...)
            result = evaluate_target(
                couplings, warm;
                source="$(source)_$(mu_source)", profile_id=id,
            )
            mu_results[key] = result
            return result
        end

        grid = collect(range(mu_lower, mu_upper; length=mu_grid_count))
        grid_results = [evaluate_mu(mu, "grid") for mu in grid]
        brackets = mu_refinement_brackets(
            grid,
            getproperty.(grid_results, :objective),
            getproperty.(grid_results, :valid);
            maximum_count=mu_refine_basins,
        )
        for (basin_index, bracket) in enumerate(brackets)
            objective(mu) = evaluate_mu(mu, "refine_$(basin_index)").objective
            try
                result = optimize(
                    objective, bracket.lower, bracket.upper, Brent();
                    abs_tol=mu_absolute_tolerance,
                    iterations=mu_refine_iterations,
                )
                evaluate_mu(
                    Optim.minimizer(result), "refine_$(basin_index)_minimum",
                )
            catch err
                @warn "mu basin refinement failed" id basin_index exception=(err, catch_backtrace())
            end
        end

        valid_results = filter(
            result -> result.valid && isfinite(result.objective),
            collect(Base.values(mu_results)),
        )
        if isempty(valid_results)
            summary = (
                id=id,
                values=physical,
                valid=false,
                reason="no_valid_mu",
                best_mu=NaN,
                objective=penalty,
                q6=Inf,
                q5=Inf,
                worst_residual=Inf,
                minimum_overlap=0.0,
                mu_at_boundary=true,
                best_result=nothing,
                source=String(source),
            )
        else
            best_mu_result = valid_results[argmin(getproperty.(valid_results, :objective))]
            six = best_mu_result.assessment.six
            summary = (
                id=id,
                values=physical,
                valid=true,
                reason="ok",
                best_mu=best_mu_result.couplings.mu,
                objective=best_mu_result.objective,
                q6=six.q,
                q5=best_mu_result.assessment.five.q,
                worst_residual=finite_maximum(six.residuals),
                minimum_overlap=best_mu_result.minimum_overlap,
                mu_at_boundary=min(
                    best_mu_result.couplings.mu - mu_lower,
                    mu_upper - best_mu_result.couplings.mu,
                ) <= mu_boundary_tolerance,
                best_result=best_mu_result,
                source=String(source),
            )
        end
        record_profile!(summary, source, time() - started)
        store_search && (search_profiles[id] = summary)
        println("profile=$id source=$source objective=$(summary.objective) " *
                "mu=$(summary.best_mu) valid=$(summary.valid)")
        flush(stdout)
        return summary
    end

    for point in scan_points
        profile_mu(point.values; source=point.source)
    end
    length(search_profiles) == expected_profile_count || error(
        "completed $(length(search_profiles)) search profiles; " *
        "expected $expected_profile_count",
    )

    valid_profiles = sort(
        filter(
            row -> row.valid && isfinite(row.objective),
            collect(values(search_profiles)),
        );
        by=row -> row.objective,
    )
    isempty(valid_profiles) && error("linked scan produced no valid profile")
    best_recorded = first(valid_profiles)
    best_profile = profile_mu(
        best_recorded.values;
        source="best_recheck", force=true, store_search=false,
    )
    best_profile.valid || error("best grid candidate failed its full mu recheck")

    neighbor_rows = NamedTuple[]
    robustness_tolerance_passed = true
    for (index, parameter) in enumerate(parameter_names)
        for direction in (-1, 1)
            neighbor_values = copy(best_profile.values)
            neighbor_values[index] += direction * robustness_steps[parameter]
            inside = all((allowed_lower .<= neighbor_values) .&
                         (neighbor_values .<= allowed_upper))
            if !inside
                robustness_tolerance_passed = false
                push!(neighbor_rows, (
                    parameter=String(parameter),
                    direction=direction,
                    value=neighbor_values[index],
                    valid=false,
                    objective=Inf,
                    delta_objective=Inf,
                    best_mu=NaN,
                ))
                continue
            end
            neighbor = profile_mu(
                neighbor_values;
                source="robustness_$(parameter)_$(direction > 0 ? "plus" : "minus")",
                force=true, store_search=false,
            )
            acceptable = neighbor.valid &&
                         neighbor.objective >= best_profile.objective - robustness_tolerance
            robustness_tolerance_passed &= acceptable
            push!(neighbor_rows, (
                parameter=String(parameter),
                direction=direction,
                value=neighbor_values[index],
                valid=neighbor.valid,
                objective=neighbor.objective,
                delta_objective=neighbor.objective - best_profile.objective,
                best_mu=neighbor.best_mu,
            ))
        end
    end
    MottJainED.atomic_csv(robustness_path, DataFrame(neighbor_rows))

    top_rows = NamedTuple[]
    for (rank, candidate) in enumerate(first(
        valid_profiles, min(top_candidate_count, length(valid_profiles)),
    ))
        couplings = profile_couplings(candidate.values; mu=candidate.best_mu)
        push!(top_rows, (
            rank=rank,
            source=candidate.source,
            profile_id=candidate.id,
            objective=candidate.objective,
            q6=candidate.q6,
            q5=candidate.q5,
            worst_residual=candidate.worst_residual,
            minimum_expected_overlap=candidate.minimum_overlap,
            mu_at_boundary=candidate.mu_at_boundary,
            Uf=couplings.Uf,
            U0=couplings.U0,
            Uf0=couplings.Uf0,
            Vf0=couplings.Vf0,
            V0=couplings.V0,
            mu=candidate.best_mu,
        ))
    end
    MottJainED.atomic_csv(joinpath(output, "top_candidates.csv"), DataFrame(top_rows))

    grid_lower = Float64[
        minimum(uf_values), minimum(uf0_values), minimum(vf0_values), minimum(v0_values),
    ]
    grid_upper = Float64[
        maximum(uf_values), maximum(uf0_values), maximum(vf0_values), maximum(v0_values),
    ]
    boundary_parameters = String[]
    for (index, parameter) in enumerate(parameter_names)
        value = best_profile.values[index]
        tolerance = 1e-10 * max(1.0, abs(grid_lower[index]), abs(grid_upper[index]))
        if value <= grid_lower[index] + tolerance || value >= grid_upper[index] - tolerance
            push!(boundary_parameters, String(parameter))
        end
    end

    q5_guard_passed = best_profile.q5 <= maximum_q5_ratio * anchor_five.q
    scan_conclusive = best_profile.valid && !best_profile.mu_at_boundary && q5_guard_passed &&
                      robustness_tolerance_passed && isempty(boundary_parameters)
    best_result = best_profile.best_result
    six = best_result.assessment.six
    best_couplings = profile_couplings(best_profile.values; mu=best_profile.best_mu)
    best_dict = Dict{String,Any}(
        "completed_at" => string(now()),
        "diagnostic_complete" => true,
        "scan_conclusive" => scan_conclusive,
        "signature" => signature,
        "nm1" => nm,
        "design" => "deterministic_linked_uf_cartesian_grid",
        "u0_over_uf" => u0_over_uf,
        "objective" => best_profile.objective,
        "q6" => best_profile.q6,
        "q5" => best_profile.q5,
        "anchor_q6" => anchor_six.q,
        "anchor_q5" => anchor_five.q,
        "q5_guard_passed" => q5_guard_passed,
        "mu_at_boundary" => best_profile.mu_at_boundary,
        "grid_boundary_parameters" => boundary_parameters,
        "robustness_tolerance_passed" => robustness_tolerance_passed,
        "factor" => six.factor,
        "delta_s" => six.delta_s,
        "delta_o" => six.delta_o,
        "worst_residual" => best_profile.worst_residual,
        "minimum_expected_overlap" => best_profile.minimum_overlap,
        "training_terms" => String.(terms),
        "labels" => six.labels,
        "residuals" => six.residuals,
        "profile_count" => length(search_profiles),
        "point_evaluation_count" => evaluation[],
        "parameters" => Dict(
            "Uf" => best_couplings.Uf,
            "Uf0" => best_couplings.Uf0,
            "U0" => best_couplings.U0,
            "Vf" => best_couplings.Vf,
            "Vf0" => best_couplings.Vf0,
            "V0" => best_couplings.V0,
            "t" => best_couplings.t,
            "mu" => best_couplings.mu,
        ),
        "config_source_sha256" => sha256_file(config_path),
        "so3lver_source_sha256" => sha256_file(so3_path),
        "search_source_sha256" => sha256_file(search_path),
        "driver_source_sha256" => sha256_file(driver_path),
        "project_git_revision" => MottJainED.git_revision(PROJECT_ROOT),
        "fuzzified_version" => string(Base.pkgversion(FuzzifiED)),
    )
    MottJainED.atomic_toml(best_path, best_dict)

    println("diagnostic_complete=true")
    println("scan_conclusive=$scan_conclusive")
    println("q5_guard_passed=$q5_guard_passed")
    println("robustness_tolerance_passed=$robustness_tolerance_passed")
    println("grid_boundary_parameters=$(join(boundary_parameters, ','))")
    println("best_objective=$(best_profile.objective) q6=$(best_profile.q6) q5=$(best_profile.q5)")
    println("best_parameters=Uf=$(best_couplings.Uf),U0=$(best_couplings.U0)," *
            "Uf0=$(best_couplings.Uf0),Vf0=$(best_couplings.Vf0)," *
            "V0=$(best_couplings.V0),mu=$(best_couplings.mu)")
    println("scan_result=$best_path")
    return nothing
end

main()
