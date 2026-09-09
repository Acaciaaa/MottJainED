#!/usr/bin/env julia

# 独立的两尺寸调参流程。
#
# N=3,4 只参与宽范围 muc 搜索和 size continuation；最终 matching、零点插值
# 和调试图严格只使用 N=5,6。这个脚本复用 MottJainED 的 FSS 求谱/找 muc 代码，
# 但不修改也不调用原来的 fss-all 拟合与作图入口。

import Pkg

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(PROJECT_ROOT; io=devnull)

using CSV
using CairoMakie
using DataFrames
using FuzzifiED
using LinearAlgebra
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

option_bool(options, key, default=false) = haskey(options, key) ?
    lowercase(options[key]) in ("1", "true", "yes", "on") : default

project_path(path::AbstractString) = isabspath(path) ? normpath(path) :
    normpath(joinpath(PROJECT_ROOT, path))

section(config, name) = get(config, String(name), Dict{String,Any}())
getvalue(values, name, default) = get(values, String(name), default)

function tuning_settings(config)
    tuning = section(config, :two_size_tuning)
    isempty(tuning) && throw(ArgumentError(
        "Profile must define a [two_size_tuning] section",
    ))

    guide_sizes = sort(unique(Int.(getvalue(tuning, :guide_nm_values, [3, 4]))))
    match_sizes = sort(unique(Int.(getvalue(tuning, :match_nm_values, [5, 6]))))
    length(match_sizes) == 2 || throw(ArgumentError(
        "two_size_tuning.match_nm_values must contain exactly two sizes",
    ))
    isempty(guide_sizes) && throw(ArgumentError(
        "two_size_tuning.guide_nm_values cannot be empty",
    ))
    maximum(guide_sizes) < minimum(match_sizes) || throw(ArgumentError(
        "Guide sizes must be smaller than both matching sizes",
    ))
    all_sizes = sort(unique(vcat(guide_sizes, match_sizes)))
    all(diff(all_sizes) .> 0) || throw(ArgumentError("System sizes must be unique"))

    scan_parameter = Symbol(getvalue(tuning, :scan_parameter, "V0"))
    scan_parameter in MottJainED.HAMILTONIAN_FIELDS || throw(ArgumentError(
        "Unknown scan parameter $scan_parameter",
    ))
    scan_parameter == :mu && throw(ArgumentError(
        "two-size tuning cannot scan mu; mu is optimized separately for every size",
    ))
    scan_values = unique(Float64.(getvalue(tuning, :scan_values, Float64[])))
    isempty(scan_values) && throw(ArgumentError(
        "two_size_tuning.scan_values cannot be empty",
    ))

    definition, terms, metric = MottJainED._score_options(
        tuning; default_definition="fss7", default_metric="q",
    )
    base_solver = MottJainED._solver(config)
    solver = MottJainED._with_k(
        base_solver, Int(getvalue(tuning, :k, base_solver.k)),
    )
    anchor_nm = Int(getvalue(tuning, :optimize_anchor_nm, maximum(guide_sizes)))
    anchor_nm == maximum(guide_sizes) || throw(ArgumentError(
        "optimize_anchor_nm must equal the largest guide size so every guide size " *
        "uses a wide anchor search and the first matching size continues from it",
    ))

    fss = FSSSettings(
        nm_values=all_sizes,
        scan_parameter=scan_parameter,
        scan_values=scan_values,
        mu_min=Float64(getvalue(tuning, :mu_min, -0.095)),
        mu_max=Float64(getvalue(tuning, :mu_max, 0.305)),
        mu_count=Int(getvalue(tuning, :mu_count, 41)),
        methods=[:optimize],
        score_definition=definition,
        score_terms=terms,
        score_metric=metric,
        optimize_strategy=:size_continuation,
        optimize_anchor_nm=anchor_nm,
        optimize_local_half_width=Float64(getvalue(
            tuning, :optimize_local_half_width, 0.02,
        )),
        optimize_local_count=Int(getvalue(tuning, :optimize_local_count, 9)),
        optimize_max_expansions=Int(getvalue(tuning, :optimize_max_expansions, 3)),
        optimize_abs_tol=Float64(getvalue(tuning, :optimize_abs_tol, 1e-4)),
        optimize_max_iterations=Int(getvalue(tuning, :optimize_max_iterations, 100)),
        optimize_wide_mode=Symbol(lowercase(replace(
            String(getvalue(tuning, :optimize_wide_mode, "always")), '-' => '_',
        ))),
        optimize_wide_adaptive_nm=Int(getvalue(tuning, :optimize_wide_adaptive_nm, 6)),
        optimize_wide_audit_first=Bool(getvalue(tuning, :optimize_wide_audit_first, false)),
        optimize_wide_jump_tol=Float64(getvalue(tuning, :optimize_wide_jump_tol, 0.03)),
        optimize_wide_mu_tol=Float64(getvalue(tuning, :optimize_wide_mu_tol, 5e-3)),
        optimize_wide_objective_tol=Float64(getvalue(
            tuning, :optimize_wide_objective_tol, 1e-4,
        )),
    )
    return (
        tuning=tuning, guide_sizes=guide_sizes, match_sizes=match_sizes,
        all_sizes=all_sizes, fss=fss, solver=solver,
    )
end

function case_output(config, tuning; output_override=nothing)
    configured_root = output_override === nothing ?
        String(getvalue(tuning, :output_root, "output/two_size_tuning")) :
        String(output_override)
    root = project_path(configured_root)
    label = String(getvalue(tuning, :run_name, "two_size_tuning"))
    return MottJainED._numbered_case_output(root, label, config; feature=:all)
end

function valid_result_rows(path)
    isfile(path) || throw(ArgumentError("Search result not found: $path"))
    data = MottJainED.latest_rows(CSV.read(path, DataFrame))
    required = ["nm1", "scan_value", "muc", "q", "factor", "delta_s", "delta_o"]
    missing_columns = setdiff(required, names(data))
    isempty(missing_columns) || throw(ArgumentError(
        "Search result is missing columns: $(join(missing_columns, ", "))",
    ))
    if "status" in names(data)
        data = filter(row -> String(row.status) == "ok", data)
    end
    if "score_valid" in names(data)
        data = filter(row -> Bool(row.score_valid), data)
    end
    return data
end

function row_for(data, nm1, scan_value)
    selected = filter(row ->
        Int(row.nm1) == nm1 &&
        isapprox(Float64(row.scan_value), scan_value; atol=1e-12, rtol=0.0),
        data,
    )
    nrow(selected) == 1 || return nothing
    return selected[1, :]
end

function build_matching_table(result_path, scan_parameter, scan_values, match_sizes)
    data = valid_result_rows(result_path)
    small_nm, large_nm = match_sizes
    rows = NamedTuple[]
    missing_pairs = String[]
    for scan_value in scan_values
        small = row_for(data, small_nm, scan_value)
        large = row_for(data, large_nm, scan_value)
        if small === nothing || large === nothing
            push!(missing_pairs, "$(scan_parameter)=$(scan_value)")
            continue
        end
        delta_s_small = Float64(small.delta_s)
        delta_s_large = Float64(large.delta_s)
        delta_o_small = Float64(small.delta_o)
        delta_o_large = Float64(large.delta_o)
        q_small = Float64(small.q)
        q_large = Float64(large.q)
        push!(rows, (
            scan_parameter=String(scan_parameter), scan_value=Float64(scan_value),
            small_nm=small_nm, large_nm=large_nm,
            muc_small=Float64(small.muc), muc_large=Float64(large.muc),
            delta_s_small=delta_s_small, delta_s_large=delta_s_large,
            delta_s_drift=delta_s_large - delta_s_small,
            delta_s_mean=(delta_s_small + delta_s_large) / 2,
            delta_o_small=delta_o_small, delta_o_large=delta_o_large,
            delta_o_drift=delta_o_large - delta_o_small,
            delta_o_mean=(delta_o_small + delta_o_large) / 2,
            q_small=q_small, q_large=q_large, q_mean=(q_small + q_large) / 2,
            factor_small=Float64(small.factor), factor_large=Float64(large.factor),
        ))
    end
    isempty(missing_pairs) || @warn "Matching output skipped incomplete size pairs" pairs=missing_pairs
    isempty(rows) && throw(ArgumentError(
        "No complete N=$small_nm/$large_nm pairs are available for matching",
    ))
    return sort(DataFrame(rows), :scan_value)
end

function zero_crossings(matching)
    output = DataFrame(
        observable=String[], lower_scan_value=Float64[], upper_scan_value=Float64[],
        estimated_scan_value=Float64[], estimated_observable=Float64[],
        q_mean=Float64[], interpolation=String[],
    )
    for (observable, drift_column, mean_column) in (
        ("delta_s", :delta_s_drift, :delta_s_mean),
        ("delta_o", :delta_o_drift, :delta_o_mean),
    )
        for index in 1:nrow(matching)
            drift = Float64(matching[index, drift_column])
            if iszero(drift)
                push!(output, (
                    observable=observable,
                    lower_scan_value=Float64(matching.scan_value[index]),
                    upper_scan_value=Float64(matching.scan_value[index]),
                    estimated_scan_value=Float64(matching.scan_value[index]),
                    estimated_observable=Float64(matching[index, mean_column]),
                    q_mean=Float64(matching.q_mean[index]), interpolation="exact",
                ))
            end
        end
        for index in 1:(nrow(matching) - 1)
            left = Float64(matching[index, drift_column])
            right = Float64(matching[index + 1, drift_column])
            (iszero(left) || iszero(right) || signbit(left) == signbit(right)) && continue
            weight = -left / (right - left)
            lerp(column) = (1 - weight) * Float64(matching[index, column]) +
                           weight * Float64(matching[index + 1, column])
            push!(output, (
                observable=observable,
                lower_scan_value=Float64(matching.scan_value[index]),
                upper_scan_value=Float64(matching.scan_value[index + 1]),
                estimated_scan_value=lerp(:scan_value),
                estimated_observable=lerp(mean_column),
                q_mean=lerp(:q_mean), interpolation="linear_sign_change",
            ))
        end
    end
    return output
end

function plot_matching(matching, scan_parameter, output)
    small_nm = first(matching.small_nm)
    large_nm = first(matching.large_nm)
    x = matching.scan_value
    figure = Figure(size=(900, 1050))

    delta_axis = Axis(
        figure[1, 1]; xlabel=String(scan_parameter), ylabel="scaling dimension",
        title="trusted-size matching (guide sizes excluded)",
    )
    scatterlines!(delta_axis, x, matching.delta_s_small; label="DeltaS N=$small_nm")
    scatterlines!(delta_axis, x, matching.delta_s_large; label="DeltaS N=$large_nm")
    scatterlines!(delta_axis, x, matching.delta_o_small; label="DeltaO N=$small_nm")
    scatterlines!(delta_axis, x, matching.delta_o_large; label="DeltaO N=$large_nm")
    axislegend(delta_axis; position=:rt)

    drift_axis = Axis(
        figure[2, 1]; xlabel=String(scan_parameter), ylabel="large - small",
        title="two-size drift",
    )
    hlines!(drift_axis, [0.0]; color=:black, linestyle=:dash)
    scatterlines!(drift_axis, x, matching.delta_s_drift; label="DeltaS drift")
    scatterlines!(drift_axis, x, matching.delta_o_drift; label="DeltaO drift")
    axislegend(drift_axis; position=:rt)

    q_axis = Axis(
        figure[3, 1]; xlabel=String(scan_parameter), ylabel="q",
        title="CFT score guardrail",
    )
    scatterlines!(q_axis, x, matching.q_small; label="q N=$small_nm")
    scatterlines!(q_axis, x, matching.q_large; label="q N=$large_nm")
    axislegend(q_axis; position=:rt)

    save(output, figure)
    return figure
end

function analyze_results(case_directory, settings)
    result_path = joinpath(case_directory, "search", "fss_optimize_results.csv")
    matching = build_matching_table(
        result_path, settings.fss.scan_parameter, settings.fss.scan_values,
        settings.match_sizes,
    )
    crossings = zero_crossings(matching)
    MottJainED.atomic_csv(joinpath(case_directory, "two_size_matching.csv"), matching)
    MottJainED.atomic_csv(
        joinpath(case_directory, "two_size_zero_crossings.csv"), crossings,
    )
    plot_matching(
        matching, settings.fss.scan_parameter,
        joinpath(case_directory, "two_size_matching.png"),
    )
    return matching, crossings
end

function print_plan(config, settings, output_override)
    tuning = settings.tuning
    configured_root = output_override === nothing ?
        String(getvalue(tuning, :output_root, "output/two_size_tuning")) :
        String(output_override)
    label = String(getvalue(tuning, :run_name, "two_size_tuning"))
    println("Two-size tuning plan (no ED has been run)")
    println("  guide sizes       = $(settings.guide_sizes)  (muc navigation only)")
    println("  matching sizes    = $(settings.match_sizes)  (only these enter analysis)")
    println("  scan parameter    = $(settings.fss.scan_parameter)")
    println("  scan values       = $(settings.fss.scan_values)")
    println("  fixed Hamiltonian = $(MottJainED.coupling_namedtuple(MottJainED._couplings(config)))")
    println("  muc wide range    = [$(settings.fss.mu_min), $(settings.fss.mu_max)]")
    println("  wide anchor sizes = N <= $(settings.fss.optimize_anchor_nm)")
    println("  wide challenger   = $(settings.fss.optimize_wide_mode)")
    if settings.fss.optimize_wide_mode == :adaptive
        println("  adaptive from     = N >= $(settings.fss.optimize_wide_adaptive_nm)")
        println("  first-point audit = $(settings.fss.optimize_wide_audit_first)")
        println("  guard tolerances  = jump $(settings.fss.optimize_wide_jump_tol), " *
                "mu $(settings.fss.optimize_wide_mu_tol), " *
                "objective $(settings.fss.optimize_wide_objective_tol)")
    end
    println("  k per sector      = $(settings.solver.k)")
    println("  output prefix     = $(joinpath(project_path(configured_root), label * "_XX"))")
    return nothing
end

function main(args=ARGS)
    options = parse_options(args)
    config_path = project_path(get(options, "config", "config/my_run.toml"))
    profile_path = project_path(get(
        options, "profile", "config/two_size_tuning/v0_stage1.toml",
    ))
    config = MottJainED.load_config(config_path; override=profile_path)
    settings = tuning_settings(config)
    output_override = get(options, "output", nothing)
    if option_bool(options, "plan", false)
        print_plan(config, settings, output_override)
        return 0
    end

    solver_section = section(config, :solver)
    BLAS.set_num_threads(Int(getvalue(solver_section, :blas_threads, 1)))
    configured_threads = Int(getvalue(solver_section, :fuzzified_threads, 0))
    FuzzifiED.NumThreads = configured_threads > 0 ? configured_threads : Threads.nthreads()

    output = case_output(config, settings.tuning; output_override=output_override)
    logging_state = MottJainED.start_task_logging(output; filename="run.log")
    println("Detailed progress: $(logging_state.path)")
    try
        MottJainED.write_resolved_config(
            output, config; base_config=config_path, override_config=profile_path,
        )
        MottJainED.write_run_metadata(
            output; command="two-size-tuning", config_path=profile_path,
        )
        MottJainED.atomic_toml(joinpath(output, "two_size_manifest.toml"), Dict(
            "command" => "two-size-tuning",
            "guide_nm_values" => settings.guide_sizes,
            "match_nm_values" => settings.match_sizes,
            "scan_parameter" => String(settings.fss.scan_parameter),
            "scan_values" => settings.fss.scan_values,
            "optimize_wide_mode" => String(settings.fss.optimize_wide_mode),
            "optimize_wide_adaptive_nm" => settings.fss.optimize_wide_adaptive_nm,
            "optimize_wide_audit_first" => settings.fss.optimize_wide_audit_first,
            "optimize_wide_jump_tol" => settings.fss.optimize_wide_jump_tol,
            "optimize_wide_mu_tol" => settings.fss.optimize_wide_mu_tol,
            "optimize_wide_objective_tol" => settings.fss.optimize_wide_objective_tol,
            "raw_search_directory" => "search",
            "matching_file" => "two_size_matching.csv",
            "zero_crossings_file" => "two_size_zero_crossings.csv",
            "plot_file" => "two_size_matching.png",
        ))
        unless_analysis_only = !option_bool(options, "analyze-only", false)
        if unless_analysis_only
            MottJainED.run_fss_scan(
                MottJainED._couplings(config), settings.fss, settings.solver;
                output=joinpath(output, "search"),
                force=option_bool(options, "force", false),
            )
        end
        matching, crossings = analyze_results(output, settings)
        println("Two-size matching rows: $(nrow(matching))")
        println("Sign-change estimates: $(nrow(crossings))")
        println("Results: $output")
    finally
        MottJainED.stop_task_logging(logging_state)
    end
    return 0
end

exit(main())
