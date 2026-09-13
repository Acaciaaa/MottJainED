module FastEDPipeline

using CSV
using DataFrames
using Dates
using LinearAlgebra
using MottJainED
using SHA
using TOML
using ..FastED

const STATE_SCHEMA_VERSION = 1
const FIXED_TERMS = [:ds_s, :j, :curlj, :dj_rank1, :t_rank1]
const REFERENCE_STATES = (
    (label="S", l2=0, c2=0, raw_rank=2),
    (label="dS", l2=2, c2=0, raw_rank=1),
    (label="J", l2=2, c2=3, raw_rank=1),
    (label="curlJ", l2=2, c2=3, raw_rank=3),
    (label="dJ", l2=6, c2=3, raw_rank=1),
    (label="T", l2=6, c2=0, raw_rank=1),
    (label="O", l2=0, c2=3, raw_rank=1),
)

key(mu::Real) = round(Float64(mu); digits=14)

function pipeline_options(spec)
    cfg = get(spec.config, "pipeline", nothing)
    cfg isa AbstractDict || throw(ArgumentError("Profile must define [pipeline]"))
    name = FastED.safe_label(String(cfg["name"]))
    root = FastED.project_path(String(get(cfg, "state_root", "output/fast_ed/pipelines")))
    archive_root = FastED.project_path(String(get(cfg, "archive_root", "output/fast_ed/archives")))
    fine_count = Int(get(cfg, "fine_count", 7))
    isodd(fine_count) && fine_count >= 3 || throw(ArgumentError("pipeline.fine_count must be odd and >=3"))
    opt = (
        name=name,
        directory=joinpath(root, name),
        archive_root=archive_root,
        scan_parameter=String(cfg["scan_parameter"]),
        scan_value=Float64(cfg["scan_value"]),
        n5_mu=Float64(cfg["n5_mu"]),
        n6_mu=Float64(cfg["n6_mu"]),
        n5_q=Float64(cfg["n5_q"]),
        n6_q=Float64(cfg["n6_q"]),
        n5_factor=Float64(cfg["n5_factor"]),
        n6_factor=Float64(cfg["n6_factor"]),
        n5_delta_s=Float64(cfg["n5_delta_s"]),
        n6_delta_s=Float64(cfg["n6_delta_s"]),
        n5_delta_o=Float64(cfg["n5_delta_o"]),
        n6_delta_o=Float64(cfg["n6_delta_o"]),
        scout_step=Float64(get(cfg, "scout_step", 0.005)),
        fine_step=Float64(get(cfg, "fine_step", 0.00025)),
        fit_snap_step=Float64(get(cfg, "fit_snap_step", 0.000125)),
        fine_count=fine_count,
        max_adaptive_rounds=Int(get(cfg, "max_adaptive_rounds", 3)),
        max_total_mus=Int(get(cfg, "max_total_mus", 16)),
        mu_min=Float64(get(cfg, "mu_min", minimum(spec.mus)-0.05)),
        mu_max=Float64(get(cfg, "mu_max", maximum(spec.mus)+0.05)),
        accept_neighbor_distance=Float64(get(cfg, "accept_neighbor_distance", 1.01*Float64(get(cfg, "fine_step", 0.00025)))),
        q_fit_improvement_tol=Float64(get(cfg, "q_fit_improvement_tol", 1e-4)),
        max_quantum_error=Float64(get(cfg, "max_quantum_error", 5e-4)),
        max_copy_split=Float64(get(cfg, "max_copy_split", spec.solver.degeneracy_tol)),
        max_factor_relative_jump=Float64(get(cfg, "max_factor_relative_jump", 0.05)),
        max_gap_relative_jump=Float64(get(cfg, "max_gap_relative_jump", 0.05)),
        guard_mus=Float64.(get(cfg, "guard_mus", Float64[])),
        prepare_cpus=Int(get(cfg, "prepare_cpus", 8)),
        prepare_threads=Int(get(cfg, "prepare_threads", 8)),
        solve_cpus=Int(get(cfg, "solve_cpus", 8)),
        solve_threads=Int(get(cfg, "solve_threads", 8)),
        max_concurrent=Int(get(cfg, "max_concurrent", 4)),
        control_cpus=Int(get(cfg, "control_cpus", 1)),
        auto_release=Bool(get(cfg, "auto_release", false)),
    )
    spec.nm1 == 7 || throw(ArgumentError("The bounded pilot is restricted to N=7"))
    spec.solver.k == 20 && !spec.solver.warm_start || throw(ArgumentError(
        "The bounded pilot requires the validated cold-start k=20 solver",
    ))
    spec.terms == FIXED_TERMS && spec.metric == :q || throw(ArgumentError(
        "The bounded pilot requires the audited five fixed q terms",
    ))
    opt.scan_parameter in ("Uf0", "Vf0") || throw(ArgumentError(
        "pipeline.scan_parameter must be Uf0 or Vf0",
    ))
    getproperty(spec.couplings, Symbol(opt.scan_parameter)) == opt.scan_value ||
        throw(ArgumentError("pipeline.scan_value does not match the Hamiltonian"))
    all(isfinite, (opt.n5_mu, opt.n6_mu, opt.n5_q, opt.n6_q, opt.n5_factor,
                   opt.n6_factor, opt.n5_delta_s, opt.n6_delta_s,
                   opt.n5_delta_o, opt.n6_delta_o)) || throw(ArgumentError(
        "pipeline N5/N6 reference values must be finite",
    ))
    continuation = opt.n6_mu + 0.5*(opt.n6_mu-opt.n5_mu)
    expected_scout = continuation .+ opt.scout_step .* collect(-2:2)
    length(spec.mus) == 5 && all(isapprox.(spec.mus, expected_scout; atol=1e-14, rtol=0)) ||
        throw(ArgumentError("Initial scout must be the five-point damped N5/N6 continuation grid"))
    opt.mu_min < minimum(spec.mus) <= maximum(spec.mus) < opt.mu_max ||
        throw(ArgumentError("Initial scout must lie strictly inside pipeline.mu_min/mu_max"))
    opt.scout_step > 0 && opt.fine_step > 0 && opt.fit_snap_step > 0 ||
        throw(ArgumentError("Pipeline spacings must be positive"))
    opt.max_adaptive_rounds >= 1 || throw(ArgumentError("max_adaptive_rounds must be positive"))
    opt.max_total_mus >= length(spec.mus)+opt.fine_count ||
        throw(ArgumentError("max_total_mus is too small for scout plus refinement"))
    opt.prepare_cpus >= opt.prepare_threads >= 1 || throw(ArgumentError("prepare CPU/thread settings are invalid"))
    opt.solve_cpus >= opt.solve_threads >= 1 || throw(ArgumentError("solve CPU/thread settings are invalid"))
    opt.max_concurrent >= 1 && opt.control_cpus >= 1 || throw(ArgumentError("pipeline concurrency settings are invalid"))
    all(mu -> mu in spec.mus, opt.guard_mus) ||
        throw(ArgumentError("Every pipeline.guard_mus value must be part of the initial scout"))
    return opt
end

state_path(opt) = joinpath(opt.directory, "pipeline_state.toml")
profiles_directory(opt) = joinpath(opt.directory, "profiles")
final_directory(opt) = joinpath(opt.directory, "final")

function read_state(opt)
    path = state_path(opt)
    isfile(path) || throw(ArgumentError("Pipeline state not found: $path"))
    state = TOML.parsefile(path)
    Int(get(state, "schema_version", -1)) == STATE_SCHEMA_VERSION ||
        throw(ArgumentError("Unsupported pipeline-state schema"))
    return state
end

write_state(opt, state) = MottJainED.atomic_toml(state_path(opt), state)

function same_problem(left, right)
    return left.nm1 == right.nm1 && left.cache_id == right.cache_id &&
           left.settings_id == right.settings_id &&
           FastED.coupling_dict(left.couplings; include_mu=false) ==
               FastED.coupling_dict(right.couplings; include_mu=false) &&
           FastED.solver_dict(left.solver) == FastED.solver_dict(right.solver) &&
           left.terms == right.terms && left.metric == right.metric
end

function initialize(config_path::AbstractString)
    spec = FastED.load_spec(config_path)
    opt = pipeline_options(spec)
    mkpath(profiles_directory(opt)); mkpath(final_directory(opt)); mkpath(opt.archive_root)
    path = state_path(opt)
    if isfile(path)
        state = read_state(opt)
        String(state["cache_id"]) == spec.cache_id || throw(ArgumentError("Existing pipeline state has another cache identity"))
        String(state["settings_id"]) == spec.settings_id || throw(ArgumentError("Existing pipeline state has another solver identity"))
        return state
    end
    state = Dict{String,Any}(
        "schema_version" => STATE_SCHEMA_VERSION,
        "name" => opt.name,
        "base_config" => abspath(config_path),
        "cache_id" => spec.cache_id,
        "settings_id" => spec.settings_id,
        "created_at" => string(now()),
        "updated_at" => string(now()),
        "action" => "ready",
        "stage" => "scout",
        "adaptive_round" => 0,
        "profiles" => [abspath(config_path)],
        "submission_log" => String[],
        "accepted" => false,
        "complete" => false,
        "cache_released" => false,
    )
    write_state(opt, state)
    return state
end

function claim_launch(config_path::AbstractString)
    spec = FastED.load_spec(config_path); opt = pipeline_options(spec)
    state = initialize(config_path)
    String(state["action"]) == "ready" || throw(ArgumentError(
        "Pipeline is already claimed or running (action=$(state["action"]))",
    ))
    state["action"] = "launching"
    state["updated_at"] = string(now())
    write_state(opt, state)
    return state
end

function record_submission(config_path::AbstractString; stage::AbstractString,
                           solve_job::AbstractString, controller_job::AbstractString,
                           prepare_job::AbstractString="")
    spec = FastED.load_spec(config_path); opt = pipeline_options(spec); state = read_state(opt)
    String(state["action"]) in ("launching", "solve") || throw(ArgumentError(
        "Cannot record submission while action=$(state["action"])",
    ))
    entry = "$(now()) stage=$(FastED.safe_label(stage))"
    isempty(prepare_job) || (entry *= " prepare=$(FastED.safe_label(prepare_job))")
    entry *= " solve=$(FastED.safe_label(solve_job)) controller=$(FastED.safe_label(controller_job))"
    submission_log = String[String(value) for value in get(state, "submission_log", Any[])]
    push!(submission_log, entry)
    state["submission_log"] = submission_log
    state["action"] = "awaiting_results"
    state["updated_at"] = string(now())
    write_state(opt, state)
    return state
end

function parabolic_fit(xs, ys)
    length(xs) == length(ys) == 3 || throw(ArgumentError("A three-point fit is required"))
    matrix = hcat(Float64.(xs).^2, Float64.(xs), ones(3))
    a, b, c = matrix \ Float64.(ys)
    isfinite(a) && isfinite(b) && isfinite(c) && a > 0 || return nothing
    vertex = -b/(2a)
    value = a*vertex^2+b*vertex+c
    return (vertex=vertex, value=value, curvature=a)
end

snap(mu::Real, step::Real) = key(round(Float64(mu)/Float64(step))*Float64(step))

function score_triplet(rows, index)
    1 < index < length(rows) || return nothing
    xs = [rows[i].mu for i in index-1:index+1]
    qs = [rows[i].q for i in index-1:index+1]
    fit = parabolic_fit(xs, qs.^2)
    fit === nothing && return nothing
    minimum(xs) < fit.vertex < maximum(xs) || return nothing
    predicted_q = sqrt(max(fit.value, 0.0))
    return (xs=xs, qs=qs, vertex=fit.vertex, predicted_q=predicted_q,
            improvement=qs[2]-predicted_q)
end

function fresh_mus(candidates, existing, opt)
    seen = Set(key.(existing)); output = Float64[]
    for candidate in candidates
        mu = key(candidate)
        opt.mu_min < mu < opt.mu_max || continue
        mu in seen && continue
        push!(output, mu); push!(seen, mu)
    end
    return sort!(output)
end

"""Bounded state-machine decision; it never performs or submits an ED solve."""
function decide_next(rows, latest_mus, stage::AbstractString, round::Integer, opt)
    isempty(rows) && return (action="review", kind="none", mus=Float64[], reason="no_scored_points")
    ordered = sort(collect(rows); by=row -> row.mu)
    all(row -> row.valid && isfinite(row.q) && row.factor > 0, ordered) ||
        return (action="review", kind="none", mus=Float64[], reason="invalid_score")
    best_index = argmin(getproperty.(ordered, :q)); best = ordered[best_index]
    existing = getproperty.(ordered, :mu)
    latest_set = Set(key.(latest_mus))
    latest = filter(row -> key(row.mu) in latest_set, ordered)
    latest_index = findfirst(row -> key(row.mu) == key(best.mu), latest)

    if stage in ("scout", "scout_extension")
        if best_index in (1, length(ordered))
            round < opt.max_adaptive_rounds ||
                return (action="review", kind="none", mus=Float64[], reason="scout_boundary_round_limit")
            direction = best_index == 1 ? -1.0 : 1.0
            candidates = [best.mu+direction*opt.scout_step, best.mu+direction*2*opt.scout_step]
            mus = fresh_mus(candidates, existing, opt)
            isempty(mus) && return (action="review", kind="none", mus=mus, reason="scout_boundary_hard_limit")
            return (action="solve", kind="scout_extension", mus=mus, reason="extend_boundary_scout")
        end
        fit = score_triplet(ordered, best_index)
        center = fit === nothing ? best.mu : fit.vertex
        center = snap(clamp(center, ordered[best_index-1].mu, ordered[best_index+1].mu), opt.fit_snap_step)
        half = (opt.fine_count-1)÷2
        mus = fresh_mus([center+i*opt.fine_step for i in -half:half], existing, opt)
        isempty(mus) && return (action="review", kind="none", mus=mus, reason="no_new_refinement_points")
        return (action="solve", kind="refine", mus=mus, reason="refine_bracketed_scout_minimum")
    end

    triplet = score_triplet(ordered, best_index)
    if triplet !== nothing
        left, right = ordered[best_index-1], ordered[best_index+1]
        close = best.mu-left.mu <= opt.accept_neighbor_distance &&
                right.mu-best.mu <= opt.accept_neighbor_distance
        rises = left.q > best.q && right.q > best.q
        factor_jump = max(abs(left.factor-best.factor), abs(right.factor-best.factor))/best.factor
        gap_jump = maximum(vcat(
            abs.((left.raw_gaps .- best.raw_gaps) ./ best.raw_gaps),
            abs.((right.raw_gaps .- best.raw_gaps) ./ best.raw_gaps),
        ))
        continuous = factor_jump <= opt.max_factor_relative_jump && gap_jump <= opt.max_gap_relative_jump
        if close && rises && triplet.improvement <= opt.q_fit_improvement_tol && continuous
            return (action="accept", kind="final", mus=Float64[], reason="measured_minimum_bracketed",
                    best=best, triplet=triplet, factor_jump=factor_jump, gap_jump=gap_jump)
        end
    end

    round < opt.max_adaptive_rounds ||
        return (action="review", kind="none", mus=Float64[], reason="adaptive_round_limit")
    candidates = Float64[]
    if latest_index !== nothing && latest_index in (1, length(latest)) && length(latest) >= 3
        direction = latest_index == 1 ? -1.0 : 1.0
        push!(candidates, best.mu+direction*opt.fine_step)
        near = latest_index == 1 ? latest[1:3] : latest[end-2:end]
        fit = parabolic_fit(getproperty.(near, :mu), getproperty.(near, :q).^2)
        if fit !== nothing && minimum(getproperty.(near, :mu)) < fit.vertex < maximum(getproperty.(near, :mu))
            push!(candidates, snap(fit.vertex, opt.fit_snap_step))
        end
    else
        if best_index == 1 || best.mu-ordered[best_index-1].mu > opt.accept_neighbor_distance
            push!(candidates, best.mu-opt.fine_step)
        elseif ordered[best_index-1].q <= best.q
            push!(candidates, best.mu-opt.fine_step)
        end
        if best_index == length(ordered) || ordered[best_index+1].mu-best.mu > opt.accept_neighbor_distance
            push!(candidates, best.mu+opt.fine_step)
        elseif ordered[best_index+1].q <= best.q
            push!(candidates, best.mu+opt.fine_step)
        end
        triplet === nothing || push!(candidates, snap(triplet.vertex, opt.fit_snap_step))
    end
    mus = fresh_mus(candidates, existing, opt)
    remaining = opt.max_total_mus-length(existing)
    remaining > 0 || return (action="review", kind="none", mus=Float64[], reason="mu_budget_exhausted")
    length(mus) > remaining && (mus = mus[1:remaining])
    isempty(mus) && return (action="review", kind="none", mus=mus, reason="no_new_followup_points")
    return (action="solve", kind="followup", mus=mus, reason="tighten_or_bracket_measured_minimum")
end

function generated_profile(base, opt, kind, round, mus)
    config = deepcopy(base.config)
    config["fast_ed"]["mus"] = Float64.(mus)
    config["fast_ed"]["run_name"] = "$(opt.name)_$(kind)_$(round)"
    config["fast_ed"]["allow_cache_release"] = false
    config["pipeline_generated"] = Dict{String,Any}(
        "kind" => kind, "adaptive_round" => round,
        "generated_at" => string(now()), "base_config" => base.config_path,
    )
    return config
end

function validate_and_score(profile_paths, opt)
    base = FastED.load_spec(first(profile_paths)); manifest = FastED.load_cache_manifest(base)
    summaries = NamedTuple[]; relations = NamedTuple[]; selected = NamedTuple[]; files = NamedTuple[]
    seen = Set{Float64}(); issues = String[]
    max_quantum_error = 0.0; max_copy_split = 0.0; max_selected_rank = 0
    spectrum_paths = Dict{Float64,String}()
    for profile_path in profile_paths
        spec = FastED.load_spec(profile_path)
        same_problem(base, spec) || throw(ArgumentError("Pipeline profile changes cache, solver, Hamiltonian, or score"))
        FastED.collect_all(spec)
        for mu in spec.mus
            mu_key = key(mu); mu_key in seen && throw(ArgumentError("Duplicate evaluated mu=$mu_key")); push!(seen, mu_key)
            result = FastED.collect_mu(spec, mu)
            score = result.score
            raw_gaps = score.valid ? Float64.(score.raw_gaps) : fill(NaN, 5)
            push!(summaries, merge(FastED.score_summary_row(spec, mu, score, length(result.states)),
                                   (valid=score.valid, raw_gaps=raw_gaps, profile=spec.config_path)))
            score.valid && append!(relations, FastED.relation_rows(spec, mu, score))
            spectrum_paths[mu_key] = joinpath(FastED.mu_directory(spec, mu), "merged_spectrum.csv")
            data = result.merged
            max_quantum_error = max(max_quantum_error,
                maximum(abs.(data.l2-round.(data.l2))), maximum(abs.(data.c2-round.(data.c2))))
            ground = data[argmin(data.energy), :]
            abs(ground.l2) <= spec.solver.quantum_tol && abs(ground.c2) <= spec.solver.quantum_tol ||
                push!(issues, "mu=$mu ground_not_singlet")
            for entry in manifest["sectors"]
                z, r = Int(entry["z"]), Int(entry["r"])
                path = FastED.sector_result_path(spec, mu, z, r)
                FastED.result_is_current(path, spec, mu, z, r) || push!(issues, "mu=$mu invalid_sector_$(z)_$(r)")
                meta_path = replace(path, r"\.csv$" => ".toml")
                isfile(meta_path) || (push!(issues, "mu=$mu missing_metadata_$(z)_$(r)"); continue)
                meta = TOML.parsefile(meta_path)
                expected_job = MottJainED.stable_id("fast-ed-sector-v$(FastED.RESULT_SCHEMA_VERSION)",
                    spec.cache_id, spec.settings_id, Float64(mu), z, r)
                get(meta, "job_id", "") == expected_job || push!(issues, "mu=$mu wrong_job_identity_$(z)_$(r)")
                Int(get(meta, "julia_threads", -1)) == opt.solve_threads || push!(issues, "mu=$mu wrong_threads_$(z)_$(r)")
                get(meta, "solver", Dict()) == FastED.solver_dict(spec.solver) || push!(issues, "mu=$mu wrong_solver_$(z)_$(r)")
                Int(get(meta, "sector_dimension", -1)) == Int(entry["dimension"]) || push!(issues, "mu=$mu wrong_dimension_$(z)_$(r)")
                for file in (path, meta_path)
                    push!(files, (file=relpath(file, FastED.PROJECT_ROOT), sha256=bytes2hex(sha256(read(file))), bytes=filesize(file)))
                end
            end
            for ref in REFERENCE_STATES
                states = filter(row -> abs(row.l2-ref.l2) <= spec.solver.quantum_tol &&
                                       abs(row.c2-ref.c2) <= spec.solver.quantum_tol, data)
                nrow(states) >= ref.raw_rank || (push!(issues, "mu=$mu missing_$(ref.label)"); continue)
                row = states[ref.raw_rank, :]
                push!(selected, (mu=Float64(mu), label=ref.label, l2=ref.l2, c2=ref.c2,
                    raw_rank=ref.raw_rank, z=Int(row.z), r=Int(row.r), sector_rank=Int(row.rank),
                    energy=Float64(row.energy), gap=Float64(row.energy-minimum(data.energy)),
                    dimension=score.valid ? Float64((row.energy-minimum(data.energy))/score.factor) : NaN))
                max_selected_rank = max(max_selected_rank, Int(row.rank))
            end
            for (l2, c2) in ((2,3), (6,3), (0,3))
                states = filter(row -> abs(row.l2-l2) <= spec.solver.quantum_tol &&
                                       abs(row.c2-c2) <= spec.solver.quantum_tol, data)
                nrow(states) >= 2 || (push!(issues, "mu=$mu missing_symmetry_copy_$(l2)_$(c2)"); continue)
                max_copy_split = max(max_copy_split, abs(Float64(states.energy[2]-states.energy[1])))
            end
        end
    end
    max_quantum_error <= opt.max_quantum_error || push!(issues, "quantum_number_error")
    max_copy_split <= opt.max_copy_split || push!(issues, "multiplet_copy_split")
    rows = sort!(summaries; by=row -> row.mu)
    return (base=base, manifest=manifest, rows=rows, relations=relations, selected=selected,
            files=files, issues=unique(issues), max_quantum_error=max_quantum_error,
            max_copy_split=max_copy_split, max_selected_rank=max_selected_rank,
            spectrum_paths=spectrum_paths)
end

function write_combined(opt, analysis)
    summary_columns = [:cache_id, :settings_id, :nm1, :mu, :k, :score_valid,
        :definition, :score_terms, :metric, :objective, :q, :cost, :factor,
        :delta_s, :delta_o, :state_count, :reason]
    summary = select(DataFrame(analysis.rows), summary_columns)
    MottJainED.atomic_csv(joinpath(final_directory(opt), "combined_scan_summary.csv"), summary)
    MottJainED.atomic_csv(joinpath(final_directory(opt), "relations.csv"), DataFrame(analysis.relations))
    MottJainED.atomic_csv(joinpath(final_directory(opt), "selected_states.csv"), DataFrame(analysis.selected))
    MottJainED.atomic_csv(joinpath(final_directory(opt), "source_hashes.csv"), DataFrame(analysis.files))
    return summary
end

function completed_state(state, archive_path::AbstractString)
    output = deepcopy(state)
    output["accepted"] = true
    output["complete"] = true
    output["action"] = "complete"
    output["stage"] = "final"
    output["bundle_path"] = abspath(archive_path)
    output["updated_at"] = string(now())
    for field in ("next_mus", "next_profile", "next_tasks", "bundle_list")
        delete!(output, field)
    end
    return output
end

function finalise!(state, opt, analysis, decision)
    summary = write_combined(opt, analysis)
    best = decision.best; best_mu = key(best.mu)
    best_index = findfirst(==(best_mu), key.(Float64.(summary.mu)))
    MottJainED.atomic_csv(joinpath(final_directory(opt), "best_summary.csv"), summary[best_index:best_index, :])
    relation_data = DataFrame(analysis.relations)
    MottJainED.atomic_csv(joinpath(final_directory(opt), "best_relations.csv"),
        filter(row -> key(row.mu) == best_mu, relation_data))
    selected_data = DataFrame(analysis.selected)
    MottJainED.atomic_csv(joinpath(final_directory(opt), "best_selected_states.csv"),
        filter(row -> key(row.mu) == best_mu, selected_data))
    cp(analysis.spectrum_paths[best_mu], joinpath(final_directory(opt), "best_spectrum.csv"); force=true)
    fss = DataFrame([
        (scan_parameter=opt.scan_parameter, scan_value=opt.scan_value, nm1=5,
         mu=opt.n5_mu, q=opt.n5_q, factor=opt.n5_factor,
         delta_s=opt.n5_delta_s, delta_o=opt.n5_delta_o, source="trusted_n56"),
        (scan_parameter=opt.scan_parameter, scan_value=opt.scan_value, nm1=6,
         mu=opt.n6_mu, q=opt.n6_q, factor=opt.n6_factor,
         delta_s=opt.n6_delta_s, delta_o=opt.n6_delta_o, source="trusted_n56"),
        (scan_parameter=opt.scan_parameter, scan_value=opt.scan_value, nm1=7,
         mu=best.mu, q=best.q, factor=best.factor,
         delta_s=best.delta_s, delta_o=best.delta_o, source="bounded_n7_pipeline"),
    ])
    MottJainED.atomic_csv(joinpath(final_directory(opt), "fss_n567.csv"), fss)
    audit = Dict{String,Any}(
        "complete" => true, "accepted" => true, "global_minimum_proven" => false,
        "thermodynamic_criticality_proven" => false, "cache_id" => analysis.base.cache_id,
        "settings_id" => analysis.base.settings_id, "evaluated_mu_count" => length(analysis.rows),
        "best_sampled_mu" => best.mu, "best_sampled_q" => best.q,
        "best_factor" => best.factor, "best_delta_s" => best.delta_s, "best_delta_o" => best.delta_o,
        "neighbor_mus" => decision.triplet.xs, "neighbor_q" => decision.triplet.qs,
        "quadratic_mu_estimate" => decision.triplet.vertex,
        "quadratic_predicted_q_improvement" => decision.triplet.improvement,
        "max_quantum_error" => analysis.max_quantum_error,
        "max_lowest_multiplet_split" => analysis.max_copy_split,
        "max_selected_sector_rank" => analysis.max_selected_rank,
        "factor_neighbor_relative_jump" => decision.factor_jump,
        "gap_neighbor_relative_jump" => decision.gap_jump,
        "issues" => analysis.issues, "fixed_raw_rank_scoring" => true,
        "scan_parameter" => opt.scan_parameter, "scan_value" => opt.scan_value,
        "logs_required" => false, "sacct_required" => false,
        "qualification" => "Accepted measured finite-grid minimum with close valid neighbors. Not a proof of a global or thermodynamic critical point.",
        "finished_at" => string(now()),
    )
    MottJainED.atomic_toml(joinpath(final_directory(opt), "pipeline_audit.toml"), audit)
    archive_path = joinpath(opt.archive_root, "$(opt.name)_final.tar.gz")
    state["accepted"] = true; state["complete"] = false; state["action"] = "bundle"
    state["stage"] = "final"; state["best_mu"] = best.mu; state["best_q"] = best.q
    state["bundle_path"] = archive_path
    for field in ("next_mus", "next_profile", "next_tasks")
        delete!(state, field)
    end
    packaged_state = completed_state(state, archive_path)
    packaged_state["packaged_snapshot"] = true
    packaged_state["bundle_sha256_sidecar"] = basename(archive_path) * ".sha256"
    MottJainED.atomic_toml(joinpath(final_directory(opt), "pipeline_state.toml"), packaged_state)
    paths = String[relpath(String(state["base_config"]), FastED.PROJECT_ROOT)]
    append!(paths, relpath.(String.(state["profiles"]), Ref(FastED.PROJECT_ROOT)))
    append!(paths, relpath.([FastED.load_spec(path).result_directory for path in state["profiles"]], Ref(FastED.PROJECT_ROOT)))
    push!(paths, relpath(FastED.cache_manifest_path(analysis.base), FastED.PROJECT_ROOT))
    push!(paths, relpath(final_directory(opt), FastED.PROJECT_ROOT))
    paths = sort!(unique(paths))
    all(path -> !startswith(path, "..") && !isabspath(path) && !occursin('\n', path), paths) ||
        error("Bundle paths must stay inside the project root")
    list_path = joinpath(opt.directory, "bundle_file_list.txt")
    temporary = list_path * ".tmp-$(getpid())"
    try
        open(temporary, "w") do io
            foreach(path -> println(io, path), paths)
        end
        mv(temporary, list_path; force=true)
    finally
        isfile(temporary) && rm(temporary; force=true)
    end
    state["bundle_list"] = list_path
    state["updated_at"] = string(now())
    write_state(opt, state)
    return state
end

function advance(config_path::AbstractString)
    base = FastED.load_spec(config_path); opt = pipeline_options(base); state = read_state(opt)
    String(state["action"]) == "awaiting_results" || throw(ArgumentError(
        "advance requires action=awaiting_results, got $(state["action"])",
    ))
    profiles = String.(state["profiles"])
    analysis = validate_and_score(profiles, opt)
    write_combined(opt, analysis)
    if !isempty(analysis.issues)
        state["action"] = "review"; state["review_reason"] = join(analysis.issues, ",")
        state["updated_at"] = string(now()); write_state(opt, state); return state
    end
    latest = FastED.load_spec(last(profiles))
    decision = decide_next(analysis.rows, latest.mus, String(state["stage"]),
                           Int(state["adaptive_round"]), opt)
    state["last_decision_reason"] = decision.reason
    if decision.action == "accept"
        return finalise!(state, opt, analysis, decision)
    elseif decision.action == "review"
        state["action"] = "review"; state["review_reason"] = decision.reason
        state["updated_at"] = string(now()); write_state(opt, state); return state
    end
    total = length(analysis.rows)+length(decision.mus)
    total <= opt.max_total_mus || begin
        state["action"] = "review"; state["review_reason"] = "mu_budget_exhausted"
        state["updated_at"] = string(now()); write_state(opt, state); return state
    end
    next_round = Int(state["adaptive_round"])+1
    config = generated_profile(base, opt, decision.kind, next_round, decision.mus)
    profile_path = joinpath(profiles_directory(opt), "round_$(next_round)_$(decision.kind).toml")
    MottJainED.atomic_toml(profile_path, config)
    generated = FastED.load_spec(profile_path)
    same_problem(base, generated) || error("Generated profile changed the physical problem")
    append!(state["profiles"], [profile_path])
    state["stage"] = decision.kind; state["adaptive_round"] = next_round
    state["action"] = "solve"; state["next_profile"] = profile_path
    state["next_tasks"] = length(decision.mus)*length(FastED.SECTOR_ORDER)
    state["next_mus"] = decision.mus; state["updated_at"] = string(now())
    write_state(opt, state)
    return state
end

function mark_bundled(config_path::AbstractString, archive_path::AbstractString,
                      expected_sha::AbstractString)
    base = FastED.load_spec(config_path); opt = pipeline_options(base); state = read_state(opt)
    String(state["action"]) == "bundle" || throw(ArgumentError("Pipeline is not awaiting a bundle"))
    real_archive = abspath(archive_path)
    real_archive == abspath(String(state["bundle_path"])) || throw(ArgumentError("Unexpected bundle path"))
    isfile(real_archive) || throw(ArgumentError("Bundle does not exist: $real_archive"))
    actual = bytes2hex(sha256(read(real_archive)))
    lowercase(expected_sha) == actual || throw(ArgumentError("Bundle SHA-256 mismatch"))
    finished = completed_state(state, real_archive)
    finished["bundle_sha256"] = actual
    finished["bundle_bytes"] = filesize(real_archive)
    write_state(opt, finished)
    return finished
end

function action_fields(config_path::AbstractString)
    base = FastED.load_spec(config_path); opt = pipeline_options(base); state = read_state(opt)
    action = String(state["action"])
    if action == "ready"
        return (action, base.config_path, string(length(base.mus)*length(FastED.SECTOR_ORDER)), state_path(opt))
    elseif action == "solve"
        return (action, String(state["next_profile"]), string(state["next_tasks"]), state_path(opt))
    elseif action == "bundle"
        return (action, String(state["bundle_path"]), String(state["bundle_list"]), state_path(opt))
    end
    return (action, "-", "-", state_path(opt))
end

function resource_fields(config_path::AbstractString)
    opt = pipeline_options(FastED.load_spec(config_path))
    return (opt.prepare_cpus, opt.prepare_threads, opt.solve_cpus, opt.solve_threads,
            opt.max_concurrent, opt.control_cpus)
end

function mark_review(config_path::AbstractString, reason::AbstractString)
    base = FastED.load_spec(config_path); opt = pipeline_options(base); state = read_state(opt)
    state["action"] = "review"
    state["review_reason"] = FastED.safe_label(reason)
    state["updated_at"] = string(now())
    write_state(opt, state)
    return state
end

function reset_launch(config_path::AbstractString)
    base = FastED.load_spec(config_path); opt = pipeline_options(base); state = read_state(opt)
    String(state["action"]) == "launching" || throw(ArgumentError(
        "Only an unrecorded launch can be reset",
    ))
    state["action"] = "ready"
    state["updated_at"] = string(now())
    write_state(opt, state)
    return state
end

export pipeline_options, initialize, claim_launch, record_submission, parabolic_fit,
       decide_next, validate_and_score, advance, mark_bundled, action_fields,
       resource_fields, mark_review, reset_launch, state_path

end
