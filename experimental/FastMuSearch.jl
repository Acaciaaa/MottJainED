# Historical automatic-search prototype. The active N7 workflow uses FastED
# prepare/sector-array/collect with n7_uf0_165_k20_scout_restart.toml.
module FastMuSearch

using CSV, DataFrames, Dates, LinearAlgebra, MottJainED, TOML
using ..FastED

const ROOT = FastED.PROJECT_ROOT

function options(spec)
    cfg = spec.config["mu_search"]
    opt = (
        lower=Float64(cfg["mu_min"]), upper=Float64(cfg["mu_max"]),
        wide_count=Int(cfg["wide_count"]),
        guard_lower=Float64(cfg["guard_min"]), guard_upper=Float64(cfg["guard_max"]),
        guard_count=Int(cfg["guard_count"]),
        seed=Float64(cfg["seed_mu"]), half_width=Float64(cfg["local_half_width"]),
        local_count=Int(cfg["local_count"]),
        tolerance=Float64(cfg["mu_abs_tol"]), iterations=Int(cfg["max_iterations"]),
        budget=Int(cfg["max_evaluations"]),
        mu_agreement=Float64(cfg["mu_agreement_tol"]),
        q_agreement=Float64(cfg["q_agreement_tol"]),
        neighbor_step=Float64(cfg["neighbor_step"]),
        cold_tolerance=Float64(cfg["cold_repeat_tol"]),
    )
    opt.lower < opt.guard_lower < opt.seed < opt.guard_upper < opt.upper ||
        error("Require mu_min < guard_min < seed_mu < guard_max < mu_max")
    min(opt.wide_count, opt.guard_count, opt.local_count) >= 3 || error("Grids need at least 3 points")
    opt.tolerance > 0 && opt.half_width > 0 && opt.neighbor_step > opt.tolerance ||
        error("Invalid search widths/tolerances")
    opt.budget > 0 && opt.iterations > 0 || error("Invalid search budget")
    opt.mu_agreement > 0 && opt.q_agreement >= 0 && opt.cold_tolerance > 0 ||
        error("Invalid audit tolerances")
    return opt
end

"""Use the established local/Brent search, then independently discover and refine valleys.

The two discovery grids are evaluated irrespective of the local outcome. Invalid scores
remain explicit in records; they are never interpreted as evidence of a high q.
This is a finite-resolution audit, not a proof of a global minimum over real μ.
"""
function search_mu(score_at, opt; on_evaluation=(mu, score, source)->nothing)
    scores = Dict{Float64,Any}()
    sources = Dict{Float64,String}()
    function measured(mu)
        key = MottJainED._mu_key(mu)
        if !haskey(scores, key)
            length(scores) < opt.budget || error("μ evaluation budget exhausted; keep results/cache and review")
            scores[key] = score_at(key)
        end
        return scores[key]
    end
    function record(mu, score, source)
        key = MottJainED._mu_key(mu)
        if !haskey(sources, key)
            sources[key] = String(source)
            on_evaluation(key, score, String(source))
        end
    end
    local_result = MottJainED._continuation_grid_refine(
        measured; center=opt.seed, mu_min=opt.lower, mu_max=opt.upper,
        local_half_width=opt.half_width, local_count=opt.local_count,
        max_expansions=3, abs_tol=opt.tolerance, max_iterations=opt.iterations,
        wide_mode=:always, wide_mu_tol=opt.mu_agreement,
        wide_objective_tol=opt.q_agreement, on_evaluation=record,
    )
    wide_search = MottJainED._new_scalar_search(measured; on_evaluation=record)
    intervals = Any[]
    for (lower, upper, count, prefix) in (
        (opt.lower, opt.upper, opt.wide_count, "wide_discovery"),
        (opt.guard_lower, opt.guard_upper, opt.guard_count, "critical_guard"),
    )
        push!(intervals, MottJainED._refine_grid_interval!(
            wide_search; lower=lower, upper=upper, count=count,
            abs_tol=opt.tolerance, max_iterations=opt.iterations, source_prefix=prefix,
        ))
    end
    wide_result = MottJainED._scalar_search_summary(
        wide_search; mu_min=opt.lower, mu_max=opt.upper, abs_tol=opt.tolerance,
    )
    finite_keys() = sort([mu for (mu, score) in scores if MottJainED._finite_objective(score)])
    function best_key()
        keys = finite_keys()
        isempty(keys) && error("No valid q anywhere in the evaluated search")
        return keys[argmin([scores[mu].objective for mu in keys])]
    end
    # Explicitly evaluate both sides of the winner; retain actual evaluated energies.
    neighbor_converged = false
    best = best_key()
    for attempt in 1:3
        best = best_key()
        lo, hi = max(opt.lower, best-opt.neighbor_step), min(opt.upper, best+opt.neighbor_step)
        interval = MottJainED._refine_grid_interval!(
            wide_search; lower=lo, upper=hi, count=3, abs_tol=opt.tolerance,
            max_iterations=opt.iterations, source_prefix="winner_check_$attempt",
        )
        push!(intervals, interval)
        best = best_key()
        if lo < best < hi && all(MottJainED._finite_objective, [measured(lo), measured(hi)])
            neighbor_converged = true
            break
        end
    end
    invalid = sort([mu for (mu, score) in scores if !MottJainED._finite_objective(score)])
    issues = String[]
    local_result.completed || push!(issues, "local_or_wide_Brent_not_converged")
    local_result.wide_disagreement && push!(issues, "local_and_wide_Brent_disagree")
    all(i.refined_count == i.refinements_converged for i in intervals) ||
        push!(issues, "discovered_valley_not_converged")
    neighbor_converged || push!(issues, "winner_not_bracketed_by_valid_neighbors")
    opt.guard_lower < best < opt.guard_upper || push!(issues, "winner_outside_dense_guard")
    any(opt.guard_lower <= mu <= opt.guard_upper for mu in invalid) &&
        push!(issues, "unscored_point_inside_dense_guard")
    local_result.score === nothing && push!(issues, "local_search_has_no_valid_score")
    wide_result.score === nothing && push!(issues, "independent_grid_has_no_valid_score")
    if local_result.score !== nothing && wide_result.score !== nothing
        abs(local_result.mu-wide_result.mu) <= opt.mu_agreement &&
            abs(local_result.score.objective-wide_result.score.objective) <= opt.q_agreement ||
            push!(issues, "local_and_independent_grid_disagree")
    end
    return (best_mu=best, score=scores[best], scores=scores, sources=sources,
            local_result=local_result, wide_result=wide_result, invalid=invalid, issues=issues,
            neighbor_converged=neighbor_converged, intervals=intervals)
end

"""Solve one sector at a time in this process, reusing disk cache and complete CSVs.

The large matrices/vectors from a completed solve are released before loading the
next sector. No extra worker or idle CPU allocation is held for another sector.
"""
function evaluate_serial(spec, count, mu; force=false)
    for index in 1:count
        started = time()
        @info "Sequential sector" mu index total=count threads=Threads.nthreads()
        try
            FastED.solve_sector(spec, Float64(mu), index; force=force)
        finally
            GC.gc()
        end
        @info "Sector finished" mu index elapsed_seconds=time()-started
    end
    return FastED.collect_mu(spec, Float64(mu))
end

function run(spec)
    opt = options(spec)
    spec.solver.warm_start && error("Set solver.warm_start=false for independent eigensolver starts")
    BLAS.set_num_threads(1)
    FastED.FuzzifiED.NumThreads = Threads.nthreads()
    manifest = FastED.prepare_caches(spec)
    count = length(manifest["sectors"])
    state_count = sum(entry["dimension"] <= spec.solver.dense_cutoff ?
        min(spec.solver.k, entry["dimension"]) : min(spec.solver.k, entry["dimension"]-2)
        for entry in manifest["sectors"])
    mkpath(spec.result_directory)
    audit_path = joinpath(spec.result_directory, "mu_search_audit.toml")
    audit = Dict{String,Any}(
        "complete" => false, "audit_passed" => false, "global_minimum_proven" => false,
        "cache_id" => spec.cache_id, "settings_id" => spec.settings_id,
        "project_git_revision" => spec.identity["project_git_revision"],
        "julia_version" => string(VERSION), "sector_count" => count,
        "execution_mode" => "serial_sectors", "solver_processes" => 1,
        "julia_threads" => Threads.nthreads(),
        "started_at" => string(now()), "search_options" => spec.config["mu_search"],
        "search_source_hash" => FastED.files_digest([
            @__FILE__, joinpath(ROOT, "experimental", "FastED.jl"),
            joinpath(ROOT, "src", "Workflows.jl"),
        ]; root=ROOT),
    )
    MottJainED.atomic_toml(audit_path, audit)
    GC.gc()
    trace = NamedTuple[]
    try
        function record(mu, score, source)
            push!(trace, merge(FastED.score_summary_row(spec, mu, score, state_count),
                              (evaluation=length(trace)+1, source=source, timestamp=string(now()))))
            MottJainED.atomic_csv(joinpath(spec.result_directory, "mu_search_evaluations.csv"), DataFrame(trace))
            @info "N7 μ search" evaluation=length(trace) mu source valid=score.valid q=score.q reason=score.reason
        end
        result = search_mu(mu -> evaluate_serial(spec, count, mu).score, opt; on_evaluation=record)
        # Independent cold re-solve at the selected point, using the same matrices and k.
        cold = evaluate_serial(spec, count, result.best_mu; force=true)
        differences = [abs(getproperty(cold.score, key)-getproperty(result.score, key))
                       for key in (:q, :delta_s, :delta_o)]
        issues = copy(result.issues)
        cold.score.valid && all(isfinite, differences) && maximum(differences) <= opt.cold_tolerance ||
            push!(issues, "independent_cold_repeat_failed")
        MottJainED.atomic_csv(joinpath(spec.result_directory, "cold_repeat.csv"), DataFrame([
            (mu=result.best_mu, q_before=result.score.q, q_after=cold.score.q,
             delta_s_error=differences[2], delta_o_error=differences[3], q_error=differences[1]),
        ]))
        # Preserve the complete scan, including unscored outer points. Do not filter them away
        # to satisfy FastED.validate_collection, and do not release any matrix cache here.
        mus = sort(collect(keys(result.scores)))
        final_spec = merge(spec, (mus=mus,))
        FastED.collect_all(final_spec)
        profile = deepcopy(spec.config)
        profile["fast_ed"]["mus"] = mus
        profile["fast_ed"]["allow_cache_release"] = false
        profile_path = joinpath(spec.result_directory, "evaluated_profile.toml")
        MottJainED.atomic_toml(profile_path, profile)
        summary = CSV.read(joinpath(spec.result_directory, "best_summary.csv"), DataFrame)
        best_mu = Float64(summary.mu[1])
        abs(best_mu-result.best_mu) <= opt.tolerance || push!(issues, "winner_changed_after_cold_repeat")
        merge!(audit, Dict{String,Any}(
            "complete"=>true, "audit_passed"=>isempty(issues), "issues"=>issues,
            "finished_at"=>string(now()), "best_mu"=>best_mu,
            "local_mu"=>result.local_result.local_mu, "grid_mu"=>result.wide_result.mu,
            "local_q"=>result.local_result.local_objective,
            "local_phase_best_mu"=>result.local_result.mu,
            "wide_brent_mu"=>result.local_result.wide_brent_mu,
            "wide_brent_q"=>result.local_result.wide_brent_objective,
            "local_wide_disagreement"=>result.local_result.wide_disagreement,
            "grid_q"=>result.wide_result.score === nothing ? Inf : result.wide_result.score.q,
            "unscored_mus"=>result.invalid, "all_evaluated_scores_valid"=>isempty(result.invalid),
            "evaluated_mu_count"=>length(mus), "cold_repeat_performed"=>true,
            "neighbor_check_passed"=>result.neighbor_converged,
            "cache_released"=>false, "evaluated_profile"=>profile_path,
            "qualification"=>"Finite-grid and local-Brent agreement in the dense guard; unscored outer regions are not excluded minima. No proof of a global or thermodynamic critical point.",
        ))
        MottJainED.atomic_toml(audit_path, audit)
        println("Results: $(spec.result_directory)")
        println("Best sampled mu: $best_mu; local/grid audit passed: $(isempty(issues))")
        isempty(result.invalid) || println("Unscored outer mu values retained for review: $(result.invalid)")
        isempty(issues) || error("μ search needs review; cache retained: $(join(issues, ", "))")
        return audit
    catch err
        audit["error"] = sprint(showerror, err)
        MottJainED.atomic_toml(audit_path, audit)
        rethrow()
    end
end

end
