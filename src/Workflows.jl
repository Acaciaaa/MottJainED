# 本文件是项目的“业务流程层”：把 Model/Spectrum/CFT/Storage 中的基础函数
# 组合成用户真正会运行的 spectrum、gap、density、critical、optimization、FSS。

function _summary_row(job_id, status, nm1, mu, couplings, settings, score; error="")
    return (
        job_id=String(job_id), status=String(status), timestamp=string(now()),
        nm1=Int(nm1), mu=Float64(mu), q=score.q, factor=score.factor,
        delta_s=score.delta_s, delta_o=score.delta_o,
        score_valid=score.valid, score_reason=score.reason,
        k=settings.k, error=String(error), coupling_namedtuple(couplings)...,
    )
end

"""
对一个系统大小和一列 μ 求低能谱，并保存每点的原始谱与任务摘要。

输出 `summary.csv` 以及 `spectra/<job_id>.csv`；`force=false` 时跳过已成功
完成的 job。只有确实需要本征向量时才设 `keep_vectors=true`，否则文件很大。
这个基础功能不再自动套用尚未确认的 CFT tower 标准。
"""
function run_spectrum_scan(
    nm1::Int,
    mus,
    couplings::Couplings=Couplings(),
    settings::SolverSettings=SolverSettings();
    output::AbstractString=joinpath(PACKAGE_ROOT, "output", "spectrum"),
    force::Bool=false,
    keep_vectors::Bool=false,
)
    output = ensure_output(output)
    write_run_metadata(output; command="spectrum")
    summary_path = joinpath(output, "summary.csv")
    completed = force ? Set{String}() : completed_job_ids(summary_path)
    # model/cache 在整列 μ 上只建立一次；循环内部只重新组合 H0+μNf 并求谱。
    model = build_model(nm1=nm1)
    cache = prepare_spectrum(model, couplings, settings)

    for mu in Float64.(collect(mus))
        job_id = stable_id("spectrum-v1", nm1, mu, coupling_vector(couplings), settings.k)
        job_id in completed && continue
        @info "spectrum" nm1 mu job_id
        try
            states = solve_spectrum(cache, mu; keep_vectors=keep_vectors)
            score = CFTScore(q=NaN, reason="not evaluated by spectrum")
            table = spectrum_dataframe(
                states; mu=mu, nm1=nm1, quantum_tol=settings.quantum_tol,
            )
            insertcols!(table, 1, :job_id => fill(job_id, nrow(table)))
            atomic_csv(joinpath(output, "spectra", "$job_id.csv"), table)
            if keep_vectors
                jldsave(joinpath(output, "spectra", "$job_id.jld2"); states=states)
            end
            append_csv(summary_path, _summary_row(job_id, "ok", nm1, mu, couplings, settings, score))
        catch err
            @error "spectrum job failed" nm1 mu exception=(err, catch_backtrace())
            append_csv(summary_path, _summary_row(
                job_id, "error", nm1, mu, couplings, settings, CFTScore();
                error=sanitize_error(err),
            ))
        end
    end
    return latest_rows(CSV.read(summary_path, DataFrame))
end

function _legacy_singlet_gap(states::Vector{SpectrumState}, quantum_tol::Real)
    isempty(states) && return NaN
    singlets = _states_in_sector(states, 0, 0; quantum_tol=quantum_tol)
    isempty(singlets) && return NaN
    reference = length(singlets) >= 2 ? singlets[2] : singlets[1]
    return reference.energy - first(states).energy
end

function _legacy_j_gap(states::Vector{SpectrumState}, quantum_tol::Real)
    isempty(states) && return NaN
    currents = _states_in_sector(states, 2, 3; quantum_tol=quantum_tol)
    isempty(currents) && return NaN
    return currents[1].energy - first(states).energy
end

function _plot_gap_curves(data::DataFrame, column::Symbol, ylabel, title, path)
    valid = filter(row -> row.status == "ok" && isfinite(row[column]), data)
    nrow(valid) > 0 || return nothing
    fig = Figure(size=(650, 650))
    axis = Axis(
        fig[1, 1]; xlabel="μ", ylabel=ylabel, title=title, aspect=1,
        xgridcolor=(:gray, 0.20), ygridcolor=(:gray, 0.20),
        xminorgridcolor=(:gray, 0.12),
    )
    ylims!(axis, 0.0, 1.0)
    axis.xminorticks = IntervalsBetween(10)
    axis.xminorgridvisible = true
    for nm1 in sort(unique(valid.nm1))
        sub = sort(filter(row -> row.nm1 == nm1, valid), :mu)
        lines!(axis, sub.mu, sub[!, column]; linewidth=2, label="nm1 = $nm1")
        scatter!(axis, sub.mu, sub[!, column]; markersize=10)
    end
    axislegend(axis; position=:rt)
    save(path, fig)
    return fig
end

"""
扫描多个系统大小和 μ，同时提取 scalar gap 与 J gap。

scalar 沿用旧 `ES_mu.jl`：原始 `(0,0)` 列表第二项减基态；J 沿用旧
CFT 分析：最低 `(2,3)` 态减基态。两者都乘 `sqrt(nm1)` 后画图。
"""
function run_gap_scan(
    nm_values::AbstractVector{<:Integer},
    mus,
    couplings::Couplings=Couplings(),
    settings::SolverSettings=SolverSettings();
    output::AbstractString=joinpath(PACKAGE_ROOT, "output", "gap"),
    force::Bool=false,
)
    output = ensure_output(output)
    write_run_metadata(output; command="gap")
    path = joinpath(output, "gap_results.csv")
    completed = force ? Set{String}() : completed_job_ids(path)
    for nm1 in Int.(nm_values)
        model = build_model(nm1=nm1)
        cache = prepare_spectrum(model, couplings, settings)
        for mu in Float64.(collect(mus))
            job_id = stable_id("gap-v1", nm1, mu, coupling_vector(couplings), settings.k)
            job_id in completed && continue
            try
                states = solve_spectrum(cache, mu)
                # 与旧 ES_mu.jl 相同：在保留对称性副本的原始 (0,0) 列表中取第二项。
                scalar_gap = _legacy_singlet_gap(states, settings.quantum_tol)
                j_gap = _legacy_j_gap(states, settings.quantum_tol)
                append_csv(path, (
                    job_id=job_id, status="ok", timestamp=string(now()), nm1=nm1,
                    mu=mu, scalar_gap=scalar_gap,
                    scaled_scalar_gap=scalar_gap*sqrt(nm1),
                    j_gap=j_gap, scaled_j_gap=j_gap*sqrt(nm1), error="",
                    coupling_namedtuple(couplings)...,
                ))
            catch err
                append_csv(path, (
                    job_id=job_id, status="error", timestamp=string(now()), nm1=nm1,
                    mu=mu, scalar_gap=NaN, scaled_scalar_gap=NaN,
                    j_gap=NaN, scaled_j_gap=NaN, error=sanitize_error(err),
                    coupling_namedtuple(couplings)...,
                ))
            end
        end
    end
    data = latest_rows(CSV.read(path, DataFrame))
    _plot_gap_curves(
        data, :scaled_scalar_gap, "ΔE_S · √nm1", "Scalar gap vs μ",
        joinpath(output, "scalar_gap.png"),
    )
    _plot_gap_curves(
        data, :scaled_j_gap, "ΔE_J · √nm1", "J gap vs μ",
        joinpath(output, "j_gap.png"),
    )
    return data
end

function _plot_density_curve(data::DataFrame, path)
    valid = filter(row -> row.status == "ok", data)
    nrow(valid) > 0 || return nothing
    sort!(valid, :mu)
    fig = Figure(size=(650, 650))
    axis = Axis(
        fig[1, 1]; xlabel="μ", ylabel="⟨n₁₂₃⟩",
        title="light fermion density vs μ", aspect=1,
    )
    ylims!(axis, 0.0, 3.0)
    lines!(axis, valid.mu, valid.nf_per_orbital; linewidth=2)
    scatter!(axis, valid.mu, valid.nf_per_orbital)
    hlines!(axis, [0.0, 1.0, 2.0]; color=:gray, linestyle=:dash, linewidth=1)
    save(path, fig)
    return fig
end

"""
扫描 μ 并计算基态中的 charge-1/charge-3 平均粒子数及每轨道密度。

这里只求全局基态并保留它的向量；由总电荷约束
`Nf+3N0=3nm1` 推出 N0。
"""
function run_density_scan(
    nm1::Int,
    mus,
    couplings::Couplings=Couplings(),
    settings::SolverSettings=SolverSettings(k=3);
    output::AbstractString=joinpath(PACKAGE_ROOT, "output", "density"),
    force::Bool=false,
)
    output = ensure_output(output)
    write_run_metadata(output; command="density")
    path = joinpath(output, "density.csv")
    completed = force ? Set{String}() : completed_job_ids(path)
    model = build_model(nm1=nm1)
    cache = prepare_spectrum(model, couplings, settings)
    sector_lookup = Dict(sector.key => sector for sector in cache.sectors)
    for mu in Float64.(collect(mus))
        job_id = stable_id("density-v1", nm1, mu, coupling_vector(couplings))
        job_id in completed && continue
        try
            ground = solve_ground_state(cache, mu)
            sector = sector_lookup[ground.sector]
            # Nf=<ground|number_f|ground>；N0 不再另造算符，直接由总电荷推出。
            nf = real(dot(ground.vector, hermitian_opmat(sector.number_f) * ground.vector))
            n0 = (3nm1 - nf) / 3
            append_csv(path, (
                job_id=job_id, status="ok", timestamp=string(now()), nm1=nm1, mu=mu,
                energy=ground.energy, Nf=nf, nf_per_orbital=nf/nm1,
                N0=n0, n0_per_orbital=n0/model.nm0, z=ground.sector.z,
                r=ground.sector.r, error="", coupling_namedtuple(couplings)...,
            ))
        catch err
            append_csv(path, (
                job_id=job_id, status="error", timestamp=string(now()), nm1=nm1, mu=mu,
                energy=NaN, Nf=NaN, nf_per_orbital=NaN, N0=NaN,
                n0_per_orbital=NaN, z=0, r=0, error=sanitize_error(err),
                coupling_namedtuple(couplings)...,
            ))
        end
    end
    data = latest_rows(CSV.read(path, DataFrame))
    _plot_density_curve(data, joinpath(output, "density.png"))
    return data
end

"""
固定除 μ 外的参数，只在给定 μ 网格上选取 tower score 最小点。

不调用优化器。输出所有网格点 `critical_scan.csv`、最佳点
`critical_point.csv` 和所选 tower score 的残差表。
"""
function run_critical_search(
    nm1::Int,
    mus,
    couplings::Couplings=Couplings(),
    settings::SolverSettings=SolverSettings();
    definition=:critical5,
    terms=nothing,
    metric=:q,
    output::AbstractString=joinpath(PACKAGE_ROOT, "output", "critical"),
)
    output = ensure_output(output)
    write_run_metadata(output; command="critical")
    definition = normalize_score_definition(definition)
    selected_terms = resolve_score_terms(definition, terms)
    metric = normalize_score_metric(metric, definition)
    definition_tag = score_tag(definition, metric, selected_terms)
    terms_tag = join(String.(selected_terms), ",")
    model = build_model(nm1=nm1)
    cache = prepare_spectrum(model, couplings, settings)
    result = scan_mu(
        cache, mus; definition=definition, terms=selected_terms, metric=metric,
    )
    scan_rows = [(
        score_definition=definition_tag, score_terms=terms_tag, index=index, mu=mu,
        objective=score.objective, q=score.q, cost=score.cost, factor=score.factor,
        delta_s=score.delta_s, delta_o=score.delta_o,
        valid=score.valid, reason=score.reason,
    ) for (index, (mu, score)) in enumerate(zip(result.mus, result.scores))]
    atomic_csv(joinpath(output, "critical_scan.csv"), DataFrame(scan_rows))
    score = result.score
    summary = DataFrame([(
        score_definition=definition_tag, score_terms=terms_tag, timestamp=string(now()),
        nm1=nm1, mu=result.mu, objective=score.objective, q=score.q, cost=score.cost,
        factor=score.factor, delta_s=score.delta_s, delta_o=score.delta_o,
        valid=score.valid, reason=score.reason, completed=result.completed,
        at_boundary=result.at_boundary, evaluations=result.evaluations,
        invalid_evaluations=result.invalid, coupling_namedtuple(couplings)...,
    )])
    atomic_csv(joinpath(output, "critical_point.csv"), summary)
    score.valid && atomic_csv(joinpath(output, "tower_residuals.csv"), score_dataframe(score))
    return result
end

mutable struct _LinearSector
    key::SectorKey
    basis::Any
    fixed::SparseMatrixCSC{Float64,Int64}
    derivatives::Dict{Symbol,SparseMatrixCSC{Float64,Int64}}
    l2::Any
    c2::Any
    warm::Vector{Float64}
end

# 多参数优化的线性缓存：H = H_fixed + Σ pᵢ (∂H/∂pᵢ)。
# 这样每次 objective evaluation 只线性组合稀疏矩阵，不重建 FuzzifiED Terms。
struct _LinearFamily
    model::ModelParameters
    base::Couplings
    free::Vector{Symbol}
    settings::SolverSettings
    sectors::Vector{_LinearSector}
    u0_over_uf::Union{Nothing,Float64}
end

function _prepare_linear_family(model, base, free, settings; u0_over_uf=nothing)
    all(name -> name in HAMILTONIAN_FIELDS, free) || throw(ArgumentError("Unknown free Hamiltonian parameter"))
    u0_over_uf !== nothing && :U0 in free && throw(ArgumentError(
        "U0 cannot be free when u0_over_uf ties U0 to Uf",
    ))
    effective_base = u0_over_uf === nothing ? base :
                     with_coupling(base, :U0, Float64(u0_over_uf) * base.Uf)
    fixed_couplings = foldl(
        (c, name) -> with_coupling(c, name, 0.0), free; init=effective_base,
    )
    if u0_over_uf !== nothing && :Uf in free
        fixed_couplings = with_coupling(fixed_couplings, :U0, 0.0)
    end
    fixed_terms = hamiltonian_terms(model, fixed_couplings; include_mu=true)
    sectors = _LinearSector[]
    for z in (1, -1), r in (1, -1)
        key = SectorKey(z, r)
        basis = Basis(model.cfs[0], [z, r], model.qnf)
        basis.dim == 0 && continue
        fixed = lower_sparse(float_opmat(Operator(basis, fixed_terms)))
        # Hamiltonian 对每个线性耦合参数的“导数”就是对应分量矩阵。
        derivatives = Dict{Symbol,SparseMatrixCSC{Float64,Int64}}()
        for name in free
            derivative = lower_sparse(float_opmat(Operator(basis, getproperty(model.components, name))))
            if name == :Uf && u0_over_uf !== nothing
                derivative += Float64(u0_over_uf) * lower_sparse(
                    float_opmat(Operator(basis, model.components.U0)),
                )
            end
            derivatives[name] = derivative
        end
        push!(sectors, _LinearSector(
            key, basis, fixed, derivatives, float_opmat(Operator(basis, model.l2)),
            float_opmat(Operator(basis, model.c2)), Float64[],
        ))
    end
    return _LinearFamily(
        model, effective_base, collect(free), settings, sectors,
        u0_over_uf === nothing ? nothing : Float64(u0_over_uf),
    )
end

function _solve_linear(family::_LinearFamily, values::AbstractVector{<:Real})
    length(values) == length(family.free) || throw(DimensionMismatch("parameter count mismatch"))
    temporary = SectorCache[]
    for source in family.sectors
        matrix = copy(source.fixed)
        for (name, value) in zip(family.free, values)
            matrix += Float64(value) * source.derivatives[name]
        end
        zero_number = spzeros(Float64, size(matrix, 1), size(matrix, 2))
        push!(temporary, SectorCache(
            source.key, source.basis, matrix, zero_number, source.l2, source.c2,
            copy(source.warm),
        ))
    end
    cache = ModelCache(family.model, family.base, family.settings, temporary)
    states = solve_spectrum(cache, 0.0)
    for (source, updated) in zip(family.sectors, cache.sectors)
        source.warm = copy(updated.warm)
    end
    return states
end

function _couplings_from_values(family::_LinearFamily, values)
    couplings = foldl(
        (c, item) -> with_coupling(c, item[1], item[2]),
        zip(family.free, values); init=family.base,
    )
    if family.u0_over_uf !== nothing
        couplings = with_coupling(couplings, :U0, family.u0_over_uf * couplings.Uf)
    end
    return couplings
end

"""
连续优化单个 μ：保留旧 FSS1.jl 的有界 Brent 逻辑。

FSS 的 `optimize` 方法和单参数优化共用这个底层函数。
"""
function optimize_mu_with_score(
    cache::ModelCache;
    mu_min::Real,
    mu_max::Real,
    definition=:fss7,
    terms=nothing,
    metric=nothing,
    abs_tol::Real=1.0e-4,
    max_iterations::Int=60,
    penalty::Real=1.0e3,
)
    mu_min < mu_max || throw(ArgumentError("mu_min must be smaller than mu_max"))
    definition = normalize_score_definition(definition)
    selected_terms = resolve_score_terms(definition, terms)
    metric = normalize_score_metric(metric, definition)
    evaluations = Dict{Float64,CFTScore}()
    function objective(mu)
        key = round(Float64(mu); digits=14)
        if !haskey(evaluations, key)
            evaluations[key] = try
                cft_score(
                    solve_spectrum(cache, key); settings=cache.settings,
                    definition=definition, terms=selected_terms, metric=metric,
                )
            catch err
                _invalid_score(
                    definition, metric, sanitize_error(err); terms=selected_terms,
                )
            end
        end
        score = evaluations[key]
        @info "mu optimization" mu=key valid=score.valid objective=score.objective
        return score.valid && isfinite(score.objective) ? score.objective : Float64(penalty)
    end
    result = optimize(
        objective, Float64(mu_min), Float64(mu_max), Brent();
        abs_tol=Float64(abs_tol), rel_tol=1.0e-5, iterations=max_iterations,
        show_trace=false,
    )
    mu = Float64(Optim.minimizer(result))
    objective(mu)
    score = evaluations[round(mu; digits=14)]
    tolerance = max(Float64(abs_tol), 1.0e-8)
    return (
        mu=mu, score=score, completed=Optim.converged(result),
        at_boundary=abs(mu-mu_min) <= tolerance || abs(mu-mu_max) <= tolerance,
        evaluations=length(evaluations),
        invalid=count(value -> !value.valid, values(evaluations)), result=result,
    )
end

"""
按旧 `optimization.jl` 的方式同时优化若干 Hamiltonian 参数。

使用 NelderMead，越界点返回 penalty；默认使用旧八项 cost，并保留
`U0 = u0_over_uf * Uf` 约束。线性矩阵缓存只改善速度，不改变这个目标函数。
"""
function run_parameter_optimization(
    nm1::Int,
    base::Couplings,
    free::Vector{Symbol},
    bounds::Dict{Symbol,Tuple{Float64,Float64}},
    settings::SolverSettings=SolverSettings();
    max_iterations::Int=200,
    algorithm=:auto,
    abs_tol::Real=1.0e-4,
    definition=:optimization8,
    terms=nothing,
    metric=:cost,
    u0_over_uf=9.0,
    penalty::Real=1.0e6,
    output::AbstractString=joinpath(PACKAGE_ROOT, "output", "optimize"),
)
    output = ensure_output(output)
    write_run_metadata(output; command="optimize")
    definition = normalize_score_definition(definition)
    selected_terms = resolve_score_terms(definition, terms)
    metric = normalize_score_metric(metric, definition)
    definition_tag = score_tag(definition, metric, selected_terms)
    terms_tag = join(String.(selected_terms), ",")
    isempty(free) && throw(ArgumentError("optimization.free must contain at least one parameter"))
    length(unique(free)) == length(free) || throw(ArgumentError(
        "optimization.free must not contain duplicate parameters",
    ))
    all(haskey(bounds, name) for name in free) || throw(ArgumentError("Every free parameter needs a bound"))
    algorithm = Symbol(lowercase(replace(String(algorithm), '-' => '_')))
    algorithm == :neldermead && (algorithm = :nelder_mead)
    algorithm in (:auto, :brent, :nelder_mead) || throw(ArgumentError(
        "optimization.algorithm must be auto, brent, or nelder_mead",
    ))
    algorithm = algorithm == :auto ? (length(free) == 1 ? :brent : :nelder_mead) : algorithm
    algorithm == :brent && length(free) != 1 && throw(ArgumentError(
        "Brent can optimize exactly one free parameter; use auto or nelder_mead",
    ))
    free_tag = join(String.(free), ",")
    model = build_model(nm1=nm1)
    family = _prepare_linear_family(
        model, base, free, settings; u0_over_uf=u0_over_uf,
    )
    initial = Float64[getfield(base, name) for name in free]
    lower = Float64[bounds[name][1] for name in free]
    upper = Float64[bounds[name][2] for name in free]
    all((lower .<= initial) .& (initial .<= upper)) ||
        throw(ArgumentError("Initial parameters must lie inside bounds"))
    optimization_signature = stable_id(
        "optimization-v3", nm1, definition_tag, free, coupling_vector(family.base),
        lower, upper, settings, algorithm, u0_over_uf, penalty,
    )
    trace_path = joinpath(output, "evaluations.csv")
    evaluation = Ref(0)
    if isfile(trace_path)
        previous = CSV.read(trace_path, DataFrame)
        "score_definition" in names(previous) || error(
            "Existing optimization trace has no score label; use a new optimization case directory",
        )
        all(String(value) == definition_tag for value in previous.score_definition) ||
            error("Existing optimization trace uses a different score; use a new optimization case directory")
        "optimization_signature" in names(previous) || error(
            "Existing optimization trace has no full configuration signature; use a new optimization case directory",
        )
        all(String(value) == optimization_signature for value in previous.optimization_signature) ||
            error("Existing optimization trace uses different nm1/k/values/bounds; use a new optimization case directory")
        "free_parameters" in names(previous) || error(
            "Existing optimization trace has no exact free-parameter label; use a new optimization case directory",
        )
        all(String(value) == free_tag for value in previous.free_parameters) ||
            error("Existing optimization trace uses different free parameters; use a new optimization case directory")
        "algorithm" in names(previous) &&
            !all(String(value) == String(algorithm) for value in previous.algorithm) &&
            error("Existing optimization trace uses a different algorithm; use a new optimization case directory")
        valid_indices = [
            i for i in 1:nrow(previous)
            if Bool(previous.valid[i]) && isfinite(previous.objective[i])
        ]
        if !isempty(valid_indices)
            best_index = valid_indices[argmin(previous.objective[valid_indices])]
            resumed = Float64[previous[best_index, name] for name in free]
            if all((lower .<= resumed) .& (resumed .<= upper))
                initial = resumed
                @info "resuming optimization from best checkpoint" objective=previous.objective[best_index] parameters=initial
            end
        end
        evaluation[] = nrow(previous)
    end

    # objective 的一次调用通常意味着四个 sector 各做一次低能稀疏对角化。
    function objective(values)
        evaluation[] += 1
        inside = all((lower .<= values) .& (values .<= upper))
        score = if !inside
            _invalid_score(
                definition, metric, "parameters outside configured bounds";
                terms=selected_terms,
            )
        else
            try
                cft_score(
                    _solve_linear(family, values); settings=settings,
                    definition=definition, terms=selected_terms, metric=metric,
                )
            catch err
                @error "parameter evaluation failed" values exception=(err, catch_backtrace())
                _invalid_score(
                    definition, metric, sanitize_error(err); terms=selected_terms,
                )
            end
        end
        row_parameters = (; (name => Float64(value) for (name, value) in zip(free, values))...)
        effective = _couplings_from_values(family, values)
        effective_parameters = (;
            (Symbol("effective_", name) => getfield(effective, name)
             for name in HAMILTONIAN_FIELDS)...,
        )
        append_csv(trace_path, (
            score_definition=definition_tag, score_terms=terms_tag,
            optimization_signature=optimization_signature,
            free_parameters=free_tag, algorithm=String(algorithm),
            evaluation=evaluation[],
            timestamp=string(now()), objective=score.objective, q=score.q,
            cost=score.cost, factor=score.factor,
            valid=score.valid, reason=score.reason, row_parameters...,
            effective_parameters...,
        ))
        @info "optimization evaluation" evaluation=evaluation[] objective=score.objective parameters=row_parameters
        return score.valid && isfinite(score.objective) ? score.objective : Float64(penalty)
    end

    result = if algorithm == :brent
        scalar_objective(value) = objective([Float64(value)])
        optimize(
            scalar_objective, lower[1], upper[1], Brent();
            abs_tol=Float64(abs_tol), rel_tol=1.0e-5,
            iterations=max_iterations, show_trace=false,
        )
    else
        optimize(
            objective, initial, NelderMead(),
            Optim.Options(
                iterations=max_iterations, g_tol=1.0e-7,
                show_trace=false,
            ),
        )
    end
    values = algorithm == :brent ? [Float64(Optim.minimizer(result))] :
             Float64.(Optim.minimizer(result))
    best = _couplings_from_values(family, values)
    score = cft_score(
        _solve_linear(family, values); settings=settings,
        definition=definition, terms=selected_terms, metric=metric,
    )
    atomic_csv(joinpath(output, "best.csv"), DataFrame([(
        score_definition=definition_tag, score_terms=terms_tag,
        optimization_signature=optimization_signature,
        free_parameters=free_tag, nm1=nm1, objective=score.objective,
        q=score.q, cost=score.cost, factor=score.factor, delta_s=score.delta_s,
        delta_o=score.delta_o, algorithm=String(algorithm),
        converged=Optim.converged(result),
        evaluations=evaluation[], coupling_namedtuple(best)...,
    )]))
    return (couplings=best, score=score, result=result, evaluations=evaluation[])
end

"""
完整 finite-size scaling 数据生成流程。

对每个 `nm1` 和外层 `scan_value`，`grid` 在固定 μ 网格上选最小点；
`optimize` 沿用旧 FSS1.jl 的 Brent 连续寻找 μc。两种结果分文件保存。
"""
function run_fss_scan(
    base::Couplings,
    fss::FSSSettings,
    settings::SolverSettings=SolverSettings();
    output::AbstractString=joinpath(PACKAGE_ROOT, "output", "fss"),
    force::Bool=false,
)
    fss.scan_parameter in HAMILTONIAN_FIELDS || throw(ArgumentError("Unknown scan parameter $(fss.scan_parameter)"))
    fss.scan_parameter == :mu && throw(ArgumentError("FSS scan_parameter cannot be mu; mu already has its own grid"))
    fss.mu_count > 0 || throw(ArgumentError("FSS mu_count must be positive"))
    methods = unique(fss.methods)
    isempty(methods) && throw(ArgumentError("FSS methods must not be empty"))
    all(method -> method in (:grid, :optimize), methods) ||
        throw(ArgumentError("FSS methods must contain only grid and/or optimize"))
    definition = normalize_score_definition(fss.score_definition)
    selected_terms = resolve_score_terms(
        definition, isempty(fss.score_terms) ? nothing : fss.score_terms,
    )
    metric = normalize_score_metric(fss.score_metric, definition)
    definition_tag = score_tag(definition, metric, selected_terms)
    terms_tag = join(String.(selected_terms), ",")
    output = ensure_output(output)
    write_run_metadata(output; command="fss")
    paths = Dict(
        method => joinpath(output, "fss_$(method)_results.csv") for method in methods
    )
    completed = Dict{Symbol,Set{String}}()
    for method in methods
        path = paths[method]
        if isfile(path)
            previous = CSV.read(path, DataFrame)
            "score_definition" in names(previous) || error(
                "Existing $method FSS data has no score label; use a new FSS case directory",
            )
            all(String(value) == definition_tag for value in previous.score_definition) ||
                error("Existing $method FSS data uses a different score; use a new FSS case directory")
        end
        completed[method] = force ? Set{String}() : completed_job_ids(path)
    end

    for nm1 in fss.nm_values
        model = build_model(nm1=nm1)
        cache = nothing
        for scan_value in fss.scan_values
            couplings = with_coupling(base, fss.scan_parameter, scan_value)
            job_ids = Dict(method => stable_id(
                "fss-$method-v1", definition_tag, nm1, fss.scan_parameter, scan_value,
                coupling_vector(couplings), settings.k, fss.mu_min, fss.mu_max,
                method == :grid ? fss.mu_count : fss.optimize_abs_tol,
            ) for method in methods)
            pending = [method for method in methods if !(job_ids[method] in completed[method])]
            isempty(pending) && continue
            if cache === nothing
                cache = prepare_spectrum(model, couplings, settings)
            else
                retune_spectrum!(cache, couplings)
            end

            for method in pending
                job_id = job_ids[method]
                path = paths[method]
                @info "FSS job" method nm1 parameter=fss.scan_parameter scan_value job_id
                try
                    result = if method == :grid
                        mus = range(fss.mu_min, fss.mu_max; length=fss.mu_count)
                        scan_mu(
                            cache, mus; definition=definition,
                            terms=selected_terms, metric=metric,
                        )
                    else
                        optimize_mu_with_score(
                            cache; mu_min=fss.mu_min, mu_max=fss.mu_max,
                            definition=definition, terms=selected_terms, metric=metric,
                            abs_tol=fss.optimize_abs_tol,
                            max_iterations=fss.optimize_max_iterations,
                        )
                    end
                    score = result.score
                    append_csv(path, (
                        score_definition=definition_tag, score_terms=terms_tag,
                        method=String(method),
                        job_id=job_id, status="ok", timestamp=string(now()), nm1=nm1,
                        x=nm1^(-0.5), scan_parameter=String(fss.scan_parameter),
                        scan_value=scan_value, muc=result.mu, objective=score.objective,
                        q=score.q, cost=score.cost, factor=score.factor,
                        delta_s=score.delta_s, delta_o=score.delta_o,
                        score_valid=score.valid, score_reason=score.reason,
                        completed=result.completed, at_boundary=result.at_boundary,
                        evaluations=result.evaluations, invalid=result.invalid,
                        error="", coupling_namedtuple(couplings)...,
                    ))
                catch err
                    @error "FSS job failed" method nm1 scan_value exception=(err, catch_backtrace())
                    append_csv(path, (
                        score_definition=definition_tag, score_terms=terms_tag,
                        method=String(method),
                        job_id=job_id, status="error", timestamp=string(now()), nm1=nm1,
                        x=nm1^(-0.5), scan_parameter=String(fss.scan_parameter),
                        scan_value=scan_value, muc=NaN, objective=Inf, q=Inf,
                        cost=Inf, factor=NaN, delta_s=NaN, delta_o=NaN,
                        score_valid=false, score_reason="", completed=false,
                        at_boundary=false, evaluations=0, invalid=0,
                        error=sanitize_error(err), coupling_namedtuple(couplings)...,
                    ))
                end
            end
        end
    end
    return Dict(
        method => latest_rows(CSV.read(paths[method], DataFrame)) for method in methods
    )
end

function _valid_fss(data::DataFrame, y::Symbol)
    String(y) in names(data) || throw(ArgumentError("Column $y is absent"))
    "score_definition" in names(data) || error("FSS data has no score definition")
    length(unique(String.(data.score_definition))) == 1 ||
        error("FSS file mixes different score definitions")
    return filter(row -> row.status == "ok" && row.score_valid && isfinite(row[y]), latest_rows(data))
end

"""读取一份 `fss_<method>_results.csv`，把指定观测量随 `x=nm1^(-1/2)` 画出。"""
function plot_fss(
    source::AbstractString;
    y::Symbol=:delta_s,
    output::AbstractString=joinpath(dirname(source), "$(y)_fss.png"),
)
    data = _valid_fss(CSV.read(source, DataFrame), y)
    nrow(data) > 0 || error("No valid FSS rows to plot")
    fig = Figure(size=(600, 600))
    scan_name = first(data.scan_parameter)
    axis = Axis(
        fig[1, 1]; xlabel="nm1^(-1/2)", ylabel=String(y),
        title="$(String(y)) vs nm1^(-1/2), $scan_name scan", aspect=1,
    )
    # 旧图只画 nm1=4:7，因此上限固定为 0.52；配置含 nm1=3 时也要显示 x=1/√3。
    xlims!(axis, 0.0, max(0.52, 1.05 * maximum(data.x)))
    y == :delta_s && ylims!(axis, 1.0, 2.0)
    for value in sort(unique(data.scan_value))
        sub = sort(filter(row -> row.scan_value == value, data), :x)
        scatter!(axis, sub.x, sub[!, y]; markersize=10, label=string(value))
    end
    axislegend(axis; title=scan_name, position=:lt, nbanks=2, framevisible=true)
    save(output, fig)
    return fig
end

"""
对所有 scan_value 联合拟合 `Δ(x,p)=Δ∞+a(p)x^ω`。

`Δ∞` 和修正指数 `ω` 由所有曲线共享，每条参数曲线有自己的振幅 `a(p)`。
至少需要三个系统大小，返回外推值并保存拟合表和图。
"""
function fit_fss(
    source::AbstractString;
    y::Symbol=:delta_s,
    output::AbstractString=dirname(source),
    label::AbstractString="",
)
    data = _valid_fss(CSV.read(source, DataFrame), y)
    values = sort(unique(Float64.(data.scan_value)))
    length(unique(data.nm1)) >= 3 || error("Joint FSS fit requires at least three system sizes")
    count = nrow(data)
    parameter_count = 2 + length(values)
    count > parameter_count || error("Need more data points than fit parameters ($count <= $parameter_count)")
    value_index = Dict(value => i for (i, value) in enumerate(values))
    x = Float64.(data.x)
    observed = Float64.(data[!, y])
    groups = [value_index[Float64(value)] for value in data.scan_value]

    # 固定非线性参数 ω 后，Δ∞ 与所有 amplitude 都可一次线性最小二乘求出。
    function linear_fit(logomega)
        omega = exp(logomega)
        design = zeros(Float64, count, parameter_count - 1)
        design[:, 1] .= 1.0
        for i in eachindex(x)
            design[i, 1 + groups[i]] = x[i]^omega
        end
        coefficients = design \ observed
        residual = design * coefficients - observed
        return coefficients, sum(abs2, residual)
    end
    objective(logomega) = last(linear_fit(logomega))
    log_grid = collect(range(log(0.05), log(10.0); length=81))
    costs = objective.(log_grid)
    best_index = argmin(costs)
    # fss-fit 也只在明确的固定网格上选最小值，不调用 Optim。
    logomega = log_grid[best_index]
    coefficients, rss = linear_fit(logomega)
    delta_inf = coefficients[1]
    omega = exp(logomega)
    amplitudes = coefficients[2:end]
    fit_table = DataFrame(
        scan_value=values, amplitude=amplitudes, delta_inf=fill(delta_inf, length(values)),
        omega=fill(omega, length(values)), rss=fill(rss, length(values)),
        npoints=fill(count, length(values)), grid_points=fill(length(log_grid), length(values)),
    )
    mkpath(output)
    suffix = isempty(label) ? "" : "_$(label)"
    atomic_csv(joinpath(output, "$(y)$(suffix)_fit.csv"), fit_table)

    fig = Figure(size=(680, 560))
    axis = Axis(fig[1, 1], xlabel="nm1^(-1/2)", ylabel=String(y), title="joint FSS fit")
    xline = range(0.0, maximum(x); length=200)
    for (index, value) in enumerate(values)
        sub = filter(row -> row.scan_value == value, data)
        scatter!(axis, sub.x, sub[!, y]; label=string(value))
        lines!(axis, xline, delta_inf .+ amplitudes[index] .* xline.^omega)
    end
    scatter!(axis, [0.0], [delta_inf]; marker=:star5, markersize=16, color=:black)
    axislegend(axis; title=first(data.scan_parameter))
    save(joinpath(output, "$(y)$(suffix)_fit.png"), fig)
    return (
        delta_inf=delta_inf, omega=omega, amplitudes=Dict(zip(values, amplitudes)),
        omega_grid=exp.(log_grid), grid_costs=costs,
    )
end

function _legacy_scaling_factor(states::Vector{SpectrumState}, quantum_tol::Real)
    sector(l2, c2) = _states_in_sector(states, l2, c2; quantum_tol=quantum_tol)

    s00 = sector(0, 0)
    s20 = sector(2, 0)
    s60 = sector(6, 0)
    a23 = sector(2, 3)
    a63 = sector(6, 3)
    length(s00) >= 2 || error("Scaling factor needs two raw (L2,C2)=(0,0) states")
    length(s20) >= 1 || error("Scaling factor needs one raw (L2,C2)=(2,0) state")
    length(s60) >= 1 || error("Scaling factor needs one raw (L2,C2)=(6,0) state")
    length(a23) >= 3 || error("Scaling factor needs three raw (L2,C2)=(2,3) states")
    length(a63) >= 1 || error("Scaling factor needs one raw (L2,C2)=(6,3) state")

    ground = first(states).energy
    gaps = Float64[
        s20[1].energy - s00[2].energy,
        a23[1].energy - ground,
        a23[3].energy - ground,
        a63[1].energy - ground,
        s60[1].energy - ground,
    ]
    targets = Float64[1, 2, 3, 3, 3]
    factor = dot(gaps, targets) / dot(targets, targets)
    isfinite(factor) && factor > 0 || error("The legacy scaling factor is not positive and finite")
    return factor
end

"""
按照旧 `scaling_dimension.jl` 的规则画 scaling-dimension spectrum。

未显式给 `factor` 时，使用旧代码的五条能隙关系。量子数先在容差内认成
整数，再按 `L²≤20, C₂≤8` 一类的整数边界筛选，避免浮点误差漏掉边界态。
"""
function plot_scaling_dimensions(
    nm1::Int,
    couplings::Couplings=Couplings(),
    settings::SolverSettings=SolverSettings();
    factor::Union{Nothing,Real}=nothing,
    l2_max::Real=20.0,
    c2_max::Real=8.0,
    output::AbstractString=joinpath(PACKAGE_ROOT, "output", "scaling"),
)
    output = ensure_output(output)
    model = build_model(nm1=nm1)
    cache = prepare_spectrum(model, couplings, settings)
    states = solve_spectrum(cache, couplings.mu)
    scale = isnothing(factor) ? _legacy_scaling_factor(states, settings.quantum_tol) : Float64(factor)
    isfinite(scale) && scale > 0 || error("A valid positive factor is required")
    ground = minimum(state.energy for state in states)
    rows = NamedTuple[]
    for state in states
        labels = _state_quantum_labels(state, settings.quantum_tol)
        if isnothing(labels)
            continue
        end
        l2, c2 = labels
        l2 <= l2_max && c2 <= c2_max || continue
        ell = round(Int, (sqrt(max(0.0, 1 + 4l2)) - 1) / 2)
        push!(rows, (
            energy=state.energy, delta=(state.energy-ground)/scale, ell=ell,
            l2=l2, c2=c2, z=state.sector.z, r=state.sector.r,
        ))
    end
    data = DataFrame(rows)
    atomic_csv(joinpath(output, "scaling_nm$(nm1).csv"), data)
    colors = Dict(0 => :red, 3 => :blue, 6 => :green, 8 => :purple)
    labels = ((0, "C₂=0"), (3, "C₂=3"), (6, "C₂=6"), (8, "C₂=8"))
    fig = Figure(size=(600, 700))
    axis = Axis(
        fig[1, 1]; xlabel="l", ylabel="Δ",
        title="Nm=$nm1: scaling dimension vs l", aspect=DataAspect(),
    )
    axis.xticks = 0:1:4
    axis.yticks = 0:1:5
    xlims!(axis, -0.3, 4.3)
    ylims!(axis, -0.3, 5.3)
    for row in eachrow(data)
        color = get(colors, row.c2, :gray)
        lines!(axis, [row.ell-0.15, row.ell+0.15], [row.delta, row.delta]; color=color, linewidth=2)
        scatter!(axis, [row.ell], [row.delta]; color=color, markersize=10)
    end
    for (c2, label) in labels
        lines!(axis, [NaN, NaN], [NaN, NaN]; color=colors[c2], linewidth=2, label=label)
    end
    axislegend(axis; position=:rb, framevisible=false)
    save(joinpath(output, "scaling_nm$(nm1).png"), fig)
    return (data=data, factor=scale, figure=fig)
end
