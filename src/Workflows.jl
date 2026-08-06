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
对一个系统大小和一列 μ 求低能谱，并保存每点的谱与 CFT 摘要。

输出 `summary.csv` 以及 `spectra/<job_id>.csv`；`force=false` 时跳过已成功
完成的 job。只有确实需要本征向量时才设 `keep_vectors=true`，否则文件很大。
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
            score = cft_score(states; settings=settings)
            table = spectrum_dataframe(states; mu=mu, nm1=nm1)
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

function _level_energy(catalog, key, rank)
    levels = get(catalog, key, PhysicalLevel[])
    return length(levels) >= rank ? levels[rank].energy : NaN
end

"""
扫描多个系统大小和 μ，提取第二个 `(L²,C₂)=(0,0)` 能级相对基态的能隙。

同时保存 `gap*sqrt(nm1)`，并自动画不同大小随 μ 的交叉图。
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
    path = joinpath(output, "singlet_gap.csv")
    completed = force ? Set{String}() : completed_job_ids(path)
    for nm1 in Int.(nm_values)
        model = build_model(nm1=nm1)
        cache = prepare_spectrum(model, couplings, settings)
        for mu in Float64.(collect(mus))
            job_id = stable_id("gap-v1", nm1, mu, coupling_vector(couplings), settings.k)
            job_id in completed && continue
            try
                states = solve_spectrum(cache, mu)
                catalog, _ = level_catalog(
                    states; quantum_tol=settings.quantum_tol,
                    degeneracy_tol=settings.degeneracy_tol,
                )
                ground = minimum(state.energy for state in states)
                singlet = _level_energy(catalog, (0, 0), 2)
                gap = singlet - ground
                append_csv(path, (
                    job_id=job_id, status="ok", timestamp=string(now()), nm1=nm1,
                    mu=mu, gap=gap, scaled_gap=gap*sqrt(nm1), error="",
                    coupling_namedtuple(couplings)...,
                ))
            catch err
                append_csv(path, (
                    job_id=job_id, status="error", timestamp=string(now()), nm1=nm1,
                    mu=mu, gap=NaN, scaled_gap=NaN, error=sanitize_error(err),
                    coupling_namedtuple(couplings)...,
                ))
            end
        end
    end
    data = latest_rows(CSV.read(path, DataFrame))
    valid = filter(row -> row.status == "ok" && isfinite(row.scaled_gap), data)
    if nrow(valid) > 0
        fig = Figure(size=(720, 520))
        axis = Axis(fig[1, 1], xlabel="mu", ylabel="singlet gap * sqrt(nm1)")
        for nm1 in sort(unique(valid.nm1))
            sub = sort(filter(row -> row.nm1 == nm1, valid), :mu)
            lines!(axis, sub.mu, sub.scaled_gap; label="nm1=$nm1")
            scatter!(axis, sub.mu, sub.scaled_gap)
        end
        axislegend(axis)
        save(joinpath(output, "singlet_gap.png"), fig)
    end
    return data
end

"""
扫描 μ 并计算基态中的 charge-1/charge-3 平均粒子数及每轨道密度。

这里必须保留基态向量；由总电荷约束 `Nf+3N0=3nm1` 推出 N0。
"""
function run_density_scan(
    nm1::Int,
    mus,
    couplings::Couplings=Couplings(),
    settings::SolverSettings=SolverSettings(k=4);
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
            states = solve_spectrum(cache, mu; keep_vectors=true)
            ground = first(states)
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
    valid = filter(row -> row.status == "ok", data)
    if nrow(valid) > 0
        sort!(valid, :mu)
        fig = Figure(size=(680, 520))
        axis = Axis(fig[1, 1], xlabel="mu", ylabel="<Nf>/nm1", title="charge-1 fermion density")
        lines!(axis, valid.mu, valid.nf_per_orbital)
        scatter!(axis, valid.mu, valid.nf_per_orbital)
        save(joinpath(output, "density.png"), fig)
    end
    return data
end

"""
固定除 μ 外的参数，在一个系统大小上寻找 tower score 最小的临界 μ。

输出最佳点摘要 `critical_point.csv` 和七个 tower 关系的残差表。
"""
function run_critical_search(
    nm1::Int,
    couplings::Couplings=Couplings(),
    settings::SolverSettings=SolverSettings();
    mu_min::Real=0.0,
    mu_max::Real=0.12,
    coarse_points::Int=9,
    output::AbstractString=joinpath(PACKAGE_ROOT, "output", "critical"),
)
    output = ensure_output(output)
    write_run_metadata(output; command="critical")
    model = build_model(nm1=nm1)
    cache = prepare_spectrum(model, couplings, settings)
    result = optimize_mu(
        cache; mu_min=mu_min, mu_max=mu_max, coarse_points=coarse_points,
    )
    score = result.score
    summary = DataFrame([(
        timestamp=string(now()), nm1=nm1, mu=result.mu, q=score.q,
        factor=score.factor, delta_s=score.delta_s, delta_o=score.delta_o,
        valid=score.valid, reason=score.reason, converged=result.converged,
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
end

function _prepare_linear_family(model, base, free, settings)
    all(name -> name in HAMILTONIAN_FIELDS, free) || throw(ArgumentError("Unknown free Hamiltonian parameter"))
    fixed_couplings = foldl((c, name) -> with_coupling(c, name, 0.0), free; init=base)
    fixed_terms = hamiltonian_terms(model, fixed_couplings; include_mu=true)
    sectors = _LinearSector[]
    for z in (1, -1), r in (1, -1)
        key = SectorKey(z, r)
        basis = Basis(model.cfs[0], [z, r], model.qnf)
        basis.dim == 0 && continue
        fixed = lower_sparse(float_opmat(Operator(basis, fixed_terms)))
        # Hamiltonian 对每个线性耦合参数的“导数”就是对应分量矩阵。
        derivatives = Dict(
            name => lower_sparse(float_opmat(Operator(basis, getproperty(model.components, name))))
            for name in free
        )
        push!(sectors, _LinearSector(
            key, basis, fixed, derivatives, float_opmat(Operator(basis, model.l2)),
            float_opmat(Operator(basis, model.c2)), Float64[],
        ))
    end
    return _LinearFamily(model, base, collect(free), settings, sectors)
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

"""
在给定边界内同时优化若干 Hamiltonian 参数，使 CFT tower score 最小。

使用 Fminbox(NelderMead)，每次评价立即追加到 `evaluations.csv`；若输出目录
已有兼容记录，会从历史最好点继续。最终参数写入 `best.csv`。
"""
function run_parameter_optimization(
    nm1::Int,
    base::Couplings,
    free::Vector{Symbol},
    bounds::Dict{Symbol,Tuple{Float64,Float64}},
    settings::SolverSettings=SolverSettings();
    max_iterations::Int=200,
    output::AbstractString=joinpath(PACKAGE_ROOT, "output", "optimization"),
)
    output = ensure_output(output)
    write_run_metadata(output; command="optimize")
    all(haskey(bounds, name) for name in free) || throw(ArgumentError("Every free parameter needs a bound"))
    model = build_model(nm1=nm1)
    family = _prepare_linear_family(model, base, free, settings)
    initial = Float64[getfield(base, name) for name in free]
    lower = Float64[bounds[name][1] for name in free]
    upper = Float64[bounds[name][2] for name in free]
    all((lower .<= initial) .& (initial .<= upper)) ||
        throw(ArgumentError("Initial parameters must lie inside bounds"))
    trace_path = joinpath(output, "evaluations.csv")
    evaluation = Ref(0)
    if isfile(trace_path)
        previous = CSV.read(trace_path, DataFrame)
        all(String(name) in names(previous) for name in free) ||
            error("Existing optimization trace uses different free parameters; choose a new output.run_name")
        valid_indices = [
            i for i in 1:nrow(previous)
            if Bool(previous.valid[i]) && isfinite(previous.q[i])
        ]
        if !isempty(valid_indices)
            best_index = valid_indices[argmin(previous.q[valid_indices])]
            resumed = Float64[previous[best_index, name] for name in free]
            if all((lower .<= resumed) .& (resumed .<= upper))
                initial = resumed
                @info "resuming optimization from best checkpoint" q=previous.q[best_index] parameters=initial
            end
        end
        evaluation[] = nrow(previous)
    end

    # objective 的一次调用通常意味着四个 sector 各做一次低能稀疏对角化。
    function objective(values)
        evaluation[] += 1
        score = try
            cft_score(_solve_linear(family, values); settings=settings)
        catch err
            @error "parameter evaluation failed" values exception=(err, catch_backtrace())
            CFTScore(reason=sanitize_error(err))
        end
        row_parameters = (; (name => Float64(value) for (name, value) in zip(free, values))...)
        append_csv(trace_path, (
            evaluation=evaluation[], timestamp=string(now()), q=score.q,
            valid=score.valid, reason=score.reason, row_parameters...,
        ))
        @info "optimization evaluation" evaluation=evaluation[] q=score.q parameters=row_parameters
        return score.valid && isfinite(score.q) ? score.q : 1.0e3
    end

    result = optimize(
        objective, lower, upper, initial, Fminbox(NelderMead()),
        Optim.Options(iterations=max_iterations, show_trace=false),
    )
    values = Float64.(Optim.minimizer(result))
    best = foldl(
        (c, item) -> with_coupling(c, item[1], item[2]),
        zip(free, values); init=base,
    )
    score = cft_score(_solve_linear(family, values); settings=settings)
    atomic_csv(joinpath(output, "best.csv"), DataFrame([(
        nm1=nm1, q=score.q, factor=score.factor, delta_s=score.delta_s,
        delta_o=score.delta_o, converged=Optim.converged(result),
        evaluations=evaluation[], coupling_namedtuple(best)...,
    )]))
    return (couplings=best, score=score, result=result, evaluations=evaluation[])
end

"""
完整 finite-size scaling 数据生成流程。

对每个 `nm1` 和外层参数 `scan_value`，先改变该 Hamiltonian 系数，再在内层
优化 μ。每个点保存最佳 μ、q、ΔS、ΔO 与边界/收敛诊断到 `fss_results.csv`。
"""
function run_fss_scan(
    base::Couplings,
    fss::FSSSettings,
    settings::SolverSettings=SolverSettings();
    output::AbstractString=joinpath(PACKAGE_ROOT, "output", "fss"),
    force::Bool=false,
)
    fss.scan_parameter in HAMILTONIAN_FIELDS || throw(ArgumentError("Unknown scan parameter $(fss.scan_parameter)"))
    fss.scan_parameter == :mu && throw(ArgumentError("FSS scan_parameter cannot be mu; mu is optimized separately"))
    output = ensure_output(output)
    write_run_metadata(output; command="fss")
    path = joinpath(output, "fss_results.csv")
    completed = force ? Set{String}() : completed_job_ids(path)

    for nm1 in fss.nm_values
        model = build_model(nm1=nm1)
        # 同一 nm1 下 basis/observables 始终复用；scan_value 改变时仅 retune H0。
        cache = nothing
        for scan_value in fss.scan_values
            couplings = with_coupling(base, fss.scan_parameter, scan_value)
            job_id = stable_id(
                "fss-v2", nm1, fss.scan_parameter, scan_value,
                coupling_vector(couplings), settings.k, fss.mu_min, fss.mu_max,
            )
            job_id in completed && continue
            @info "FSS job" nm1 parameter=fss.scan_parameter scan_value job_id
            try
                if cache === nothing
                    cache = prepare_spectrum(model, couplings, settings)
                else
                    retune_spectrum!(cache, couplings)
                end
                result = optimize_mu(
                    cache; mu_min=fss.mu_min, mu_max=fss.mu_max,
                    coarse_points=fss.coarse_points, abs_tol=fss.mu_abs_tol,
                    max_iterations=fss.max_iterations,
                )
                score = result.score
                append_csv(path, (
                    job_id=job_id, status="ok", timestamp=string(now()), nm1=nm1,
                    x=nm1^(-0.5), scan_parameter=String(fss.scan_parameter),
                    scan_value=scan_value, mu=result.mu, q=score.q,
                    factor=score.factor, delta_s=score.delta_s, delta_o=score.delta_o,
                    score_valid=score.valid, score_reason=score.reason,
                    converged=result.converged, at_boundary=result.at_boundary,
                    evaluations=result.evaluations, invalid=result.invalid,
                    error="", coupling_namedtuple(couplings)...,
                ))
            catch err
                @error "FSS job failed" nm1 scan_value exception=(err, catch_backtrace())
                append_csv(path, (
                    job_id=job_id, status="error", timestamp=string(now()), nm1=nm1,
                    x=nm1^(-0.5), scan_parameter=String(fss.scan_parameter),
                    scan_value=scan_value, mu=NaN, q=Inf, factor=NaN,
                    delta_s=NaN, delta_o=NaN, score_valid=false,
                    score_reason="", converged=false, at_boundary=false,
                    evaluations=0, invalid=0, error=sanitize_error(err),
                    coupling_namedtuple(couplings)...,
                ))
            end
        end
    end
    return latest_rows(CSV.read(path, DataFrame))
end

function _valid_fss(data::DataFrame, y::Symbol)
    String(y) in names(data) || throw(ArgumentError("Column $y is absent"))
    return filter(row -> row.status == "ok" && row.score_valid && isfinite(row[y]), latest_rows(data))
end

"""读取 `fss_results.csv`，把指定观测量随 `x=nm1^(-1/2)` 画出。"""
function plot_fss(
    source::AbstractString;
    y::Symbol=:delta_s,
    output::AbstractString=joinpath(dirname(source), "$(y)_fss.png"),
)
    data = _valid_fss(CSV.read(source, DataFrame), y)
    nrow(data) > 0 || error("No valid FSS rows to plot")
    fig = Figure(size=(680, 560))
    axis = Axis(fig[1, 1], xlabel="nm1^(-1/2)", ylabel=String(y), title="finite-size scaling")
    for value in sort(unique(data.scan_value))
        sub = sort(filter(row -> row.scan_value == value, data), :x)
        lines!(axis, sub.x, sub[!, y]; label=string(value))
        scatter!(axis, sub.x, sub[!, y])
    end
    axislegend(axis; title=first(data.scan_parameter))
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
    left = log_grid[max(1, best_index-1)]
    right = log_grid[min(length(log_grid), best_index+1)]
    result = optimize(objective, left, right, Brent(); abs_tol=1e-8, rel_tol=0.0)
    logomega = Float64(Optim.minimizer(result))
    coefficients, rss = linear_fit(logomega)
    delta_inf = coefficients[1]
    omega = exp(logomega)
    amplitudes = coefficients[2:end]
    fit_table = DataFrame(
        scan_value=values, amplitude=amplitudes, delta_inf=fill(delta_inf, length(values)),
        omega=fill(omega, length(values)), rss=fill(rss, length(values)),
        npoints=fill(count, length(values)), converged=fill(Optim.converged(result), length(values)),
    )
    mkpath(output)
    atomic_csv(joinpath(output, "$(y)_fit.csv"), fit_table)

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
    save(joinpath(output, "$(y)_fit.png"), fig)
    return (delta_inf=delta_inf, omega=omega, amplitudes=Dict(zip(values, amplitudes)), result=result)
end

"""
把一个系统大小的低能能隙转换成标度维数并按角动量画 tower 图。

若未显式给 `factor`，使用 `cft_score` 拟合出的 factor；输出一行一个态的
CSV 和 scaling-dimension 图。
"""
function plot_scaling_dimensions(
    nm1::Int,
    couplings::Couplings=Couplings(),
    settings::SolverSettings=SolverSettings();
    factor::Union{Nothing,Real}=nothing,
    l2_max::Real=20.0,
    c2_max::Real=12.0,
    output::AbstractString=joinpath(PACKAGE_ROOT, "output", "scaling"),
)
    output = ensure_output(output)
    model = build_model(nm1=nm1)
    cache = prepare_spectrum(model, couplings, settings)
    states = solve_spectrum(cache, couplings.mu)
    score = cft_score(states; settings=settings)
    scale = isnothing(factor) ? score.factor : Float64(factor)
    isfinite(scale) && scale > 0 || error("A valid positive factor is required")
    ground = minimum(state.energy for state in states)
    rows = NamedTuple[]
    for state in states
        state.l2 <= l2_max && state.c2 <= c2_max || continue
        ell = (sqrt(max(0.0, 1 + 4state.l2)) - 1) / 2
        push!(rows, (
            energy=state.energy, delta=(state.energy-ground)/scale, ell=ell,
            l2=state.l2, c2=state.c2, z=state.sector.z, r=state.sector.r,
        ))
    end
    data = DataFrame(rows)
    atomic_csv(joinpath(output, "scaling_nm$(nm1).csv"), data)
    fig = Figure(size=(680, 720))
    axis = Axis(fig[1, 1], xlabel="ell", ylabel="Delta", title="nm1=$nm1, factor=$(round(scale; digits=6))")
    palette = Makie.wong_colors()
    for (index, c2value) in enumerate(sort(unique(round.(Int, data.c2))))
        sub = filter(row -> round(Int, row.c2) == c2value, data)
        color = palette[mod1(index, length(palette))]
        scatter!(axis, sub.ell, sub.delta; color=color, label="C2=$c2value")
        for row in eachrow(sub)
            lines!(axis, [row.ell-0.12, row.ell+0.12], [row.delta, row.delta]; color=color)
        end
    end
    axislegend(axis)
    save(joinpath(output, "scaling_nm$(nm1).png"), fig)
    return (data=data, score=score, factor=scale, figure=fig)
end
