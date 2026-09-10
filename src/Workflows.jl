# 本文件是项目的“业务流程层”：把 Model/Spectrum/CFT/Storage 中的基础函数
# 组合成用户真正会运行的 spectrum、gap、density、critical、optimization、FSS。

"""
从已经分类的本征态生成精简物理能级表。

不同离散 `(Z,R)` sector 中同 `(L²,C₂)`、同能量的副本先由 `level_catalog`
合并；能量不同的能级始终保留。每个请求的 `(L²,C₂)` 组合成为一列，列内是
最低 `levels_per_block` 个 `rescaled_energy=(E-E₀)/factor`。列按 `C₂` 在外、
`L²` 在内排列，方便先横向比较 singlet，再比较 adjoint；不猜测算符身份。
"""
function spectrum_level_table(
    states::Vector{SpectrumState};
    l2_values::AbstractVector{<:Integer},
    c2_values::AbstractVector{<:Integer},
    levels_per_block::Integer=7,
    factor::Real,
    quantum_tol::Real=2.0e-3,
    degeneracy_tol::Real=2.0e-6,
)
    isempty(states) && throw(ArgumentError("Cannot build a spectrum table from no states"))
    levels_per_block > 0 || throw(ArgumentError("levels_per_block must be positive"))
    isfinite(factor) && factor > 0 || throw(ArgumentError("factor must be positive and finite"))

    requested_l2 = unique(Int.(l2_values))
    requested_c2 = unique(Int.(c2_values))
    isempty(requested_l2) && throw(ArgumentError("l2_values must not be empty"))
    isempty(requested_c2) && throw(ArgumentError("c2_values must not be empty"))

    catalog, rejected = level_catalog(
        states; quantum_tol=quantum_tol, degeneracy_tol=degeneracy_tol,
    )
    ground = minimum(state.energy for state in states)
    data = DataFrame()
    for c2 in requested_c2, l2 in requested_l2
        levels = get(catalog, (l2, c2), PhysicalLevel[])
        length(levels) >= levels_per_block || throw(ArgumentError(
            "Spectrum block (L2,C2)=($l2,$c2) has only $(length(levels)) distinct " *
            "levels; need $levels_per_block. Increase [spectrum].k.",
        ))
        column = "L2=$l2 C2=$c2"
        data[!, column] = [
            (levels[index].energy - ground) / Float64(factor)
            for index in 1:Int(levels_per_block)
        ]
    end
    return (data=data, catalog=catalog, rejected=rejected, ground=ground)
end

"""
在 `[hamiltonian].mu` 的单个参数点计算精简的 rescaled spectrum。

本征向量只在求解过程中用于识别 `L²/C₂`，不会写入磁盘。结果固定写成
`spectrum.csv`，与 `generator` 的 ED 快照和 tower 输出完全分开。
"""
function run_spectrum(
    nm1::Int,
    couplings::Couplings=Couplings(),
    settings::SolverSettings=SolverSettings();
    l2_values::AbstractVector{<:Integer}=[0, 2, 6],
    c2_values::AbstractVector{<:Integer}=[0, 3],
    levels_per_block::Integer=7,
    factor::Real,
    output::AbstractString=joinpath(PACKAGE_ROOT, "output", "spectrum"),
    force::Bool=false,
)
    output = ensure_output(output)
    scale = Float64(factor)
    isfinite(scale) && scale > 0 || throw(ArgumentError("factor must be positive and finite"))
    write_run_metadata(output; command="spectrum")
    path = joinpath(output, "spectrum.csv")
    if isfile(path) && !force
        @info "reusing spectrum" path
        return CSV.read(path, DataFrame)
    end

    @info "computing spectrum" nm1 mu=couplings.mu k=settings.k
    model = build_model(nm1=nm1)
    cache = prepare_spectrum(model, couplings, settings)
    states = solve_spectrum(cache, couplings.mu)
    result = spectrum_level_table(
        states; l2_values=l2_values, c2_values=c2_values,
        levels_per_block=levels_per_block, factor=scale,
        quantum_tol=settings.quantum_tol,
        degeneracy_tol=settings.degeneracy_tol,
    )
    !isempty(result.rejected) && @warn(
        "some states had non-integer L2/C2 and were excluded",
        rejected=length(result.rejected),
    )
    atomic_csv(path, result.data)
    @info "saved spectrum" path factor=scale rows=nrow(result.data)
    return result.data
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

"""旧式单区间有界 Brent 搜索；保留作底层对照，FSS 不再依赖其全局性。"""
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

_finite_objective(score) = score.valid && isfinite(score.objective)
_mu_key(mu) = round(Float64(mu); digits=14)

"""
找出发现网格上的局部最小值，并把连续平台合并成一个候选谷底。

无效点按 `Inf` 处理，因此每段有效区间的端点也可以成为候选。全局网格最小值
始终在返回结果中，避免严格局部极小判据因数值平台而漏掉当前最佳点。
"""
function _grid_local_minimum_indices(scores)
    objectives = [
        _finite_objective(score) ? Float64(score.objective) : Inf for score in scores
    ]
    valid = findall(isfinite, objectives)
    isempty(valid) && return Int[]
    candidates = Int[]
    for index in eachindex(objectives)
        isfinite(objectives[index]) || continue
        left = index == firstindex(objectives) ? Inf : objectives[index-1]
        right = index == lastindex(objectives) ? Inf : objectives[index+1]
        objectives[index] <= left && objectives[index] <= right && push!(candidates, index)
    end

    # 平坦谷底上的相邻点只需要精修一次；取平台中点使 bracket 尽量对称。
    collapsed = Int[]
    cursor = 1
    while cursor <= length(candidates)
        stop = cursor
        while stop < length(candidates) && candidates[stop+1] == candidates[stop] + 1
            stop += 1
        end
        group = candidates[cursor:stop]
        minimum_value = minimum(objectives[group])
        minima = [index for index in group if objectives[index] == minimum_value]
        push!(collapsed, minima[cld(length(minima), 2)])
        cursor = stop + 1
    end

    global_index = valid[argmin(objectives[valid])]
    global_index in collapsed || push!(collapsed, global_index)
    sort!(collapsed)
    return collapsed
end

function _new_scalar_search(
    score_at;
    penalty::Real=1.0e3,
    on_evaluation=(mu, score, source) -> nothing,
)
    evaluations = Dict{Float64,Any}()
    sources = Dict{Float64,String}()
    function evaluate(mu, source)
        key = _mu_key(mu)
        if !haskey(evaluations, key)
            score = score_at(key)
            evaluations[key] = score
            sources[key] = String(source)
            on_evaluation(key, score, String(source))
        end
        score = evaluations[key]
        return _finite_objective(score) ? Float64(score.objective) : Float64(penalty)
    end
    return (evaluate=evaluate, evaluations=evaluations, sources=sources)
end

"""在一个区间内做发现网格，并分别精修该网格识别出的局部谷底。"""
function _refine_grid_interval!(
    search;
    lower::Real,
    upper::Real,
    count::Int,
    abs_tol::Real,
    max_iterations::Int,
    source_prefix::AbstractString,
)
    lower < upper || throw(ArgumentError("search lower bound must be below upper bound"))
    count >= 3 || throw(ArgumentError("grid refinement needs at least three points"))
    grid = collect(range(Float64(lower), Float64(upper); length=count))
    grid_scores = Any[]
    for (index, mu) in enumerate(grid)
        @info "FSS mu discovery grid" source_prefix index total=count mu lower upper
        search.evaluate(mu, "$(source_prefix)_grid")
        push!(grid_scores, search.evaluations[_mu_key(mu)])
    end
    candidates = _grid_local_minimum_indices(grid_scores)
    attempted = 0
    converged = 0
    for (candidate_number, index) in enumerate(candidates)
        # 窗口边界没有双侧 bracket；调用者会据此决定是否扩大窗口。
        (index == firstindex(grid) || index == lastindex(grid)) && continue
        bracket_lower, bracket_upper = grid[index-1], grid[index+1]
        source = "$(source_prefix)_refine_$(candidate_number)"
        attempted += 1
        @info "FSS refining discovered valley" source index bracket_lower bracket_upper grid_mu=grid[index]
        try
            result = optimize(
                mu -> search.evaluate(mu, source), bracket_lower, bracket_upper, Brent();
                abs_tol=Float64(abs_tol), rel_tol=1.0e-5,
                iterations=max_iterations, show_trace=false,
            )
            search.evaluate(Optim.minimizer(result), source)
            converged += Optim.converged(result)
        catch err
            @warn "FSS valley refinement failed; retaining evaluated points" source bracket_lower bracket_upper error=sanitize_error(err)
        end
    end
    valid = findall(_finite_objective, grid_scores)
    best_index = isempty(valid) ? nothing : valid[argmin([
        Float64(grid_scores[index].objective) for index in valid
    ])]
    return (
        grid=grid, scores=grid_scores, candidates=candidates,
        candidate_count=length(candidates), refined_count=attempted,
        refinements_converged=converged, grid_best_index=best_index,
        grid_best_mu=best_index === nothing ? NaN : grid[best_index],
        grid_best_objective=best_index === nothing ? Inf :
                            Float64(grid_scores[best_index].objective),
        grid_best_at_boundary=best_index === nothing ? false :
                              best_index == firstindex(grid) || best_index == lastindex(grid),
    )
end

function _scalar_search_summary(search; mu_min::Real, mu_max::Real, abs_tol::Real)
    ordered_keys = sort!(collect(keys(search.evaluations)))
    valid_keys = [key for key in ordered_keys if _finite_objective(search.evaluations[key])]
    records = [
        (mu=key, score=search.evaluations[key], source=search.sources[key])
        for key in ordered_keys
    ]
    if isempty(valid_keys)
        return (
            mu=NaN, score=nothing, at_boundary=false,
            evaluations=length(ordered_keys), invalid=length(ordered_keys),
            best_source="none", evaluation_records=records,
        )
    end
    best_key = valid_keys[argmin([
        Float64(search.evaluations[key].objective) for key in valid_keys
    ])]
    tolerance = max(Float64(abs_tol), 1.0e-8)
    return (
        mu=best_key, score=search.evaluations[best_key],
        at_boundary=abs(best_key-mu_min) <= tolerance || abs(best_key-mu_max) <= tolerance,
        evaluations=length(ordered_keys),
        invalid=count(key -> !_finite_objective(search.evaluations[key]), ordered_keys),
        best_source=search.sources[best_key], evaluation_records=records,
    )
end

"""宽范围 anchor 搜索：完整发现网格，并精修网格中的每一个局部谷底。"""
function _global_grid_refine(
    score_at;
    mu_min::Real,
    mu_max::Real,
    mu_count::Int,
    abs_tol::Real=1.0e-4,
    max_iterations::Int=60,
    penalty::Real=1.0e3,
    on_evaluation=(mu, score, source) -> nothing,
)
    mu_min < mu_max || throw(ArgumentError("mu_min must be smaller than mu_max"))
    mu_count >= 3 || throw(ArgumentError("global grid-refine search needs mu_count >= 3"))
    abs_tol > 0 || throw(ArgumentError("abs_tol must be positive"))
    max_iterations > 0 || throw(ArgumentError("max_iterations must be positive"))
    search = _new_scalar_search(
        score_at; penalty=penalty, on_evaluation=on_evaluation,
    )
    interval = _refine_grid_interval!(
        search; lower=mu_min, upper=mu_max, count=mu_count,
        abs_tol=abs_tol, max_iterations=max_iterations, source_prefix="anchor",
    )
    summary = _scalar_search_summary(
        search; mu_min=mu_min, mu_max=mu_max, abs_tol=abs_tol,
    )
    return merge(summary, (
        completed=interval.refinements_converged == interval.refined_count,
        grid_points=mu_count, grid_passes=1,
        candidate_count=interval.candidate_count,
        refined_count=interval.refined_count,
        refinements_converged=interval.refinements_converged,
        grid_best_mu=interval.grid_best_mu,
        grid_best_objective=interval.grid_best_objective,
        search_mode="anchor", center_mu=NaN,
        search_lower=Float64(mu_min), search_upper=Float64(mu_max),
        expansions=0, wide_brent_mu=NaN, wide_brent_objective=Inf,
        wide_brent_converged=false,
        wide_brent_ran=false, wide_guard_reason="",
        local_mu=NaN, local_objective=Inf, local_invalid=0,
        wide_disagreement=false, wide_disagreement_reason="",
        candidate_mus=interval.grid[interval.candidates],
    ))
end

"""
沿前一尺寸的 μc 做局部细搜；边界最低会自动扩窗。

`wide_mode=:always` 保留原来的宽区间 Brent 候选对照；`wide_mode=:adaptive`
只在审计点或局部搜索出现明显异常时运行宽区间 challenger。
"""
function _continuation_grid_refine(
    score_at;
    center::Real,
    mu_min::Real,
    mu_max::Real,
    local_half_width::Real,
    local_count::Int,
    max_expansions::Int,
    abs_tol::Real=1.0e-4,
    max_iterations::Int=60,
    penalty::Real=1.0e3,
    wide_mode::Symbol=:always,
    force_wide::Bool=false,
    force_wide_reason::AbstractString="",
    wide_jump_tol::Real=0.03,
    wide_mu_tol::Real=5.0e-3,
    wide_objective_tol::Real=1.0e-4,
    on_evaluation=(mu, score, source) -> nothing,
)
    mu_min < mu_max || throw(ArgumentError("mu_min must be smaller than mu_max"))
    isfinite(center) || throw(ArgumentError("continuation center must be finite"))
    local_half_width > 0 || throw(ArgumentError("local_half_width must be positive"))
    local_count >= 3 || throw(ArgumentError("local continuation needs at least three grid points"))
    max_expansions >= 0 || throw(ArgumentError("max_expansions must be nonnegative"))
    wide_mode in (:always, :adaptive) || throw(ArgumentError(
        "wide_mode must be always or adaptive",
    ))
    wide_jump_tol > 0 || throw(ArgumentError("wide_jump_tol must be positive"))
    wide_mu_tol > 0 || throw(ArgumentError("wide_mu_tol must be positive"))
    wide_objective_tol >= 0 || throw(ArgumentError(
        "wide_objective_tol must be nonnegative",
    ))
    search = _new_scalar_search(
        score_at; penalty=penalty, on_evaluation=on_evaluation,
    )
    clipped_center = clamp(Float64(center), Float64(mu_min), Float64(mu_max))
    half_width = Float64(local_half_width)
    expansion = 0
    intervals = Any[]
    while true
        lower = max(Float64(mu_min), clipped_center-half_width)
        upper = min(Float64(mu_max), clipped_center+half_width)
        interval = _refine_grid_interval!(
            search; lower=lower, upper=upper, count=local_count,
            abs_tol=abs_tol, max_iterations=max_iterations,
            source_prefix="local_$(expansion)",
        )
        push!(intervals, interval)
        full_range = lower <= mu_min && upper >= mu_max
        needs_expansion = interval.grid_best_index === nothing ||
                          interval.grid_best_at_boundary
        (!needs_expansion || full_range || expansion >= max_expansions) && break
        expansion += 1
        half_width *= 2
        @info "FSS local minimum reached window edge; expanding" center=clipped_center expansion half_width
    end

    local_summary = _scalar_search_summary(
        search; mu_min=mu_min, mu_max=mu_max, abs_tol=abs_tol,
    )
    final_interval = last(intervals)
    wide_guard_reasons = String[]
    if wide_mode == :always
        push!(wide_guard_reasons, "always")
    else
        if force_wide
            push!(
                wide_guard_reasons,
                isempty(force_wide_reason) ? "forced_audit" : String(force_wide_reason),
            )
        end
        if local_summary.score === nothing
            push!(wide_guard_reasons, "no_valid_local_score")
        end
        if final_interval.grid_best_index === nothing || final_interval.grid_best_at_boundary
            push!(wide_guard_reasons, "unresolved_local_boundary")
        end
        if local_summary.score !== nothing &&
           abs(local_summary.mu-clipped_center) > Float64(wide_jump_tol)
            push!(wide_guard_reasons, "large_muc_jump")
        end
    end
    unique!(wide_guard_reasons)
    wide_ran = !isempty(wide_guard_reasons)
    severe_guard_reasons = filter(
        reason -> reason in (
            "no_valid_local_score", "unresolved_local_boundary", "large_muc_jump",
        ),
        wide_guard_reasons,
    )
    if !isempty(severe_guard_reasons)
        @warn "FSS adaptive guard detected a suspicious local muc search; running the wide challenger" center=clipped_center reasons=join(severe_guard_reasons, ";")
    end
    wide_mu = NaN
    wide_objective = Inf
    wide_converged = false
    if wide_ran
        try
            wide = optimize(
                mu -> search.evaluate(mu, "wide_brent"),
                Float64(mu_min), Float64(mu_max), Brent();
                abs_tol=Float64(abs_tol), rel_tol=1.0e-5,
                iterations=max_iterations, show_trace=false,
            )
            wide_mu = Float64(Optim.minimizer(wide))
            search.evaluate(wide_mu, "wide_brent")
            wide_score = search.evaluations[_mu_key(wide_mu)]
            wide_objective = _finite_objective(wide_score) ? wide_score.objective : Inf
            wide_converged = Optim.converged(wide)
        catch err
            @warn "FSS wide Brent challenger failed; retaining continuation points" error=sanitize_error(err)
        end
    end

    summary = _scalar_search_summary(
        search; mu_min=mu_min, mu_max=mu_max, abs_tol=abs_tol,
    )
    disagreement_reasons = String[]
    if wide_ran && local_summary.score !== nothing && isfinite(wide_objective)
        abs(local_summary.mu-wide_mu) > Float64(wide_mu_tol) &&
            push!(disagreement_reasons, "muc")
        abs(Float64(local_summary.score.objective)-Float64(wide_objective)) >
            Float64(wide_objective_tol) && push!(disagreement_reasons, "objective")
    elseif wide_ran && (local_summary.score !== nothing) != isfinite(wide_objective)
        push!(disagreement_reasons, "validity")
    end
    refined_count = sum(interval.refined_count for interval in intervals)
    refinements_converged = sum(interval.refinements_converged for interval in intervals)
    return merge(summary, (
        completed=(!wide_ran || wide_converged) && refinements_converged == refined_count,
        grid_points=local_count, grid_passes=length(intervals),
        candidate_count=sum(interval.candidate_count for interval in intervals),
        refined_count=refined_count,
        refinements_converged=refinements_converged,
        grid_best_mu=final_interval.grid_best_mu,
        grid_best_objective=final_interval.grid_best_objective,
        search_mode="continuation", center_mu=clipped_center,
        search_lower=first(final_interval.grid), search_upper=last(final_interval.grid),
        expansions=expansion, wide_brent_mu=wide_mu,
        wide_brent_objective=wide_objective,
        wide_brent_converged=wide_converged,
        wide_brent_ran=wide_ran,
        wide_guard_reason=join(wide_guard_reasons, ";"),
        local_mu=local_summary.mu,
        local_objective=local_summary.score === nothing ? Inf :
                        Float64(local_summary.score.objective),
        local_invalid=local_summary.invalid,
        wide_disagreement=!isempty(disagreement_reasons),
        wide_disagreement_reason=join(disagreement_reasons, ";"),
        candidate_mus=final_interval.grid[final_interval.candidates],
    ))
end

"""FSS 专用：宽范围 anchor，或沿前一尺寸 μc 的 continuation 搜索。"""
function optimize_fss_mu_with_score(
    cache::ModelCache;
    mu_min::Real,
    mu_max::Real,
    mu_count::Int,
    center=nothing,
    local_half_width::Real=0.02,
    local_count::Int=9,
    max_expansions::Int=3,
    definition=:fss7,
    terms=nothing,
    metric=nothing,
    abs_tol::Real=1.0e-4,
    max_iterations::Int=60,
    penalty::Real=1.0e3,
    wide_mode::Symbol=:always,
    force_wide::Bool=false,
    force_wide_reason::AbstractString="",
    wide_jump_tol::Real=0.03,
    wide_mu_tol::Real=5.0e-3,
    wide_objective_tol::Real=1.0e-4,
    on_evaluation=(mu, score, source) -> nothing,
)
    definition = normalize_score_definition(definition)
    selected_terms = resolve_score_terms(definition, terms)
    metric = normalize_score_metric(metric, definition)
    score_at(mu) = try
        cft_score(
            solve_spectrum(cache, mu); settings=cache.settings,
            definition=definition, terms=selected_terms, metric=metric,
        )
    catch err
        @error "FSS mu evaluation failed" mu exception=(err, catch_backtrace())
        _invalid_score(definition, metric, sanitize_error(err); terms=selected_terms)
    end
    result = if center === nothing
        _global_grid_refine(
            score_at; mu_min=mu_min, mu_max=mu_max, mu_count=mu_count,
            abs_tol=abs_tol, max_iterations=max_iterations, penalty=penalty,
            on_evaluation=on_evaluation,
        )
    else
        _continuation_grid_refine(
            score_at; center=Float64(center), mu_min=mu_min, mu_max=mu_max,
            local_half_width=local_half_width, local_count=local_count,
            max_expansions=max_expansions, abs_tol=abs_tol,
            max_iterations=max_iterations, penalty=penalty,
            wide_mode=wide_mode, force_wide=force_wide,
            force_wide_reason=force_wide_reason,
            wide_jump_tol=wide_jump_tol, wide_mu_tol=wide_mu_tol,
            wide_objective_tol=wide_objective_tol,
            on_evaluation=on_evaluation,
        )
    end
    if result.score === nothing
        return merge(result, (
            score=_invalid_score(
                definition, metric, "no valid score in FSS global search";
                terms=selected_terms,
            ),
        ))
    end
    return result
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

function _fss_scan_couplings(
    base::Couplings,
    fss::FSSSettings,
    scan_value::Real;
    u0_over_uf::Union{Nothing,Real}=nothing,
)
    couplings = with_coupling(base, fss.scan_parameter, scan_value)
    if u0_over_uf !== nothing
        fss.scan_parameter == :Uf || throw(ArgumentError(
            "u0_over_uf can only be used when FSS scan_parameter is Uf",
        ))
        ratio = Float64(u0_over_uf)
        isfinite(ratio) || throw(ArgumentError("u0_over_uf must be finite"))
        couplings = with_coupling(couplings, :U0, ratio * couplings.Uf)
    end
    return validate(couplings)
end

"""
完整 finite-size scaling 数据生成流程。

对每个 `nm1` 和外层 `scan_value`，`grid` 在固定 μ 网格上选最小点；
`optimize` 对每个固定外层参数沿 size 追踪 μc：小尺寸建立 anchor，大尺寸局部续接。
两种结果分文件保存。
"""
function run_fss_scan(
    base::Couplings,
    fss::FSSSettings,
    settings::SolverSettings=SolverSettings();
    output::AbstractString=joinpath(PACKAGE_ROOT, "output", "fss"),
    force::Bool=false,
    u0_over_uf::Union{Nothing,Real}=nothing,
)
    fss.scan_parameter in HAMILTONIAN_FIELDS || throw(ArgumentError("Unknown scan parameter $(fss.scan_parameter)"))
    fss.scan_parameter == :mu && throw(ArgumentError("FSS scan_parameter cannot be mu; mu already has its own grid"))
    fss.mu_count > 0 || throw(ArgumentError("FSS mu_count must be positive"))
    methods = unique(fss.methods)
    isempty(methods) && throw(ArgumentError("FSS methods must not be empty"))
    all(method -> method in (:grid, :optimize), methods) ||
        throw(ArgumentError("FSS methods must contain only grid and/or optimize"))
    if :optimize in methods
        fss.optimize_strategy == :size_continuation || throw(ArgumentError(
            "FSS optimize_strategy must be size_continuation",
        ))
        fss.mu_count >= 3 || throw(ArgumentError(
            "FSS size_continuation needs mu_count >= 3 for anchor searches",
        ))
        fss.optimize_anchor_nm > 0 || throw(ArgumentError("FSS optimize_anchor_nm must be positive"))
        fss.optimize_local_half_width > 0 || throw(ArgumentError("FSS optimize_local_half_width must be positive"))
        fss.optimize_local_count >= 3 || throw(ArgumentError("FSS optimize_local_count must be at least 3"))
        fss.optimize_max_expansions >= 0 || throw(ArgumentError("FSS optimize_max_expansions must be nonnegative"))
        fss.optimize_wide_mode in (:always, :adaptive) || throw(ArgumentError(
            "FSS optimize_wide_mode must be always or adaptive",
        ))
        if fss.optimize_wide_mode == :adaptive
            fss.optimize_wide_adaptive_nm > fss.optimize_anchor_nm || throw(ArgumentError(
                "FSS optimize_wide_adaptive_nm must be larger than optimize_anchor_nm",
            ))
        end
        fss.optimize_wide_jump_tol > 0 || throw(ArgumentError(
            "FSS optimize_wide_jump_tol must be positive",
        ))
        fss.optimize_wide_mu_tol > 0 || throw(ArgumentError(
            "FSS optimize_wide_mu_tol must be positive",
        ))
        fss.optimize_wide_objective_tol >= 0 || throw(ArgumentError(
            "FSS optimize_wide_objective_tol must be nonnegative",
        ))
    end
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
    saved_muc = Dict{Tuple{Int,Float64},Float64}()
    wide_fallback_sizes = Set{Int}()
    for method in methods
        path = paths[method]
        if isfile(path)
            previous = CSV.read(path, DataFrame)
            "score_definition" in names(previous) || error(
                "Existing $method FSS data has no score label; use a new FSS case directory",
            )
            all(String(value) == definition_tag for value in previous.score_definition) ||
                error("Existing $method FSS data uses a different score; use a new FSS case directory")
            if method == :optimize && !force
                for row in eachrow(latest_rows(previous))
                    if String(row.status) == "ok" && Bool(row.score_valid) &&
                       isfinite(row.muc)
                        saved_muc[(Int(row.nm1), Float64(row.scan_value))] = Float64(row.muc)
                    end
                    if fss.optimize_wide_mode == :adaptive &&
                       "wide_disagreement" in names(previous) &&
                       !ismissing(row.wide_disagreement) && Bool(row.wide_disagreement)
                        push!(wide_fallback_sizes, Int(row.nm1))
                    end
                end
            end
        end
        completed[method] = force ? Set{String}() : completed_job_ids(path)
    end

    # 为复用同一 size 的 Basis/H0 缓存，执行次序仍是 size 在外；这个字典按
    # scan_value 分开保存上一 size 的 μc，物理上等价于逐条固定参数线做 continuation。
    previous_size_muc = Dict{Float64,Float64}()
    for nm1 in sort(unique(fss.nm_values))
        model = build_model(nm1=nm1)
        cache = nothing
        for (scan_index, scan_value) in enumerate(fss.scan_values)
            couplings = _fss_scan_couplings(
                base, fss, scan_value; u0_over_uf=u0_over_uf,
            )
            job_ids = Dict(method => stable_id(
                "fss-$method-v4", definition_tag, nm1, fss.scan_parameter, scan_value,
                coupling_vector(couplings), settings, fss.mu_min, fss.mu_max,
                fss.mu_count,
                method == :grid ? :fixed_grid : (
                    fss.optimize_strategy, fss.optimize_anchor_nm,
                    fss.optimize_local_half_width, fss.optimize_local_count,
                    fss.optimize_max_expansions, fss.optimize_abs_tol,
                    fss.optimize_max_iterations, fss.optimize_wide_mode,
                    fss.optimize_wide_adaptive_nm, fss.optimize_wide_audit_first,
                    fss.optimize_wide_audit_all,
                    fss.optimize_wide_jump_tol, fss.optimize_wide_mu_tol,
                    fss.optimize_wide_objective_tol,
                ),
            ) for method in methods)
            if :optimize in methods && job_ids[:optimize] in completed[:optimize]
                key = (nm1, Float64(scan_value))
                if nm1 >= fss.optimize_anchor_nm && haskey(saved_muc, key)
                    previous_size_muc[Float64(scan_value)] = saved_muc[key]
                end
            end
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
                        center = nm1 > fss.optimize_anchor_nm ?
                                 get(previous_size_muc, Float64(scan_value), nothing) : nothing
                        center === nothing && nm1 > fss.optimize_anchor_nm && @warn(
                            "previous-size muc unavailable; using a wide anchor search",
                            nm1, scan_value,
                        )
                        trace_path = joinpath(output, "fss_optimize_evaluations.csv")
                        adaptive_wide = fss.optimize_wide_mode == :adaptive &&
                                        nm1 >= fss.optimize_wide_adaptive_nm
                        fallback_active = adaptive_wide && nm1 in wide_fallback_sizes
                        audit_point = adaptive_wide && (
                            fss.optimize_wide_audit_all ||
                            (fss.optimize_wide_audit_first &&
                             scan_index == firstindex(fss.scan_values))
                        )
                        force_wide = fallback_active || audit_point
                        force_wide_reason = fallback_active ?
                            "prior_audit_disagreement" :
                            (audit_point ? "first_point_audit" : "")
                        evaluation_number = Ref(0)
                        on_evaluation = function(mu, score, source)
                            evaluation_number[] += 1
                            append_csv(trace_path, (
                                evaluation_id=stable_id("fss-mu-evaluation-v1", job_id, mu),
                                job_id=job_id, evaluation=evaluation_number[],
                                status="ok", timestamp=string(now()), nm1=nm1,
                                scan_parameter=String(fss.scan_parameter),
                                scan_value=scan_value, mu=mu, source=source,
                                center_muc=center === nothing ? NaN : center,
                                objective=score.objective, q=score.q, cost=score.cost,
                                factor=score.factor, delta_s=score.delta_s,
                                delta_o=score.delta_o, score_valid=score.valid,
                                score_reason=score.reason,
                                raw_gaps=join(score.raw_gaps, ";"),
                                target_gaps=join(score.target_gaps, ";"),
                                labels=join(score.labels, ";"),
                            ))
                        end
                        optimize_fss_mu_with_score(
                            cache; mu_min=fss.mu_min, mu_max=fss.mu_max,
                            mu_count=fss.mu_count,
                            center=center,
                            local_half_width=fss.optimize_local_half_width,
                            local_count=fss.optimize_local_count,
                            max_expansions=fss.optimize_max_expansions,
                            definition=definition, terms=selected_terms, metric=metric,
                            abs_tol=fss.optimize_abs_tol,
                            max_iterations=fss.optimize_max_iterations,
                            wide_mode=adaptive_wide ? :adaptive : :always,
                            force_wide=force_wide,
                            force_wide_reason=force_wide_reason,
                            wide_jump_tol=fss.optimize_wide_jump_tol,
                            wide_mu_tol=fss.optimize_wide_mu_tol,
                            wide_objective_tol=fss.optimize_wide_objective_tol,
                            on_evaluation=on_evaluation,
                        )
                    end
                    score = result.score
                    search_strategy = method == :grid ? "fixed_grid" : String(fss.optimize_strategy)
                    candidate_count = method == :grid ? 0 : result.candidate_count
                    refined_count = method == :grid ? 0 : result.refined_count
                    refinements_converged = method == :grid ? 0 : result.refinements_converged
                    grid_best_muc = method == :grid ? result.mu : result.grid_best_mu
                    grid_best_objective = method == :grid ? score.objective : result.grid_best_objective
                    best_source = method == :grid ? "grid" : result.best_source
                    search_mode = method == :grid ? "grid" : result.search_mode
                    center_muc = method == :grid ? NaN : result.center_mu
                    search_lower = method == :grid ? fss.mu_min : result.search_lower
                    search_upper = method == :grid ? fss.mu_max : result.search_upper
                    grid_points = method == :grid ? fss.mu_count : result.grid_points
                    grid_passes = method == :grid ? 1 : result.grid_passes
                    expansions = method == :grid ? 0 : result.expansions
                    wide_brent_muc = method == :grid ? NaN : result.wide_brent_mu
                    wide_brent_objective = method == :grid ? Inf : result.wide_brent_objective
                    wide_brent_converged = method == :grid ? false : result.wide_brent_converged
                    wide_brent_ran = method == :grid ? false : result.wide_brent_ran
                    wide_guard_reason = method == :grid ? "" : result.wide_guard_reason
                    local_muc = method == :grid ? NaN : result.local_mu
                    local_objective = method == :grid ? Inf : result.local_objective
                    local_invalid = method == :grid ? 0 : result.local_invalid
                    wide_disagreement = method == :grid ? false : result.wide_disagreement
                    wide_disagreement_reason = method == :grid ? "" :
                                               result.wide_disagreement_reason
                    append_csv(path, (
                        score_definition=definition_tag, score_terms=terms_tag,
                        method=String(method), search_strategy=search_strategy,
                        job_id=job_id, status="ok", timestamp=string(now()), nm1=nm1,
                        x=nm1^(-0.5), scan_parameter=String(fss.scan_parameter),
                        scan_value=scan_value, muc=result.mu, objective=score.objective,
                        q=score.q, cost=score.cost, factor=score.factor,
                        delta_s=score.delta_s, delta_o=score.delta_o,
                        score_valid=score.valid, score_reason=score.reason,
                        completed=result.completed, at_boundary=result.at_boundary,
                        evaluations=result.evaluations, invalid=result.invalid,
                        search_mode=search_mode, center_muc=center_muc,
                        search_lower=search_lower, search_upper=search_upper,
                        grid_points=grid_points, grid_passes=grid_passes,
                        candidate_count=candidate_count,
                        refined_count=refined_count,
                        refinements_converged=refinements_converged,
                        expansions=expansions,
                        wide_brent_muc=wide_brent_muc,
                        wide_brent_objective=wide_brent_objective,
                        wide_brent_converged=wide_brent_converged,
                        wide_brent_ran=wide_brent_ran,
                        wide_guard_reason=wide_guard_reason,
                        local_muc=local_muc,
                        local_objective=local_objective,
                        local_invalid=local_invalid,
                        wide_disagreement=wide_disagreement,
                        wide_disagreement_reason=wide_disagreement_reason,
                        grid_best_muc=grid_best_muc,
                        grid_best_objective=grid_best_objective,
                        best_source=best_source,
                        error="", coupling_namedtuple(couplings)...,
                    ))
                    if method == :optimize && score.valid && isfinite(result.mu) &&
                       nm1 >= fss.optimize_anchor_nm
                        previous_size_muc[Float64(scan_value)] = result.mu
                    end
                    if method == :optimize &&
                       fss.optimize_wide_mode == :adaptive &&
                       nm1 >= fss.optimize_wide_adaptive_nm &&
                       result.wide_disagreement
                        push!(wide_fallback_sizes, nm1)
                        @warn "FSS local/wide audit disagreed; remaining points at this size will keep the wide challenger" nm1 scan_value local_muc wide_brent_muc local_objective wide_brent_objective reason=wide_disagreement_reason
                    end
                catch err
                    @error "FSS job failed" method nm1 scan_value exception=(err, catch_backtrace())
                    append_csv(path, (
                        score_definition=definition_tag, score_terms=terms_tag,
                        method=String(method),
                        search_strategy=method == :grid ? "fixed_grid" : String(fss.optimize_strategy),
                        job_id=job_id, status="error", timestamp=string(now()), nm1=nm1,
                        x=nm1^(-0.5), scan_parameter=String(fss.scan_parameter),
                        scan_value=scan_value, muc=NaN, objective=Inf, q=Inf,
                        cost=Inf, factor=NaN, delta_s=NaN, delta_o=NaN,
                        score_valid=false, score_reason="", completed=false,
                        at_boundary=false, evaluations=0, invalid=0,
                        search_mode="", center_muc=NaN,
                        search_lower=NaN, search_upper=NaN,
                        grid_points=0, grid_passes=0, candidate_count=0,
                        refined_count=0, refinements_converged=0,
                        expansions=0, wide_brent_muc=NaN,
                        wide_brent_objective=Inf,
                        wide_brent_converged=false,
                        wide_brent_ran=false, wide_guard_reason="",
                        local_muc=NaN, local_objective=Inf, local_invalid=0,
                        wide_disagreement=false, wide_disagreement_reason="",
                        grid_best_muc=NaN, grid_best_objective=Inf,
                        best_source="",
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
