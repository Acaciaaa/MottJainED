# 本文件把有限尺寸低能谱压缩成一个“像不像目标 CFT tower”的分数，
# 并据此优化化学势 μ。它不负责构造 Hamiltonian。

function _require_levels(catalog, key::Tuple{Int,Int}, count::Int)
    levels = get(catalog, key, PhysicalLevel[])
    length(levels) >= count || return nothing
    return levels
end

"""
按照原分析中的七个低能 CFT tower 关系评价一组本征态。

函数先合并离散对称性产生的同能副本，再选 `(L²,C₂)` 能级。`factor` 是用
最小二乘得到的有限尺寸能量到标度维数的换算因子；`q` 是七个关系在标度
维数单位下的 RMS 误差，因此 `q` 越小越符合目标 CFT tower。
"""
function cft_score(states::Vector{SpectrumState}; settings::SolverSettings=SolverSettings())
    isempty(states) && return CFTScore(reason="empty spectrum")
    catalog, rejected = level_catalog(
        states; quantum_tol=settings.quantum_tol,
        degeneracy_tol=settings.degeneracy_tol,
    )
    # 每个变量是对应 (L²,C₂) sector 中最低的若干个“不同能量”物理能级。
    s00 = _require_levels(catalog, (0, 0), 3)
    s20 = _require_levels(catalog, (2, 0), 1)
    s60 = _require_levels(catalog, (6, 0), 2)
    a23 = _require_levels(catalog, (2, 3), 2)
    a63 = _require_levels(catalog, (6, 3), 1)
    required = (s00=s00, s20=s20, s60=s60, a23=a23, a63=a63)
    missing = [String(name) for (name, value) in pairs(required) if isnothing(value)]
    if !isempty(missing)
        detail = isempty(rejected) ? "" : "; $(length(rejected)) states had non-integer quantum numbers"
        return CFTScore(reason="missing levels: $(join(missing, ", "))$detail")
    end

    ground = minimum(state.energy for state in states)
    # 七个有限尺寸原始能隙，顺序必须与 targets/labels 一一对应。
    raw = Float64[
        s20[1].energy - s00[2].energy, # dS - S
        s60[2].energy - s20[1].energy, # ddS - dS
        s00[3].energy - s00[2].energy, # boxS - S
        a23[1].energy - ground,         # J
        a23[2].energy - ground,         # curl J
        a63[1].energy - ground,         # dJ
        s60[1].energy - ground,         # T
    ]
    targets = Float64[1, 1, 2, 2, 3, 3, 3]
    labels = ["dS-S", "ddS-dS", "boxS-S", "J", "curlJ", "dJ", "T"]
    all(isfinite, raw) || return CFTScore(reason="non-finite tower gap")
    # 拟合 raw ≈ factor * targets；factor 吸收非普适的速度/球半径尺度。
    factor = dot(raw, targets) / dot(targets, targets)
    factor > eps(Float64) || return CFTScore(reason="non-positive fitted energy factor")
    residual = raw ./ factor .- targets
    q = sqrt(mean(abs2, residual))
    delta_s = (s00[2].energy - ground) / factor
    a03 = get(catalog, (0, 3), PhysicalLevel[])
    delta_o = isempty(a03) ? NaN : (a03[1].energy - ground) / factor
    return CFTScore(
        valid=true, q=q, factor=factor, delta_s=delta_s, delta_o=delta_o,
        raw_gaps=raw, target_gaps=targets, labels=labels,
    )
end

"""
在一个缓存好的 Hamiltonian family 中寻找使 `cft_score.q` 最小的 μ。

先在 `[mu_min,mu_max]` 上粗网格扫描，找到最好点附近的 bracket，再用有界
Brent 方法细化。返回最佳 μ、最终分数、收敛/边界诊断及评估次数。
"""
function optimize_mu(
    cache::ModelCache;
    mu_min::Real=0.0,
    mu_max::Real=0.12,
    coarse_points::Int=9,
    abs_tol::Real=1.0e-5,
    max_iterations::Int=60,
    penalty::Real=1.0e3,
)
    mu_min < mu_max || throw(ArgumentError("mu_min must be smaller than mu_max"))
    coarse_points >= 3 || throw(ArgumentError("coarse_points must be at least 3"))
    evaluations = Dict{Float64,CFTScore}()
    errors = String[]

    # 同一个 μ 可能被优化器重复询问；Dict 避免重复对角化。
    function evaluate(mu)
        key = round(Float64(mu); digits=14)
        if !haskey(evaluations, key)
            evaluations[key] = try
                cft_score(solve_spectrum(cache, key); settings=cache.settings)
            catch err
                push!(errors, "mu=$key: $(sprint(showerror, err))")
                CFTScore(reason=sprint(showerror, err))
            end
        end
        score = evaluations[key]
        return score.valid && isfinite(score.q) ? score.q : Float64(penalty)
    end

    # 粗扫的目的不是给最终答案，而是避免 Brent 在错误的局部区间搜索。
    grid = collect(range(Float64(mu_min), Float64(mu_max); length=coarse_points))
    grid_costs = evaluate.(grid)
    best_index = argmin(grid_costs)
    left_index = max(1, best_index - 1)
    right_index = min(length(grid), best_index + 1)
    left_index == right_index && error("Could not form a local mu bracket")
    result = optimize(
        evaluate, grid[left_index], grid[right_index], Brent();
        abs_tol=Float64(abs_tol), rel_tol=0.0, iterations=max_iterations,
        show_trace=false,
    )
    mu = Float64(Optim.minimizer(result))
    final_states = solve_spectrum(cache, mu)
    score = cft_score(final_states; settings=cache.settings)
    boundary_tol = max(Float64(abs_tol), 1.0e-8)
    at_boundary = abs(mu - mu_min) <= boundary_tol || abs(mu - mu_max) <= boundary_tol
    return (
        mu=mu, score=score, converged=Optim.converged(result),
        at_boundary=at_boundary, evaluations=length(evaluations),
        invalid=count(value -> !value.valid, values(evaluations)),
        errors=errors, result=result,
    )
end

function score_dataframe(score::CFTScore)
    # 一行一个 tower relation，适合直接保存为诊断 CSV。
    return DataFrame(
        label=score.labels, raw_gap=score.raw_gaps,
        target=score.target_gaps,
        scaled_gap=score.raw_gaps ./ score.factor,
        residual=score.raw_gaps ./ score.factor .- score.target_gaps,
    )
end
