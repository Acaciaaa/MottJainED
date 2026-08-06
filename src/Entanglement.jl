# 本文件计算两种基态纠缠谱。共同步骤是先求带向量的基态，再定义 A/B 切分，
# 调用 FuzzifiED.GetEntSpec，最后统一保存 Schmidt 权重 λ 与 ξ=-log(λ)。

function _choose_sumsets(values::Vector{Int}, maximum_count::Int)
    # 动态规划：从给定单粒子 2Lz 列表中选 count 个时，所有可能的总 2Lz。
    sums = [Set{Int}() for _ in 0:maximum_count]
    push!(sums[1], 0)
    for value in values
        for count in maximum_count-1:-1:0
            for old in sums[count+1]
                push!(sums[count+2], old + value)
            end
        end
    end
    return sums
end

function _diagonal_sectors(values1::Vector{Int}, values0::Vector{Int})
    # 枚举子系统可能的 (charge,2Lz,F3,F8,0) diagonal quantum-number sectors。
    sums1 = _choose_sumsets(values1, length(values1))
    sums0 = _choose_sumsets(values0, length(values0))
    sectors = Set{NTuple{5,Int}}()
    for n11 in 0:length(values1), n12 in 0:length(values1),
        n13 in 0:length(values1), n0 in 0:length(values0)
        charge = n11 + n12 + n13 + 3n0
        f3 = n11 - n12
        f8 = n11 + n12 - 2n13
        for lz1 in sums1[n11+1], lz2 in sums1[n12+1],
            lz3 in sums1[n13+1], lz0 in sums0[n0+1]
            push!(sectors, (charge, lz1+lz2+lz3+lz0, f3, f8, 0))
        end
    end
    return sectors
end

function _ground_state(model, couplings, settings)
    # 纠缠谱必须知道基态的完整 Fock-basis coefficient，因此 keep_vectors=true。
    cache = prepare_spectrum(model, couplings, settings)
    states = solve_spectrum(cache, couplings.mu; keep_vectors=true)
    ground = first(states)
    return ground.vector, ground.basis, ground, cache
end

function _orbital_cut(model; nm1_a=model.nm1 ÷ 2, nm0_a=model.nm0 ÷ 2)
    # A 取两类粒子的前半球面轨道，B 取剩余轨道；amplitudes 是 0/1 硬切分。
    orbitals_a = Int[]
    for m in 0:nm1_a-1, flavor in 1:model.nf1
        push!(orbitals_a, m*model.nf1 + flavor)
    end
    append!(orbitals_a, model.no1 .+ collect(1:nm0_a))
    set_a = Set(orbitals_a)
    orbitals_b = [orbital for orbital in 1:model.no if orbital ∉ set_a]
    amplitudes = [orbital ∈ set_a ? 1 : 0 for orbital in 1:model.no]
    return orbitals_a, orbitals_b, amplitudes
end

function _entanglement_dataframe(entanglement)
    # GetEntSpec 的嵌套输出转成长表；λ 是 density-matrix eigenvalue。
    rows = NamedTuple[]
    for (sector, values) in entanglement
        for value in values
            value > 0 || continue
            sec = sector.secd_a
            push!(rows, (
                QA=sec[1], Lz2A=sec[2], F3A=length(sec) >= 3 ? sec[3] : 0,
                F8A=length(sec) >= 4 ? sec[4] : 0, lambda=Float64(value),
                xi=-log(Float64(value)),
            ))
        end
    end
    return DataFrame(rows)
end

function _save_entanglement(data, output, title; xi_cut=10.0)
    # 选择总 Schmidt 权重最大的 QA sector，平移最低 ξ，并统计 cutoff 以下 counting。
    nrow(data) > 0 || error("The entanglement spectrum is empty")
    weights = combine(groupby(data, :QA), :lambda => sum => :weight)
    qa = weights.QA[argmax(weights.weight)]
    selected = filter(row -> row.QA == qa, data)
    minimum_xi = minimum(selected.xi)
    selected.xi_shifted = selected.xi .- minimum_xi
    atomic_csv(joinpath(output, "spectrum.csv"), data)
    atomic_csv(joinpath(output, "qa_weights.csv"), weights)
    counting = combine(
        groupby(filter(row -> row.xi_shifted <= xi_cut, selected), :Lz2A),
        nrow => :count,
    )
    sort!(counting, :Lz2A)
    if nrow(counting) > 0
        counting.delta_L = (counting.Lz2A .- minimum(counting.Lz2A)) ./ 2
    else
        counting.delta_L = Float64[]
    end
    atomic_csv(joinpath(output, "counting.csv"), counting)
    fig = Figure(size=(720, 520))
    axis = Axis(fig[1, 1], xlabel="LzA", ylabel="xi - xi_min", title="$title, QA=$qa")
    scatter!(axis, selected.Lz2A ./ 2, selected.xi_shifted)
    hlines!(axis, [xi_cut]; linestyle=:dash)
    save(joinpath(output, "spectrum.png"), fig)
    return (QA=qa, weights=weights, counting=counting, figure=fig)
end

"""
计算 orbital entanglement spectrum (OES)。

A/B 按单粒子轨道硬切分，默认两类粒子都取各自前一半轨道。输出完整谱、
各 QA 权重、低纠缠能级 counting、图和基态摘要。
"""
function run_orbital_entanglement(
    nm1::Int,
    couplings::Couplings=Couplings(),
    settings::SolverSettings=SolverSettings(k=4);
    nm1_a::Union{Nothing,Int}=nothing,
    nm0_a::Union{Nothing,Int}=nothing,
    total_charge::Union{Nothing,Int}=nothing,
    total_lz2::Int=0,
    total_f3::Int=0,
    total_f8::Int=0,
    xi_cut::Float64=10.0,
    output::AbstractString=joinpath(PACKAGE_ROOT, "output", "oes"),
)
    output = ensure_output(output)
    write_run_metadata(output; command="oes")
    model = build_model(nm1=nm1)
    state, basis, ground, _ = _ground_state(model, couplings, settings)
    n1a = isnothing(nm1_a) ? model.nm1 ÷ 2 : nm1_a
    n0a = isnothing(nm0_a) ? model.nm0 ÷ 2 : nm0_a
    charge = isnothing(total_charge) ? 3model.nm1 : total_charge
    orbitals_a, orbitals_b, amplitudes = _orbital_cut(model; nm1_a=n1a, nm0_a=n0a)
    qnd_a = [model.qnd; GetPinOrbQNDiag(model.no, orbitals_b)]
    qnd_b = [model.qnd; GetPinOrbQNDiag(model.no, orbitals_a)]
    values1 = collect(-(model.nm1-1):2:model.nm1-1)
    values0 = collect(-(model.nm0-1):2:model.nm0-1)
    sectors_a = _diagonal_sectors(values1[1:n1a], values0[1:n0a])
    sectors_b = _diagonal_sectors(values1[n1a+1:end], values0[n0a+1:end])
    # 只保留 A/B 量子数相加恰好等于全体系量子数的 sector pair。
    sector_pairs = Vector{Vector{Int64}}[]
    for a in sectors_a
        b = (charge-a[1], total_lz2-a[2], total_f3-a[3], total_f8-a[4], 0)
        b in sectors_b && push!(sector_pairs, [Int64[a...], Int64[b...]])
    end
    sort!(sector_pairs, by=pair -> Tuple(pair[1]))
    entanglement = GetEntSpec(
        state, basis, sector_pairs, [[Int[], Int[]]];
        qnd_a=qnd_a, qnd_b=qnd_b, qnf_a=QNOffd[], amp_oa=amplitudes,
    )
    data = _entanglement_dataframe(entanglement)
    analysis = _save_entanglement(data, output, "orbital entanglement"; xi_cut=xi_cut)
    atomic_csv(joinpath(output, "summary.csv"), DataFrame([(
        nm1=nm1, mu=couplings.mu, ground_energy=ground.energy,
        trace_rho=sum(data.lambda),
        entropy=-sum(value > 0 ? value*log(value) : 0.0 for value in data.lambda),
        dominant_QA=analysis.QA,
    )]))
    return (entanglement=entanglement, data=data, ground=ground, analysis=analysis)
end

function _regularized_beta(a::Int, b::Int, x::Float64)
    value = beta_inc(a, b, x, 1-x)
    return value isa Tuple ? first(value) : value
end

_hemisphere_single(nm::Int, x::Float64) =
    [sqrt(_regularized_beta(m, nm-m+1, x)) for m in 1:nm]

function _hemisphere_amplitudes(model, x)
    # 球面单粒子轨道落在 real-space region A 的概率振幅，由正则 beta 积分给出。
    amplitudes = Float64[]
    single = _hemisphere_single(model.nm1, Float64(x))
    for value in single, _ in 1:model.nf1
        push!(amplitudes, value)
    end
    append!(amplitudes, _hemisphere_single(model.nm0, Float64(x)))
    return amplitudes
end

"""
计算 real-space entanglement spectrum (RSES)。

`x` 控制球冠面积比例（0.5 为半球）；为控制成本，只枚举中心附近的 QA 和
`|2LzA|≤lz2_cap` sector。输出格式与 OES 相同。
"""
function run_realspace_entanglement(
    nm1::Int,
    couplings::Couplings=Couplings(),
    settings::SolverSettings=SolverSettings(k=4);
    x::Float64=0.5,
    qa_half_window::Int=2,
    lz2_cap::Int=14,
    xi_cut::Float64=10.0,
    output::AbstractString=joinpath(PACKAGE_ROOT, "output", "rses"),
)
    0 < x < 1 || throw(ArgumentError("real-space cut x must lie between 0 and 1"))
    output = ensure_output(output)
    write_run_metadata(output; command="rses")
    model = build_model(nm1=nm1)
    state, basis, ground, _ = _ground_state(model, couplings, settings)
    qnd_a = QNDiag[
        pad_qn_diag(GetNeQNDiag(model.no1), 0, model.no0) +
            3pad_qn_diag(GetNeQNDiag(model.no0), model.no1, 0),
        pad_qn_diag(GetLz2QNDiag(model.nm1, model.nf1), 0, model.no0) +
            pad_qn_diag(GetLz2QNDiag(model.nm0, 1), model.no1, 0),
    ]
    total_charge = 3model.nm1
    middle = total_charge ÷ 2
    sector_pairs = Vector{Vector{Vector{Int64}}}()
    for qa in max(0, middle-qa_half_window):min(total_charge, middle+qa_half_window)
        for lz2a in -lz2_cap:2:lz2_cap
            push!(sector_pairs, [[qa, lz2a], [total_charge-qa, -lz2a]])
        end
    end
    entanglement = GetEntSpec(
        state, basis, sector_pairs, [[ComplexF64[], ComplexF64[]]];
        qnd_a=qnd_a, qnf_a=QNOffd[], amp_oa=ComplexF64.(_hemisphere_amplitudes(model, x)),
    )
    data = _entanglement_dataframe(entanglement)
    analysis = _save_entanglement(data, output, "real-space entanglement"; xi_cut=xi_cut)
    atomic_csv(joinpath(output, "summary.csv"), DataFrame([(
        nm1=nm1, mu=couplings.mu, ground_energy=ground.energy, cut_x=x,
        trace_rho=sum(data.lambda),
        entropy=-sum(value > 0 ? value*log(value) : 0.0 for value in data.lambda),
        dominant_QA=analysis.QA,
    )]))
    return (entanglement=entanglement, data=data, ground=ground, analysis=analysis)
end
