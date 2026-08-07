"""
为固定 `nm1` 和非 μ 耦合参数准备求谱缓存。

最昂贵的 basis、H₀、Nf、L²、C₂ 矩阵在这里生成一次。以后改变 μ 时只做
`H(μ)=H₀+μNf`，这是 FSS/临界点搜索能复用计算的关键。
"""
function prepare_spectrum(
    model::ModelParameters,
    couplings::Couplings,
    settings::SolverSettings=SolverSettings(),
)
    validate(couplings)
    validate(settings)
    h0terms = hamiltonian_terms(model, with_coupling(couplings, :mu, 0.0); include_mu=false)
    sectors = SectorCache[]

    # 四个 (Z,R) 离散对称 sector 分开对角化，既减小矩阵，也保留态的标签。
    for z in (1, -1), r in (1, -1)
        key = SectorKey(z, r)
        basis = Basis(model.cfs[0], [z, r], model.qnf)
        basis.dim == 0 && continue
        h0 = lower_sparse(float_opmat(Operator(basis, h0terms)))
        number_f = lower_sparse(float_opmat(Operator(basis, model.number_f)))
        l2 = float_opmat(Operator(basis, model.l2))
        c2 = float_opmat(Operator(basis, model.c2))
        push!(sectors, SectorCache(key, basis, h0, number_f, l2, c2, Float64[]))
    end
    isempty(sectors) && error("No non-empty symmetry sectors were found")
    return ModelCache(model, couplings, settings, sectors)
end

"""
修改非 μ 参数后，仅重建 H₀；basis、Nf、L²、C₂ 继续复用。

参数优化会反复调用它；`reset_warm=true` 可同时丢弃上一个参数点的初始向量。
"""
function retune_spectrum!(cache::ModelCache, couplings::Couplings; reset_warm::Bool=false)
    validate(couplings)
    h0terms = hamiltonian_terms(cache.model, with_coupling(couplings, :mu, 0.0); include_mu=false)
    for sector in cache.sectors
        sector.h0 = lower_sparse(float_opmat(Operator(sector.basis, h0terms)))
        reset_warm && empty!(sector.warm)
    end
    cache.couplings = couplings
    return cache
end

function _projected_matrix(op, vectors::Matrix{Float64})
    # 计算算符在一组近简并本征向量张成的子空间中的小矩阵 V' O V。
    count = size(vectors, 2)
    result = zeros(Float64, count, count)
    applied = [op * vectors[:, j] for j in 1:count]
    for j in 1:count, i in 1:j
        value = real(dot(vectors[:, i], applied[j]))
        result[i, j] = value
        result[j, i] = value
    end
    return result
end

function _clusters(values::AbstractVector{<:Real}, atol::Real)
    # 将已排序的数值按相对容差切成连续的近简并区块。
    isempty(values) && return UnitRange{Int}[]
    groups = UnitRange{Int}[]
    first_index = 1
    for i in 2:length(values)
        scale = max(1.0, abs(values[first_index]), abs(values[i]))
        if abs(values[i] - values[first_index]) > atol * scale
            push!(groups, first_index:i-1)
            first_index = i
        end
    end
    push!(groups, first_index:length(values))
    return groups
end

"""
Resolve accidental numerical mixing in degenerate energy subspaces.

ARPACK may return arbitrary linear combinations when states are degenerate.
This routine diagonalizes L² and then C₂ inside each energy block before states
are classified.  It prevents a valid multiplet from being lost merely because
an expectation value lies between two exact quantum numbers.
"""
function _resolve_quantum_numbers!(
    energies::Vector{Float64},
    vectors::Matrix{Float64},
    l2op,
    c2op,
    settings::SolverSettings,
)
    for block in _clusters(energies, settings.energy_tol)
        idx = collect(block)
        length(idx) == 1 && continue
        old_vectors = copy(vectors[:, idx])
        rotation = eigen(Symmetric(_projected_matrix(l2op, old_vectors))).vectors
        vectors[:, idx] .= old_vectors * rotation
        old_energies = copy(energies[idx])
        energies[idx] .= [sum(abs2.(rotation[:, j]) .* old_energies) for j in axes(rotation, 2)]

        l2values = [real(dot(vectors[:, j], l2op * vectors[:, j])) for j in idx]
        order = sortperm(l2values)
        vectors[:, idx] .= vectors[:, idx[order]]
        energies[idx] .= energies[idx[order]]
        l2values = l2values[order]

        for subgroup in _clusters(l2values, settings.quantum_tol)
            subidx = idx[collect(subgroup)]
            length(subidx) == 1 && continue
            old_subvectors = copy(vectors[:, subidx])
            crotation = eigen(Symmetric(_projected_matrix(c2op, old_subvectors))).vectors
            vectors[:, subidx] .= old_subvectors * crotation
            old_subenergies = copy(energies[subidx])
            energies[subidx] .= [sum(abs2.(crotation[:, j]) .* old_subenergies) for j in axes(crotation, 2)]
        end
    end
    return energies, vectors
end

function _eigensystem(sector::SectorCache, mu::Float64, settings::SolverSettings)
    # H0 与 Nf 都只保存 Hermitian 矩阵的一半；这里才组合当前 μ 的矩阵。
    lower = sector.h0 + mu * sector.number_f
    dropzeros!(lower)
    hamiltonian = hermitian_opmat(lower)
    n = hamiltonian.dimd

    if n == 1
        return Float64[Matrix(hamiltonian)[1, 1]], ones(Float64, 1, 1)
    elseif n <= settings.dense_cutoff
        # 小 sector 全对角化更稳定；大 sector 才使用 FuzzifiED/ARPACK 稀疏求解。
        decomp = eigen(Symmetric(Matrix(hamiltonian)))
        count = min(settings.k, n)
        order = sortperm(decomp.values)[1:count]
        return Float64.(decomp.values[order]), Matrix{Float64}(decomp.vectors[:, order])
    end

    count = min(settings.k, n - 2)
    count > 0 || error("Sector dimension $n is too small for sparse diagonalization")
    ncv = min(n - 1, max(2count, count + settings.ncv_extra))
    ncv > count || (ncv = min(n - 1, count + 1))
    kwargs = (; tol=settings.eig_tol, ncv=ncv, disp_std=false)
    if settings.warm_start && length(sector.warm) == n
        # 连续参数点的基态通常很接近，用前一点向量可显著减少 Krylov 迭代。
        energies, vectors = GetEigensystem(hamiltonian, count; initvec=sector.warm, kwargs...)
    else
        energies, vectors = GetEigensystem(hamiltonian, count; kwargs...)
    end
    return Float64.(real.(energies)), Matrix{Float64}(real.(vectors))
end

"""
在指定化学势 `mu` 下对角化所有 `(Z,R)` sector，返回按能量排序的态。

默认不保留大本征向量以节省内存；密度、共形生成元、纠缠谱需要传入
`keep_vectors=true`。
"""
function solve_spectrum(cache::ModelCache, mu::Real; keep_vectors::Bool=false)
    isfinite(mu) || throw(ArgumentError("mu must be finite"))
    all_states = SpectrumState[]
    for sector in cache.sectors
        energies, vectors = _eigensystem(sector, Float64(mu), cache.settings)
        order = sortperm(energies)
        energies = energies[order]
        vectors = vectors[:, order]
        # ARPACK 对简并态会返回任意线性组合；重新在简并子空间对角化 L²/C₂。
        _resolve_quantum_numbers!(energies, vectors, sector.l2, sector.c2, cache.settings)
        order = sortperm(energies)
        energies = energies[order]
        vectors = vectors[:, order]
        cache.settings.warm_start && (sector.warm = copy(vectors[:, 1]))

        for rank in eachindex(energies)
            vector = vectors[:, rank]
            l2 = real(dot(vector, sector.l2 * vector))
            c2 = real(dot(vector, sector.c2 * vector))
            stored_vector = keep_vectors ? copy(vector) : nothing
            push!(all_states, SpectrumState(
                energies[rank], l2, c2, sector.key, rank, stored_vector,
                keep_vectors ? sector.basis : nothing,
            ))
        end
    end
    sort!(all_states, by=state -> state.energy)
    return all_states
end

"""
只求四个 `(Z,R)` sector 中的全局基态，并保留基态向量。

密度只需要 `⟨ψ₀|Nf|ψ₀⟩`，因此这条路径不计算各个低能态的 `L²/C₂`，
也不在简并子空间中做额外的量子数分类。
"""
function solve_ground_state(cache::ModelCache, mu::Real)
    isfinite(mu) || throw(ArgumentError("mu must be finite"))
    best = nothing
    for sector in cache.sectors
        energies, vectors = _eigensystem(sector, Float64(mu), cache.settings)
        index = argmin(energies)
        vector = copy(vectors[:, index])
        cache.settings.warm_start && (sector.warm = copy(vector))
        candidate = SpectrumState(
            energies[index], NaN, NaN, sector.key, 1, vector, sector.basis,
        )
        if isnothing(best) || candidate.energy < best.energy
            best = candidate
        end
    end
    isnothing(best) && error("No ground-state candidate was found")
    return best
end

"""把数值 Casimir 在容差内识别成整数标签；否则返回 `nothing`。"""
function _integer_quantum_number(value::Real, tolerance::Real)
    rounded = round(Int, value)
    scale = max(1.0, abs(value), abs(rounded))
    return abs(value - rounded) <= tolerance * scale ? rounded : nothing
end

"""统一取得一个本征态的整数 `(L²,C₂)` 标签；无法可靠识别时返回 `nothing`。"""
function _state_quantum_labels(state::SpectrumState, quantum_tol::Real)
    l2 = _integer_quantum_number(state.l2, quantum_tol)
    c2 = _integer_quantum_number(state.c2, quantum_tol)
    return isnothing(l2) || isnothing(c2) ? nothing : (l2, c2)
end

"""按统一的容差整数标签选出一个 `(L²,C₂)` sector，保留原始能量顺序和副本。"""
function _states_in_sector(
    states::Vector{SpectrumState}, l2::Integer, c2::Integer;
    quantum_tol::Real=2.0e-3,
)
    selected = filter(states) do state
        _state_quantum_labels(state, quantum_tol) == (Int(l2), Int(c2))
    end
    return sort(selected; by=state -> state.energy)
end

"""
把不同 `(Z,R)` sector 中同能、同 `(L²,C₂)` 的副本合并为物理能级目录。

返回 `(catalog, rejected)`：`catalog[(L²,C₂)]` 是按能量排列的能级；
`rejected` 收集无法在容差内识别为整数 Casimir 的态，便于诊断数值问题。
"""
function level_catalog(
    states::Vector{SpectrumState};
    quantum_tol::Real=2.0e-3,
    degeneracy_tol::Real=2.0e-6,
)
    grouped = Dict{Tuple{Int,Int},Vector{SpectrumState}}()
    rejected = SpectrumState[]
    for state in states
        labels = _state_quantum_labels(state, quantum_tol)
        if isnothing(labels)
            push!(rejected, state)
        else
            push!(get!(grouped, labels, SpectrumState[]), state)
        end
    end

    catalog = Dict{Tuple{Int,Int},Vector{PhysicalLevel}}()
    for (quantum_numbers, members) in grouped
        sort!(members, by=state -> state.energy)
        levels = PhysicalLevel[]
        start = 1
        for i in 2:length(members) + 1
            at_end = i > length(members)
            if !at_end
                scale = max(1.0, abs(members[start].energy), abs(members[i].energy))
                same = abs(members[i].energy - members[start].energy) <= degeneracy_tol * scale
                same && continue
            end
            block = members[start:i-1]
            energy = mean(state.energy for state in block)
            push!(levels, PhysicalLevel(
                energy, quantum_numbers[1], quantum_numbers[2], length(block), copy(block),
            ))
            start = i
        end
        catalog[quantum_numbers] = levels
    end
    return catalog, rejected
end

function spectrum_dataframe(
    states::Vector{SpectrumState}; mu::Real, nm1::Int, quantum_tol::Real=2.0e-3,
)
    # l2/c2 保存可靠的整数标签；l2_raw/c2_raw 同时保留原始期望值用于数值诊断。
    labels = [_state_quantum_labels(state, quantum_tol) for state in states]
    return DataFrame(
        nm1=fill(nm1, length(states)), mu=fill(Float64(mu), length(states)),
        energy=getfield.(states, :energy),
        l2=[isnothing(label) ? missing : label[1] for label in labels],
        c2=[isnothing(label) ? missing : label[2] for label in labels],
        l2_raw=getfield.(states, :l2), c2_raw=getfield.(states, :c2),
        z=[state.sector.z for state in states],
        r=[state.sector.r for state in states], rank=getfield.(states, :rank),
    )
end
