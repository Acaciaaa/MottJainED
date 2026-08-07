# 本文件实现旧 conformal_generator.jl 的工程化版本：选出有名字的低能态，
# 构造 microscopic L=1 候选算符，用 SVD 拟合 Λ|S>≈|dS>，再检查作用结果。

"""一个被选作 CFT primary/descendant 的具体 ED 本征态。"""
Base.@kwdef mutable struct ConformalState
    label::Symbol
    name::String
    l2::Int                   # L(L+1) 的整数值，不是 L 本身
    c2::Int                   # SU(3) quadratic Casimir
    rank::Int                 # 该 (L²,C₂) 目录内第几个不同能级
    state::Vector{Float64}    # basis 中的本征向量
    basis::Any
    energy::Float64
end

"""按 `:G`, `:S`, `:dS` 等标签索引的本征态集合。"""
struct ConformalStateStore
    states::Dict{Symbol,ConformalState}
end
ConformalStateStore() = ConformalStateStore(Dict{Symbol,ConformalState}())
Base.getindex(store::ConformalStateStore, label::Symbol) = store.states[label]

"""
根据 `(l2,c2,rank)` specification 从保留向量的谱中建立有名字的态集合。

输入谱必须来自 `solve_spectrum(...; keep_vectors=true)`。
"""
function build_conformal_store(states::Vector{SpectrumState}, specifications)
    catalog, _ = level_catalog(states)
    store = ConformalStateStore()
    for specification in specifications
        label = Symbol(specification.label)
        key = (Int(specification.l2), Int(specification.c2))
        rank = Int(specification.rank)
        levels = get(catalog, key, PhysicalLevel[])
        length(levels) >= rank || error("Missing state $label at (L2,C2)=$key rank=$rank")
        member = first(levels[rank].members)
        member.vector === nothing && error("solve_spectrum must be called with keep_vectors=true")
        store.states[label] = ConformalState(
            label=label, name=String(specification.name), l2=key[1], c2=key[2],
            rank=rank, state=member.vector, basis=member.basis, energy=member.energy,
        )
    end
    return store
end

"""当前分析使用的默认 primary/descendant 标签及其有限尺寸能级位置。"""
default_conformal_specs() = [
    (label=:G, name="G", l2=0, c2=0, rank=1),
    (label=:S, name="S", l2=0, c2=0, rank=2),
    (label=:boxS, name="box S", l2=0, c2=0, rank=3),
    (label=:dS, name="dS", l2=2, c2=0, rank=1),
    (label=:T, name="T", l2=6, c2=0, rank=1),
    (label=:ddS, name="ddS", l2=6, c2=0, rank=2),
    (label=:J, name="J", l2=2, c2=3, rank=1),
    (label=:curlJ, name="curl J", l2=2, c2=3, rank=2),
    # 旧 conformal_generator.jl 在配对后的 (6,3) 物理能级中把第二项标为 dJ。
    (label=:dJ, name="dJ", l2=6, c2=3, rank=2),
]

"""
构造共形生成元的 18 个 microscopic `L=1,m=0` 候选算符。

每项同时有稳定名称和 FuzzifiED `Terms`，因此不同系统大小间可按名字比较
拟合系数，而不依赖容易变化的数组位置。
"""
function generator_candidates(model::ModelParameters)
    s = model.s
    charge1 = [GetElectronMod(model.nm1, model.nf1, flavor) for flavor in 1:3]
    charge3 = PadAngModes(GetElectronMod(model.nm0, model.nf0, 1), model.no1)

    ff_u = Dict{Tuple{Int,Int},Any}()
    ff_v = Dict{Tuple{Int,Int},Any}()
    for i in 1:3, j in i+1:3
        ff_u[(i, j)] = FilterL2(charge1[i] * charge1[j], 2s)
        ff_v[(i, j)] = FilterL2(charge1[i] * charge1[j], 2s - 2)
    end
    pair00 = charge3 * charge3
    pair00_u = FilterL2(pair00, 6s - 1)
    pair00_v = FilterL2(pair00, 6s - 3)
    mixed_u = [FilterL2(charge3 * charge1[i], 4s) for i in 1:3]
    mixed_v = [FilterL2(charge3 * charge1[i], 4s - 2) for i in 1:3]
    pair12 = FilterL2(charge1[1] * charge1[2], 2s)
    trion = FilterL2(pair12 * charge1[3], 3s)

    candidates = NamedTuple[]
    for i in 1:3, j in i+1:3
        term = ff_u[(i, j)]' * ff_u[(i, j)]
        push!(candidates, (name="Uf_$(i)$(j)", terms=GetComponent(FilterL2(term, 1), 1, 0)))
    end
    for i in 1:3, j in i+1:3
        term = ff_v[(i, j)]' * ff_v[(i, j)]
        push!(candidates, (name="Vf_$(i)$(j)", terms=GetComponent(FilterL2(term, 1), 1, 0)))
    end
    push!(candidates, (name="U0", terms=GetComponent(FilterL2(pair00_u' * pair00_u, 1), 1, 0)))
    push!(candidates, (name="V0", terms=GetComponent(FilterL2(pair00_v' * pair00_v, 1), 1, 0)))
    for i in 1:3
        push!(candidates, (name="Uf0_$i", terms=GetComponent(FilterL2(mixed_u[i]' * mixed_u[i], 1), 1, 0)))
    end
    for i in 1:3
        push!(candidates, (name="Vf0_$i", terms=GetComponent(FilterL2(mixed_v[i]' * mixed_v[i], 1), 1, 0)))
    end
    tunneling = charge3' * trion + trion' * charge3
    push!(candidates, (name="t", terms=GetComponent(FilterL2(tunneling, 1), 1, 0)))
    for i in 1:3
        density = charge1[i]' * charge1[i]
        push!(candidates, (name="mu_$i", terms=GetComponent(FilterL2(density, 1), 1, 0)))
    end
    return candidates
end

function _best_connected_target(source::ConformalState, target::ConformalState, terms_list)
    vectors = Vector{Vector{Float64}}()
    for terms in terms_list
        if isempty(terms)
            push!(vectors, zeros(Float64, length(target.state)))
        else
            operator = Operator(source.basis, target.basis, terms)
            push!(vectors, real.(operator * source.state))
        end
    end
    return vectors
end

"""
用截断 SVD 最小二乘拟合 `Λ|source> ≈ |target>`。

返回各候选系数、组合后的 `lambda_terms`、fidelity、奇异值与数值秩。
很小的奇异值按 `svd_rtol` 截掉，避免坏条件方向产生巨大系数。
"""
function fit_generator(
    source::ConformalState,
    target::ConformalState,
    candidates=generator_candidates;
    svd_rtol::Real=1.0e-10,
)
    candidate_list = candidates isa Function ? error("Pass generator_candidates(model), not the function") : candidates
    vectors = _best_connected_target(source, target, getfield.(candidate_list, :terms))
    design = hcat(vectors...)
    decomposition = svd(design; full=false)
    cutoff = Float64(svd_rtol) * maximum(decomposition.S; init=0.0)
    inverse_values = [value > cutoff ? inv(value) : 0.0 for value in decomposition.S]
    coefficients = decomposition.V * Diagonal(inverse_values) * decomposition.U' * target.state
    generated = design * coefficients
    generated_norm = real(dot(generated, generated))
    fidelity = generated_norm > 0 ? abs2(dot(target.state, generated)) / generated_norm : 0.0
    lambda_terms = SimplifyTerms(sum(
        coefficients[i] * candidate_list[i].terms for i in eachindex(candidate_list)
    ))
    return (
        coefficients=coefficients, names=String.(getfield.(candidate_list, :name)),
        terms=lambda_terms, fidelity=fidelity, singular_values=decomposition.S,
        numerical_rank=count(>(cutoff), decomposition.S), generated=generated,
    )
end

function project_angular_momentum(state, l2operator, target_l::Int, possible_l)
    # 用 Lagrange polynomial in L² 投影到 target_l，无需显式求 L² 全部本征向量。
    output = copy(state)
    target_l2 = target_l * (target_l + 1)
    for ell in possible_l
        ell == target_l && continue
        bad_l2 = ell * (ell + 1)
        output = (l2operator * output .- bad_l2 .* output) ./ (target_l2 - bad_l2)
    end
    return output
end

"""
把已拟合的 Λ 作用到一个输入态，并计算它落到指定目标态集合的权重。

可选 `target_l` 会先做角动量投影；返回逐态 overlap、总 overlap 和生成向量。
"""
function generator_overlap(
    input::ConformalState,
    targets::Vector{ConformalState},
    lambda_terms;
    target_l::Union{Nothing,Int}=nothing,
    l2_terms=nothing,
)
    isempty(targets) && throw(ArgumentError("At least one target is required"))
    basis = first(targets).basis
    all(target -> target.basis === basis, targets) ||
        throw(ArgumentError("All targets must use the same basis"))
    generated = real.(Operator(input.basis, basis, lambda_terms) * input.state)
    if target_l !== nothing
        l2_terms === nothing && throw(ArgumentError("Pass l2_terms when target_l is requested"))
        possible = input.l2 == 0 ? [1] : begin
            ell = round(Int, (sqrt(1 + 4input.l2) - 1) / 2)
            [ell - 1, ell, ell + 1]
        end
        generated = project_angular_momentum(
            generated, float_opmat(Operator(basis, l2_terms)), target_l, possible,
        )
    end
    norm2 = real(dot(generated, generated))
    overlaps = Dict(target.label => (norm2 > 0 ? abs2(dot(target.state, generated))/norm2 : 0.0) for target in targets)
    return (overlaps=overlaps, total=sum(values(overlaps)), norm2=norm2, generated=generated)
end

"""
执行默认共形生成元分析：求带向量的谱、选态、造候选、拟合 S→dS 并保存。

输出系数 CSV、摘要 CSV 和含完整拟合/态对象的 `generator.jld2`。
"""
function run_generator_analysis(
    nm1::Int,
    couplings::Couplings=Couplings(),
    settings::SolverSettings=SolverSettings(k=80);
    output::AbstractString=joinpath(PACKAGE_ROOT, "output", "generator"),
)
    output = ensure_output(output)
    write_run_metadata(output; command="generator")
    model = build_model(nm1=nm1)
    cache = prepare_spectrum(model, couplings, settings)
    states = solve_spectrum(cache, couplings.mu; keep_vectors=true)
    store = build_conformal_store(states, default_conformal_specs())
    candidates = generator_candidates(model)
    fit = fit_generator(store[:S], store[:dS], candidates)
    coefficient_table = DataFrame(name=fit.names, coefficient=fit.coefficients)
    atomic_csv(joinpath(output, "generator_coefficients.csv"), coefficient_table)
    jldsave(joinpath(output, "generator.jld2"); fit=fit, store=store, couplings=couplings)
    summary = DataFrame([(
        nm1=nm1, mu=couplings.mu, fidelity=fit.fidelity,
        numerical_rank=fit.numerical_rank, candidate_count=length(candidates),
    )])
    atomic_csv(joinpath(output, "generator_summary.csv"), summary)
    return (fit=fit, store=store, model=model, cache=cache)
end
