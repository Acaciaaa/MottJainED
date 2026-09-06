"""
哈密顿量的八个实系数。

`Uf, Uf0, U0` 对应三个无拉普拉斯的密度相互作用，
`Vf, Vf0, V0` 对应三个含拉普拉斯的相互作用，`t` 是三粒子转化项，
`mu` 是 charge-1 费米子的化学势。字段名与 `config/*.toml` 完全一致。
"""
Base.@kwdef struct Couplings
    Uf::Float64 = 0.5
    Uf0::Float64 = 1.8
    U0::Float64 = 4.5
    Vf::Float64 = 0.0
    Vf0::Float64 = 0.4
    V0::Float64 = 1.0
    t::Float64 = 0.5
    mu::Float64 = 0.05
end

function Couplings(values::AbstractDict)
    # 允许从 TOML.Dict 同时按 String 或 Symbol 读取参数；缺少的项用默认值。
    getf(name, default) = Float64(get(values, name, get(values, String(name), default)))
    return Couplings(
        Uf=getf(:Uf, 0.5), Uf0=getf(:Uf0, 1.8), U0=getf(:U0, 4.5),
        Vf=getf(:Vf, 0.0), Vf0=getf(:Vf0, 0.4), V0=getf(:V0, 1.0),
        t=getf(:t, 0.5), mu=getf(:mu, 0.05),
    )
end

with_coupling(c::Couplings, name::Symbol, value::Real) =
    Couplings(; (field => (field == name ? Float64(value) : getfield(c, field))
                 for field in fieldnames(Couplings))...)

coupling_vector(c::Couplings) =
    Float64[c.Uf, c.Uf0, c.U0, c.Vf, c.Vf0, c.V0, c.t, c.mu]

"""一次本征值计算所使用的数值精度、求解规模和加速选项。"""
Base.@kwdef struct SolverSettings
    k::Int = 20                    # 每个对称性 sector 需要的最低本征态数
    eig_tol::Float64 = 1.0e-8      # 稀疏本征值求解器收敛容差
    energy_tol::Float64 = 1.0e-7   # 合并不同 sector 中同能量态的容差
    quantum_tol::Float64 = 2.0e-3  # 将 <L²>, <C₂> 认成整数本征值的容差
    degeneracy_tol::Float64 = 2.0e-6 # 判断简并子空间的能量容差
    dense_cutoff::Int = 128        # 矩阵维数不超过它时直接 dense diagonalization
    ncv_extra::Int = 12            # Krylov 子空间相对 k 额外保留的向量数
    warm_start::Bool = true        # 参数扫描时是否用上一个点的基态作初始向量
end

"""一个守恒量 sector；`z`、`r` 是离散对称性量子数。"""
struct SectorKey
    z::Int
    r::Int
end

Base.show(io::IO, key::SectorKey) = print(io, "Z=$(key.z),R=$(key.r)")

"""
某个系统大小 `nm1` 的完整模型定义。

这里保存轨道数、Fock basis 规则、对称性量子数以及已经符号化构造好的
Hamiltonian/observable `Terms`。它们与具体参数值无关，因此只需建立一次。
"""
Base.@kwdef mutable struct ModelParameters
    name::Symbol = :MottJainSU3 # 模型标签，仅用于记录和输出
    nm1::Int                   # 每种 charge-1 flavor 的球面轨道数 2s+1
    s::Float64                 # charge-1 粒子的单粒子角动量 s=(nm1-1)/2
    nf1::Int                   # charge-1 flavor 数，固定为 3
    no1::Int                   # 三种 charge-1 粒子的总单粒子轨道数
    nm0::Int                   # charge-3 粒子的轨道数，3nm1-2
    nf0::Int                   # charge-3 flavor 数，固定为 1
    no0::Int                   # charge-3 的总单粒子轨道数
    no::Int                    # 两类粒子的总轨道数 no1+no0
    qnd::Any                   # diagonal conserved quantum numbers
    qnf::Any                   # off-diagonal/discrete symmetry definitions
    cfs::Any                   # 满足总电荷约束的 Fock configuration space
    hop::Any                   # f0† f1 f2 f3 的球面积分 Terms
    number_f::Any              # charge-1 总粒子数算符，乘以 mu
    v0::Any                    # V0 相互作用（保留兼容旧分析）
    l2::Any                    # 总角动量 Casimir L²
    lp::Any                    # 总角动量升算符 L+
    lm::Any                    # 总角动量降算符 L-
    c2::Any                    # SU(3) quadratic Casimir C₂
    n0::Any                    # charge-3 密度球面 observable
    nf::Any                    # 三种 charge-1 总密度 observable
    components::Any            # 八个 Hamiltonian 分量的 NamedTuple
end

"""一个 sector 的矩阵缓存；扫描参数时重复使用 basis、算符矩阵和 warm start。"""
mutable struct SectorCache
    key::SectorKey
    basis::Any
    h0::SparseMatrixCSC{Float64,Int64}       # 不含 mu 的 Hamiltonian 下三角矩阵
    number_f::SparseMatrixCSC{Float64,Int64} # 化学势算符 Nf 的下三角矩阵
    l2::Any
    c2::Any
    warm::Vector{Float64}                    # 上次求得的基态向量
end

"""一个系统大小在一组耦合参数下的全部 sector 缓存。"""
mutable struct ModelCache
    model::ModelParameters
    couplings::Couplings
    settings::SolverSettings
    sectors::Vector{SectorCache}
end

"""求解器直接得到的一个本征态及其能量、量子数、来源 sector。"""
struct SpectrumState
    energy::Float64
    l2::Float64
    c2::Float64
    sector::SectorKey
    rank::Int     # 此态在所在 sector 中按能量排序的编号
    vector::Union{Nothing,Vector{Float64}} # 仅 keep_vectors=true 时保存
    basis::Any
end

"""将能量相同且属于同一多重态的 `SpectrumState` 合并后的物理能级。"""
struct PhysicalLevel
    energy::Float64
    l2::Int
    c2::Int
    multiplicity::Int # 实际找到并合并的简并态数
    members::Vector{SpectrumState}
end

"""用低能谱与预期 CFT tower 比较后得到的目标函数及诊断信息。"""
Base.@kwdef struct CFTScore
    valid::Bool = false
    definition::Symbol = :unknown # 旧预设名；自由组合时为 custom
    terms::Vector{Symbol} = Symbol[] # 本次实际选择的 tower relations
    metric::Symbol = :q           # 当前 objective 使用 q 还是 cost
    objective::Float64 = Inf      # 当前功能实际最小化的数值
    q::Float64 = Inf       # scaling-dimension 单位的 RMS 误差
    cost::Float64 = Inf    # 旧 FSS/optimization 使用的方向夹角 cost
    factor::Float64 = NaN  # 将有限尺寸能隙归一到 CFT 标度维数的比例
    delta_s::Float64 = NaN # 提取出的 singlet 标度维数
    delta_o::Float64 = NaN # 提取出的 order-parameter 标度维数
    raw_gaps::Vector{Float64} = Float64[]
    target_gaps::Vector{Float64} = Float64[]
    labels::Vector{String} = String[]
    reason::String = ""
end

"""参数扫描中的一个独立任务点，可并行或断点续算。"""
Base.@kwdef struct ScanJob
    nm1::Int
    scan_name::String
    scan_parameter::Symbol
    scan_value::Float64
    couplings::Couplings
end

"""finite-size scaling 扫描 μ 与另一 Hamiltonian 参数时所需的设置。"""
Base.@kwdef struct FSSSettings
    nm_values::Vector{Int} = [4, 5]
    scan_parameter::Symbol = :Uf0
    scan_values::Vector{Float64} = collect(1.5:0.5:4.0)
    mu_min::Float64 = 0.0
    mu_max::Float64 = 0.12
    mu_count::Int = 9
    methods::Vector{Symbol} = [:grid, :optimize]
    score_definition::Symbol = :fss7
    score_terms::Vector{Symbol} = Symbol[]
    score_metric::Symbol = :cost
    optimize_strategy::Symbol = :size_continuation
    optimize_anchor_nm::Int = 4
    optimize_local_half_width::Float64 = 0.02
    optimize_local_count::Int = 9
    optimize_max_expansions::Int = 3
    optimize_abs_tol::Float64 = 1.0e-4
    optimize_max_iterations::Int = 60
end

# 所有允许从配置和优化器中修改的 Hamiltonian 参数名。
const HAMILTONIAN_FIELDS = (:Uf, :Uf0, :U0, :Vf, :Vf0, :V0, :t, :mu)

"""尽早检查 Hamiltonian 系数，避免 NaN/Inf 进入昂贵的矩阵构造。"""
function validate(c::Couplings)
    all(isfinite, coupling_vector(c)) || throw(ArgumentError("Hamiltonian coefficients must be finite"))
    return c
end

"""检查最基本的 eigensolver 参数合法性。"""
function validate(s::SolverSettings)
    s.k > 0 || throw(ArgumentError("solver.k must be positive"))
    s.eig_tol > 0 || throw(ArgumentError("solver.eig_tol must be positive"))
    s.energy_tol > 0 || throw(ArgumentError("solver.energy_tol must be positive"))
    s.quantum_tol > 0 || throw(ArgumentError("solver.quantum_tol must be positive"))
    return s
end
