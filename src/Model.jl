# 本文件只负责“物理模型 -> FuzzifiED Terms / 稀疏矩阵”的转换。
# 它不求本征值，也不做参数扫描。

# FuzzifiED 原始量子数/算符只覆盖一类粒子的轨道；pad_* 把它们嵌入
# charge-1 与 charge-3 共同组成的总 Hilbert space。
pad_qn_diag(qn::QNDiag, left::Int, right::Int) = QNDiag(
    qn.name,
    [fill(0, left); qn.charge; fill(0, right)],
    qn.modul,
)

pad_qn_offdiag(qn::QNOffd, left::Int, right::Int) = QNOffd(
    [collect(1:left); qn.perm .+ left; collect(1:right) .+ left .+ length(qn.perm)],
    [fill(0, left); qn.ph; fill(0, right)],
    [fill(ComplexF64(1), left); qn.fac; fill(ComplexF64(1), right)],
    qn.cyc,
)

pad_term(term::Term, left::Int) = Term(
    term.coeff,
    [isodd(i) ? term.cstr[i] : term.cstr[i] + left for i in eachindex(term.cstr)],
)
pad_term(terms::Terms, left::Int) = pad_term.(terms, left)
pad_sphere_observable(obs::SphereObs, left::Int) =
    SphereObs(obs.s2, obs.l2m, (l, m) -> pad_term(obs.get_comp(l, m), left))

"""
Build the charge-1 SU(3) plus charge-3 fuzzy-sphere model.

`nm1` is the number of charge-1 orbitals.  All radius-dependent terms are
constructed with an explicit `norm_r2=nm1`, avoiding accidental dependence on
FuzzifiED's global radius when several sizes are handled in one process.
"""
function build_model(; nm1::Int)
    nm1 >= 2 || throw(ArgumentError("nm1 must be at least 2"))
    # charge-1: 三个 SU(3) flavor，各有 nm1=2s+1 个轨道。
    s = (nm1 - 1) / 2
    nf1 = 3
    no1 = nm1 * nf1
    # charge-3 粒子看到三倍磁通，因此有 6s+1 = 3nm1-2 个轨道。
    nm0 = 3nm1 - 2
    nf0 = 1
    no0 = nm0
    no = no1 + no0
    radius2 = Float64(nm1)
    FuzzifiED.ObsNormRadSq = radius2

    # qnd 定义 diagonal 守恒量：总电荷、2Lz、两个 SU(3) Cartan charge。
    qnd = QNDiag[
        pad_qn_diag(GetNeQNDiag(no1), 0, no0) + 3pad_qn_diag(GetNeQNDiag(no0), no1, 0),
        pad_qn_diag(GetLz2QNDiag(nm1, nf1), 0, no0) + pad_qn_diag(GetLz2QNDiag(nm0, 1), no1, 0),
        pad_qn_diag(GetFlavQNDiag(nm1, nf1, [1, -1, 0]), 0, no0),
        pad_qn_diag(GetFlavQNDiag(nm1, nf1, [1, 1, -2]), 0, no0),
    ]
    # qnf 定义 flavor permutation 与球面 y-rotation 离散对称性。
    qnf = QNOffd[
        pad_qn_offdiag(GetFlavPermQNOffd(nm1, nf1, [2, 1, 3], [1, -1, 1]), 0, no0),
        pad_qn_offdiag(GetRotyQNOffd(nm1, nf1), 0, no0) *
            pad_qn_offdiag(GetRotyQNOffd(nm0, 1), no1, 0),
    ]
    # 固定总 charge = 3*nm1；[no1,0,0,0] 是这些 QN 的目标值。
    cfs = Dict(0 => Confs(no, [no1, 0, 0, 0], qnd))

    f0 = pad_sphere_observable(GetElectronObs(nm0, 1, 1), no1)
    f = [GetElectronObs(nm1, nf1, flavor) for flavor in 1:nf1]
    n0 = f0' * f0
    nf = GetDensityObs(nm1, nf1)

    # 局域 transition：一个 charge-3 粒子与三个不同 flavor 的 charge-1
    # 粒子相互转化。Hermitian conjugate 在 components.t 中补上。
    hop = SimplifyTerms(GetIntegral(f0' * f[1] * f[2] * f[3]; norm_r2=radius2))
    number_f = GetPolTerms(nm1, nf1)

    # 显式构造总 Lz、L+、L-，随后用 L²=Lz²-Lz+L+L- 得到 Casimir。
    lz1 = [
        begin
            m = div(o - 1, nf1)
            Term(m - s, [1, o, 0, o])
        end for o in 1:no1
    ]
    lp1 = [
        begin
            m = div(o - 1, nf1)
            Term(sqrt(m * (nm1 - m)), [1, o, 0, o - nf1])
        end for o in nf1 + 1:no1
    ]
    lz0 = [Term(m - 3s - 1, [1, m + no1, 0, m + no1]) for m in 1:nm0]
    lp0 = [Term(sqrt((m - 1) * (nm0 - m + 1)), [1, m + no1, 0, m - 1 + no1]) for m in 2:nm0]
    lz = lz1 + lz0
    lp = lp1 + lp0
    lm = lp'
    l2 = SimplifyTerms(lz * lz - lz + lp * lm)
    c2 = GetC2Terms(nm1, nf1, :SU)

    # 预先构造八个与数值系数无关的 Hamiltonian 分量。扫描参数时只需
    # 线性组合这些 Terms，无需反复进行球面积分和符号化简。
    components = (
        Uf=SimplifyTerms(GetIntegral(nf * nf; norm_r2=radius2)),
        Uf0=SimplifyTerms(GetIntegral(nf * n0; norm_r2=radius2)),
        U0=SimplifyTerms(GetIntegral(n0 * n0; norm_r2=radius2)),
        Vf=SimplifyTerms(GetIntegral(nf * Laplacian(nf; norm_r2=radius2); norm_r2=radius2)),
        Vf0=SimplifyTerms(GetIntegral(n0 * Laplacian(nf; norm_r2=radius2); norm_r2=radius2)),
        V0=SimplifyTerms(GetIntegral(n0 * Laplacian(n0; norm_r2=radius2); norm_r2=radius2)),
        t=SimplifyTerms(-(hop + hop')),
        mu=number_f,
    )

    return ModelParameters(
        nm1=nm1, s=Float64(s), nf1=nf1, no1=no1, nm0=nm0, nf0=nf0,
        no0=no0, no=no, qnd=qnd, qnf=qnf, cfs=cfs, hop=hop,
        number_f=number_f, v0=components.V0, l2=l2, lp=lp, lm=lm,
        c2=c2, n0=n0, nf=nf, components=components,
    )
end

"""Return the Hamiltonian as FuzzifiED `Terms`."""
function hamiltonian_terms(model::ModelParameters, c::Couplings; include_mu::Bool=true)
    validate(c)
    p = model.components
    terms = c.Uf*p.Uf + c.Uf0*p.Uf0 + c.U0*p.U0 +
            c.Vf*p.Vf + c.Vf0*p.Vf0 + c.V0*p.V0 + c.t*p.t
    include_mu && (terms += c.mu*p.mu)
    return SimplifyTerms(terms)
end

function lower_sparse(mat::OpMat{Float64})
    # FuzzifiED 的存储形似 CSC，但同一列内 row index 不保证有序。
    # Julia 稀疏运算要求 canonical CSC；直接包裹底层数组可能静默丢元素。
    # 通过 sparse(rows, cols, vals, ...) 重新排序并合并重复位置。
    columns = Vector{Int64}(undef, mat.nel)
    for column in 1:mat.dimd
        columns[mat.colptr[column]:mat.colptr[column+1]-1] .= column
    end
    return sparse(
        copy(mat.rowid), columns, copy(mat.elval), mat.dimf, mat.dimd,
    )
end

"""把 FuzzifiED `Operator` 真正作用到 basis 上，得到 Float64 `OpMat`。"""
float_opmat(operator::Operator) = OpMat(operator; type=Float64)

"""把只存一半的 Hermitian 稀疏矩阵重新标记为 FuzzifiED Hermitian OpMat。"""
function hermitian_opmat(lower::SparseMatrixCSC{Float64,Int64})
    mat = OpMat(lower)
    mat.sym_q = 1
    return mat
end
