# SO(3)lver 实验分支

本文件只适用于 `so3lver-experiment` 分支；SO(3)lver 功能没有合并到 `main`
或 `fast-ed-experiment`。三个分支的依赖声明均已统一到 FuzzifiED 2.0.1；
已有 N=7 稀疏矩阵缓存和既有物理结论未修改，也不得因依赖升级而删除。

## 固定的软件版本

- MottJainED 分支：`so3lver-experiment`
- FuzzifiED：2.0.1
- FuzzifiED Git commit：`877d1403f97a7044645af48df5caa8cd0e211a8a`
- Project source：上述 GitHub commit

开发和验证时可以采用两个隔离工作树：

```text
MottJainED-so3lver/   # 本实验分支
FuzzifiED-so3lver/    # 上述固定 commit
```

服务器共享源码目录 `/public/home/ruiqixu/FuzzifiED.jl` 也更新到同一 commit；
MottJainED 的各 Julia 环境则按 Project 中的 Git source 解析，不再混用 path override：

```bash
git -C /public/home/ruiqixu/FuzzifiED.jl fetch origin
git -C /public/home/ruiqixu/FuzzifiED.jl checkout --detach 877d1403f97a7044645af48df5caa8cd0e211a8a
julia --project=. -e 'using Pkg; Pkg.add(PackageSpec(url="https://github.com/FuzzifiED/FuzzifiED.jl.git", rev="877d1403f97a7044645af48df5caa8cd0e211a8a")); Pkg.instantiate()'
```

依赖升级会改变后续新缓存的环境身份；保留旧缓存和归档，不要跨依赖版本强行复用。

## 实现内容

`experimental/SO3lverED.jl` 将 Hilbert space 分成两个 segment：

1. 三个 charge-1 SU(3) flavor；
2. 一个 charge-3 SU(3) singlet。

八个 Hamiltonian 分量与 `src/Model.jl` 完全对应：

```text
Uf, Uf0, U0, Vf, Vf0, V0, t, mu
```

同一 segment 的相互作用使用 `SingleSegCouple`，交叉密度和三粒子转化使用
`ContactCouple`。`Vf0` 中由 Laplacian 产生的恒等零 `l=0` channel 在进入底层
`Operator` 前被严格删除；否则 FuzzifiED 2.0.1 会对空 `Terms` 做 reduction。

支持的精确物理 block 为：

| representation | SU(3) 选择 | C2 | L |
|---|---:|---:|---:|
| singlet | `(F3,F8)=(0,0)` | 0 | 0, 1, 2, ... |
| adjoint | 最高权 `(F3,F8)=(1,3)` | 3 | 0, 1, 2, ... |

adjoint 使用最高权而不是旧算法的 Cartan zero-weight sector，因此每个 octet
只出现一次，不再带两个 zero-weight 副本。旧流程中的 raw rank 不能直接照搬；
例如 curl-J 应按物理 adjoint block 内的 rank 解释。

`SO3Workspace` 保存昂贵的 segment spaces 和 reduced matrix elements；
`retune!` 只更新各 channel 的数值系数，不重建 workspace。这是后续参数优化必须
复用的边界。

## 正确性

运行：

```bash
julia --project=. --check-bounds=yes test/so3lver_runtests.jl
```

当前专项结果为 79/79 tests passed。测试在 N=2 上逐 block 比较传统 Fock-basis ED
与 SO(3)lver，覆盖：

- singlet 和 adjoint；
- L=0、1、2；
- 一组非规则系数下的全部八个 Hamiltonian 分量；
- 每个分量未经 `Symmetric` 包装的显式 Hermiticity；
- 强制走 matrix-free Krylov 路径的低能谱；
- workspace 不变、只改系数的 `retune!` 路径。

所有能量在 `2e-11` 绝对误差内一致。

## 本地性能结果

命令：

```bash
julia -t auto --project=. scripts/so3lver_ed.jl benchmark \
  --nm=7 --representation=adjoint --ell=2 --k=4
```

机器：Apple Silicon，24 GiB RAM；Julia 实际使用 4 threads。以下均为 adjoint
L=2，默认 Hamiltonian，单个新 Julia 进程：

| N | exact block dim | workspace | operator | steady matvec | lowest 4 |
|---:|---:|---:|---:|---:|---:|
| 3 | 6 | 7.09 s | 1.63 s | 未单列 | 3.03 s |
| 4 | 41 | 6.65 s | 1.31 s | 0.0010 s | 2.68 s |
| 5 | 360 | 6.79 s | 1.38 s | 0.0059 s | 3.63 s |
| 6 | 3,747 | 7.82 s | 1.58 s | 0.0267 s | 5.96 s |
| 7 | 45,085 | 14.78 s | 1.71 s | 0.0798 s | 41.16 s |

这张表的 N=3...7 求解时间采自修复费米奇宇称转化项相对符号之前，只用于资源
和性能估算；其中 block 维数、workspace/operator/matvec 时间仍有效，旧输出中的
能量不得用于物理分析。符号修复后的能量应在服务器重新生成。

原 N=7 方法的最大离散对称 sector 约 2,932,946 维；SO(3)lver 的目标
adjoint L=2 block 是 45,085 维。这里只比较 Hilbert-space reduction，两个实现的
算符表示不同，不能用维数比直接当作端到端加速比。

N=8 本地压力测试按用户要求提前停止，以免继续加热本机，不应记录为完成的
benchmark。停止前已确认：

- 两个 retained segment dimensions：9,240 和 16,425；
- reduced-operator 构造期间 RSS 约 3.7 GiB；
- 进程正常且未发生 OOM；
- 独立精确计数给出的 adjoint L=2 block dimension 为 606,185。

本地生成的完整 N=3...7 TOML 记录位于
`output/so3lver/benchmarks/`（该目录按项目规则不纳入 Git）。

## 服务器资源建议

单区块压力测试可以使用：

```bash
julia -t 8 --project=. scripts/so3lver_ed.jl benchmark \
  --nm=8 --representation=adjoint --ell=2 --k=4
```

实际调参性能必须使用六区块 workload，而不能把上述单区块时间当成一次 CFT score：

```bash
julia -t 8 --project=. scripts/so3lver_cft_workload.jl \
  --nm=7 --k=4 --mus=0.1216874050164 \
  --Uf=0.46 --Uf0=1.834 --U0=4.14 --Vf=0 \
  --Vf0=0.41 --V0=0.525 --t=0.5
```

该入口在同一个进程中各构造一次 singlet/adjoint light workspace，并让两种表示共享完全
相同的 charge-3 heavy segment，再复用于六个 `(representation,L)` block。这一共享是精确
重用，不截断 Hilbert space。`--mus` 可以给逗号分隔的多个化学势；从第二点开始每个
block 默认使用上一点基态 warm start。输出直接包含固定五条
`dS-S,J,curlJ,dJ(rank1),T(rank1)` 的 `q/factor/DeltaS/DeltaO`。adjoint 最高权 block
每个 octet 只出现一次，所以旧 raw rank 3 的 curl-J 在这里是物理 rank 2。

正式参数扫描必须在同一 Julia 进程内保留 workspace，依次 `retune!`；不能对每个参数
点重新执行单区块 benchmark 脚本。与旧 FastED 不同，SO(3)lver 的 reduced operators
允许复用全部八个 Hamiltonian 系数，不只复用 `mu`。

六个当前需要的 block 维数为：

| N | singlet L=0 | singlet L=1 | singlet L=2 | adjoint L=0 | adjoint L=1 | adjoint L=2 |
|---:|---:|---:|---:|---:|---:|---:|
| 7 | 4,327 | 12,290 | 20,585 | 9,189 | 27,404 | 45,085 |
| 8 | 49,920 | 146,802 | 244,016 | 122,826 | 366,925 | 606,185 |
| 9 | 655,293 | 1,951,649 | 3,238,238 | 1,798,342 | 5,378,406 | 8,908,622 |

不要直接提交 N=9。当前 SO(3)lver 的 segment 构造仍会 dense diagonalize 最大约
75k--111k 的 representative block；即使最终 composite Hamiltonian 是 matrix-free，
这个前处理也会成为内存和 O(n^3) 时间瓶颈。N=9 需要先实现 segment workspace
持久缓存或改进 segment 对角化算法，再做单 block pilot。

## 尚未覆盖

- orbital/real-space entanglement 仍需要 Fock-basis coefficients，不能直接复用
  SO(3)lver 的 coupled basis；
- workspace 尚未持久化到磁盘；
- 六区块 workload 已能计算现行五项 CFT score；当前 N=6 调参使用去掉不稳定
  `boxS` 后的六项 score，并按 `docs/SO3LVER_PARAMETER_SEARCH.md` 对每个外层点
  完整重做 `mu` profile；旧七项 generator audit 只保留作历史诊断；
- N=8 六区块、两个 mu 的服务器实测已完成：峰值约 14.6 GiB，workspace 约 69 分钟，
  workspace 建成后的单个完整 CFT 点约 4--5 分钟；
- N=6 的 `Uf0,Vf0,V0` 联合搜索已由
  `slurm/so3lver_n6_nested_optimize.sbatch` 提交；`ddS` 保留固定 rank/overlap 门，
  `boxS` 不进入目标函数或身份门。
