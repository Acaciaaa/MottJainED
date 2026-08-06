# MottJainED 使用手册

## 1. 项目目标

`MottJainED` 是旧 `mott_jain` 工作流的独立重构。它保留当前哈密顿量

\[
H=U_f\!\int n_f^2+U_{f0}\!\int n_fn_0+U_0\!\int n_0^2
 +V_f\!\int n_f\nabla^2n_f+V_{f0}\!\int n_0\nabla^2n_f
 +V_0\!\int n_0\nabla^2n_0-t(f_0^\dagger f_1f_2f_3+h.c.)+\mu N_f,
\]

但把模型、数值求解、物理诊断、数据存取和画图分离。项目不会修改
`../mott_jain` 或 `../FuzzifiED.jl`。

| 旧文件 | 新入口 |
|---|---|
| `pad_su3.jl` | `src/Model.jl`、`src/Spectrum.jl` |
| `ES_mu.jl` | `spectrum`、`gap` |
| `find_critical_point.jl` | `critical` |
| `optimization.jl` | `optimize` |
| `FSS1.jl`、`FSS2.jl` | `fss`、`fss-plot`、`fss-fit` |
| `particle_density_mu.jl` | `density` |
| `scaling_dimension.jl` | `scaling` |
| `conformal_generator.jl` | `generator`、`src/Conformal.jl` |
| `orbital_entangle.jl` | `oes` |
| `realspace_entangle.jl` | `rses` |

## 2. 安装与快速检查

```bash
cd /Users/ruiqi/Documents/hkust/research/fuzzysphere/MottJainED
julia --project=. scripts/setup.jl
julia --project=. bin/mottjain.jl plan
julia --project=. -e 'using Pkg; Pkg.test()'
```

项目通过相对路径绑定本地 `../FuzzifiED.jl`。复制到服务器时，最好保持
`MottJainED` 与 `FuzzifiED.jl` 为同级目录；若布局不同，执行：

```julia
using Pkg
Pkg.activate("/path/to/MottJainED")
Pkg.develop(path="/actual/path/to/FuzzifiED.jl")
Pkg.instantiate()
```

## 3. 用配置改变任务

先复制默认配置：

```bash
cp config/default.toml config/my_run.toml
```

最常改的区域是：

```toml
[model]
nm1 = 6
nm_values = [4, 5, 6, 7]

[hamiltonian]
Uf = 0.5
Uf0 = 1.8
U0 = 4.5
Vf = 0.0
Vf0 = 0.4
V0 = 1.0
t = 0.5
mu = 0.05

[solver]
k = 30
fuzzified_threads = 16
blas_threads = 1
```

服务器上建议令 Julia 线程数与申请的 CPU 数一致，并保持 BLAS 为 1，避免
FuzzifiED、Julia 和 BLAS 三层线程互相超额占用。

## 4. 命令

统一调用形式：

```bash
julia --startup-file=no --project=. bin/mottjain.jl COMMAND --config=config/my_run.toml
```

### 谱、singlet gap 与密度

```bash
julia --project=. bin/mottjain.jl spectrum --config=config/my_run.toml
julia --project=. bin/mottjain.jl gap     --config=config/my_run.toml
julia --project=. bin/mottjain.jl density --config=config/my_run.toml
```

`spectrum` 为每个 μ 写一个独立 CSV，最后才登记到 `summary.csv`。中断后重跑
会根据稳定 `job_id` 跳过已完成点；加入 `--force` 可覆盖重算。

### 临界 μ

```bash
julia --project=. bin/mottjain.jl critical --config=config/my_run.toml
```

算法先在 `[mu_min,mu_max]` 做 coarse grid，再在最优网格点相邻区间进行 Brent
优化。这比直接在整个区间做一次局部搜索更不容易落入坏点。

### 多参数优化

```toml
[optimization]
free = ["Uf", "Uf0", "Vf0", "V0", "mu"]
max_iterations = 200

[optimization.bounds]
Uf = [0.0, 10.0]
Uf0 = [0.0, 20.0]
Vf0 = [-20.0, 20.0]
V0 = [-20.0, 20.0]
mu = [-1.0, 1.0]
```

```bash
julia --project=. bin/mottjain.jl optimize --config=config/my_run.toml
```

程序只对自由方向各构造一次稀疏矩阵。每次 objective evaluation 会立即追加到
`evaluations.csv`，包括无效点和错误原因，不再静默吞掉异常。中断后重跑会从
已有 trace 中最好的有效参数重新启动（优化器的整个 simplex 不会序列化）。

### Finite-size scaling

```toml
[model]
nm_values = [4, 5, 6, 7]

[fss]
scan_parameter = "Uf0"
scan_values = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
mu_min = 0.0
mu_max = 0.12
coarse_points = 9
```

```bash
julia --project=. bin/mottjain.jl fss      --config=config/my_run.toml
julia --project=. bin/mottjain.jl fss-plot --config=config/my_run.toml --y=delta_s
julia --project=. bin/mottjain.jl fss-fit  --config=config/my_run.toml --y=delta_s
```

也可用 `fss-all` 顺序完成三步。联合拟合采用

\[
\Delta(N,g)=\Delta_\infty+a_gN^{-\omega/2}
              =\Delta_\infty+a_gx^\omega,\qquad x=N^{-1/2},
\]

所有扫描参数共享 `Delta_inf` 与 `omega`，每条曲线有独立振幅 `a_g`。至少需要
三个系统大小，且数据点必须多于参数数目；条件不足时会明确报错，不会生成一个
看似正常但实际上欠定的拟合。

### Scaling dimension、共形生成元和纠缠谱

```bash
julia --project=. bin/mottjain.jl scaling   --config=config/my_run.toml
julia --project=. bin/mottjain.jl generator --config=config/my_run.toml
julia --project=. bin/mottjain.jl oes       --config=config/my_run.toml
julia --project=. bin/mottjain.jl rses      --config=config/my_run.toml
```

这些命令同时保存机器可读 CSV 与 PNG；向量、生成元等对象保存为 JLD2。

## 5. CFT tower score

默认七个关系为：

1. `dS-S = 1`
2. `ddS-dS = 1`
3. `boxS-S = 2`
4. `J = 2`
5. `curlJ = 3`
6. `dJ = 3`
7. `T = 3`

程序先最小二乘拟合能量因子 `factor`，再计算 scaling-dimension 单位的 RMS 残差
`q`。`q` 越小表示这些特定 tower 关系越接近，但它不是“存在 CFT”的充分判据；
仍需结合跨尺寸收敛、量子数稳定性、生成元 overlap 与其它观测量。

新版本有两个重要的分类步骤：

- 能量简并子空间中重新对角化 `L2`，再在同一 `L2` 子空间对角化 `C2`；
- 同一物理 multiplet 在不同 `(Z,R)` 扇区的等能副本先合并，再按 distinct energy
  排 rank。副本数作为 `multiplicity` 保留。

容差由 `energy_tol`、`quantum_tol`、`degeneracy_tol` 控制。若结论对这些容差
非常敏感，应把它视为诊断信号，而不是通过放宽容差隐藏。

## 6. 输出、恢复与可追溯性

```text
output/<run_name>/
├── critical/
├── density/
├── fss/
├── generator/
├── oes/
├── rses/
├── scaling/
└── spectrum/
```

每个目录含 `run_metadata.toml`，记录 Julia 版本、线程数、本项目与 FuzzifiED git
revision。长扫描按 job 追加 checkpoint；失败也会保存错误原因。

不同物理方案请修改 `[output].run_name`，不要把多个方案写入同一目录：

```toml
[output]
root = "/scratch/username/mottjain-output"
run_name = "uf0_scan_v3"
```

## 7. 服务器运行

模板是 `scripts/run_slurm.sh`。按集群修改 partition、内存和 wall time后运行：

```bash
sbatch scripts/run_slurm.sh
```

不要在同一进程并行构造不同 `nm1`。FuzzifiED 的球面半径是全局设置；当前
workflow 按系统大小顺序执行，并为每个 radius-dependent term 显式传
`norm_r2=nm1`。若以后做进程级并行，最安全的是“一种 nm1 一个 Julia 进程”。

## 8. 性能和内存

- μ 扫描：每组 Hamiltonian 只构造一次 `H0`、`Nf`、`L2`、`C2`，随后只做
  `H0 + μNf`。
- FSS：同一 `nm1` 的四个 basis 和观测量跨扫描点复用；每个新扫描值只重建
  一次 `H0`。
- 多参数优化：固定矩阵和每个自由参数方向各构造一次。自由参数越多内存越大；
  大尺寸内存不足时优先减少自由参数。
- 每个离散扇区保留上一次基态，作为下一次 ARPACK 的 warm start。

大尺寸仍可能很慢，因为 Hilbert 空间增长是真实的；重构消除重复构造成本，不能
改变 Hilbert 空间的指数增长。

## 9. 常见问题

**`score_valid=false` 或 missing levels**：先提高 `solver.k`，再检查被拒绝态的
`L2/C2` 偏差，不要直接扩大 μ 优化范围。

**最优 μ 落在边界**：结果会写 `at_boundary=true`。扩大区间后重跑，并检查附近
是否发生 level crossing。

**FSS fit 不可识别**：至少准备三个系统大小；扫描曲线数增加时，每条曲线也增加
一个振幅参数，因此总数据点必须同步增加。

**服务器找不到 FuzzifiED**：运行 `Pkg.develop(path="实际路径")`，随后
`Pkg.instantiate()`，不要在源文件里硬编码服务器路径。

**想直接使用 Julia API**：

```julia
using MottJainED

c = Couplings(Uf0=2.5, Vf0=0.4, mu=0.05)
s = SolverSettings(k=30)
model = build_model(nm1=5)
cache = prepare_spectrum(model, c, s)
states = solve_spectrum(cache, c.mu)
score = cft_score(states; settings=s)
```

加载包只定义功能，不会自动开始任何扫描。
