# MottJainED 使用手册

## 1. 项目目标

`MottJainED` 是旧 `mott_jain` 工作流的独立重构。它保留当前哈密顿量

\[
H=U_f\!\int n_f^2+U_{f0}\!\int n_fn_0+U_0\!\int n_0^2
 +V_f\!\int n_f\nabla^2n_f+V_{f0}\!\int n_0\nabla^2n_f
 +V_0\!\int n_0\nabla^2n_0-t(f_0^\dagger f_1f_2f_3+h.c.)+\mu N_f,
\]

但把模型、数值求解、物理诊断、数据存取和画图分离。项目不会修改旧的
`../mott_jain` 或上游 FuzzifiED 源码。

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
cd MottJainED
julia --project=. scripts/setup.jl
julia --project=. bin/mottjain.jl plan
julia --project=. -e 'using Pkg; Pkg.test()'
```

FuzzifiED 已在 `Project.toml`中固定到核对过的 Git commit。
新服务器只需 clone 本项目，然后运行：

```bash
cd MottJainED
julia --project=. scripts/setup.jl
```

Julia 会下载固定版本的 FuzzifiED，并为服务器平台安装它的 JLL 二进制依赖；
不要求服务器上已有相邻的 `FuzzifiED.jl`。本地Julia 1.11使用仓库中的
`Manifest-v1.11.toml`；服务器Julia 1.12第一次安装时生成自己的
`Manifest-v1.12.toml`。不能让1.12读取1.11的通用Manifest，否则Pkg自身的
OpenSSL/LibSSH2依赖就可能冲突。第一次完整预编译CairoMakie等依赖可能需要
十几分钟；看到持续出现 `✓ package` 是正常进度。

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
fuzzified_threads = 0
blas_threads = 1
```

`fuzzified_threads=0` 表示采用启动 Julia 时的线程数。服务器的专用 Slurm 文件
会用自己的 `threads` 设置它，并保持 BLAS 为 1，避免 FuzzifiED、Julia 和 BLAS
三层线程互相超额占用。由于 `sdicnormal` 的内存随申请 CPU 数增长，实际线程数
可以小于申请 CPU 数。
本地直接运行数值任务时建议加 `--threads=auto`；如果不加，Julia 通常只有
1 个线程，当 `fuzzified_threads=0` 时 FuzzifiED 也会只使用 1 个线程。

## 4. 命令

统一调用形式：

```bash
julia --startup-file=no --threads=auto --project=. bin/mottjain.jl COMMAND --config=config/my_run.toml
```

下面的例子就是本项目目前全部功能的参考，不需要使用程序里的 `--help`。
所有命令都应先 `cd` 到 `MottJainED` 根目录再运行。

### 全部命令速查

| 功能 | `COMMAND` | 主要读取的配置 | 主要结果 |
|---|---|---|---|
| 只检查任务 | `plan` | 全部配置的任务规模 | 只打印预览，不计算、不写数值结果 |
| 单点筛选低能谱 | `spectrum` | `[model] [hamiltonian].mu [solver] [spectrum]` | `spectrum.csv` |
| scalar/J gap | `gap` | `[model].nm_values [hamiltonian] [gap]` | `gap_results.csv`、两张 gap 图 |
| particle density | `density` | `[model].nm1 [hamiltonian] [density]` | `density.csv`、`density.png` |
| 固定 μ 网格找临界点 | `critical` | `[critical]` | `critical_scan.csv`、`critical_point.csv` |
| 优化任意参数组合 | `optimize` | `[optimization]`，可加 profile | `evaluations.csv`、`best.csv` |
| 计算 FSS 数据 | `fss` | `[model].nm_values [fss]` | grid/optimize 两套 CSV |
| FSS 数据和两种标度维数图 | `fss-all` | `[fss]` | FSS CSV、`delta_s`/`delta_o` 图 |
| 只画已有 FSS | `fss-plot` | 已有 FSS CSV | PNG，不做 ED |
| 只拟合已有 FSS | `fss-fit` | 已有 FSS CSV | fit CSV/PNG，不做 ED |
| scaling-dimension 图 | `scaling` | `[model] [hamiltonian] [solver] [scaling]` | scaling CSV/PNG |
| 登记候选参数 | `generator-register` | optimization 的 `best.csv` | `config/generator_points.csv` 新增一行 |
| 固定 ED 与生成元 | `generator` | 全局候选点、`config/generator/<point>/generator_fit.toml` | ED 快照与固定 Lambda |
| 试 tower/overlap | `tower` | 已有 ED/Lambda、`config/generator/<point>/tower.toml` | 每次选态和 overlap 结果 |
| 轨道纠缠谱 | `oes` | `[model] [hamiltonian] [solver] [entanglement]` | `oes/` 下 CSV/PNG |
| 实空间纠缠谱 | `rses` | 同上 | `rses/` 下 CSV/PNG |

#### 通用写法和临时参数

最普通的写法：

```bash
julia --threads=auto --project=. bin/mottjain.jl COMMAND \
  --config=config/my_run.toml
```

常见的临时参数如下。不是每个参数都适用于每个功能：

| 参数 | 意义 | 适用情况 |
|---|---|---|
| `--config=FILE` | 本次读取的完整任务配置 | 所有数值功能 |
| `--override=FILE` | 在完整配置上叠加一个小 TOML | 给 `critical/fss/optimize` 切换案例 profile |
| `--nm1=N` | 临时替换 `[model].nm1` | `plan/optimize/spectrum/density/critical/scaling/oes/rses` |
| `--k=N` | 临时替换当前功能自己的 k | `plan/spectrum/optimize/gap/density/critical/fss/generator/tower` |
| `--run-name=NAME` | 在当前功能目录下建立 `NAME_01/NAME_02` | 普通功能需要区分多组输入时 |
| `--output=DIR` | 临时替换 `[output].root` | 本地与服务器输出根目录不同时 |
| `--force` | 覆盖已有 checkpoint 或固定结果 | `spectrum/gap/density/fss/generator/tower`，慎用 |

`spectrum` 使用自己的 `[spectrum].k`，也可临时传 `--k=N`。`scaling/oes/rses`
使用通用 `[solver].k`；`gap/density/critical/fss/optimize` 也各有独立 k。

只有少数功能有专用命令参数：

| 功能 | 专用参数 | 意义 |
|---|---|---|
| `fss/fss-all` | `--method=grid\|optimize\|both` | 临时选择本次 μc 方法，不改 TOML |
| `fss-plot/fss-fit` | `--method=...`、`--y=列名`、`--source=CSV` | 选择已有数据和要画/拟合的列 |
| `generator-register` | `--point=ID`、`--from=CSV`、`--notes=文字` | 指定新候选 ID、optimization 结果和备注 |
| `generator-register` | `--replace`、`--registry=CSV` | 明确替换同名行，或改用另一份全局表 |
| `generator` | `--point=ID`、`--ed-only`、`--refit` | 指定候选点，以及只做 ED/重拟合的特殊模式 |
| `generator` | `--registry=CSV`、`--data-root=DIR`、`--fit-config=TOML` | 临时改候选表、快照根目录或 Lambda 训练定义 |
| `tower` | `--point=ID`、`--tower-config=TOML` | 指定已有候选点和本次 tower 选态定义 |
| `tower` | `--registry=CSV`、`--data-root=DIR` | 必须与生成 ED 时使用的位置相同 |

表中没有列出的物理参数都应在 TOML 中修改；例如 Hamiltonian 系数不作为长串
命令参数传入，这样每次运行都能把完整配置留档。

#### plan：只确认配置，不开始计算

```bash
julia --project=. bin/mottjain.jl plan --config=config/my_run.toml
```

预览某个 optimization profile、系统大小和 k：

```bash
julia --project=. bin/mottjain.jl plan \
  --config=config/my_run.toml \
  --override=config/optimization_profiles/mu_uf_uf0_vf0.toml \
  --nm1=5 --k=10
```

#### spectrum、gap、density、critical

```bash
# 单个 nm1，在 [hamiltonian].mu 输出筛选并 rescale 后的物理能级
julia --threads=auto --project=. bin/mottjain.jl spectrum \
  --config=config/my_run.toml

# 对 [model].nm_values 的每个尺寸同时画 scalar gap 和 J gap
julia --threads=auto --project=. bin/mottjain.jl gap \
  --config=config/my_run.toml

# 单尺寸扫描 ⟨Nf⟩
julia --threads=auto --project=. bin/mottjain.jl density \
  --config=config/my_run.toml

# 硬算 [critical] 给出的全部 μ；score 从子 TOML 读取
julia --threads=auto --project=. bin/mottjain.jl critical \
  --config=config/my_run.toml \
  --override=config/critical_profiles/critical5.toml
```

改 μ 范围、点数或 k 时编辑 `my_run.toml` 对应 section；`score_terms` 和
`score_metric` 编辑本次 critical profile。`critical` 不调用优化器。

#### optimize：五种常用自由参数组合

```bash
# 只优化 muc
julia --threads=auto --project=. bin/mottjain.jl optimize \
  --config=config/my_run.toml \
  --override=config/optimization_profiles/mu_only.toml \
  --nm1=5 --k=10

# muc + Uf0
julia --threads=auto --project=. bin/mottjain.jl optimize \
  --config=config/my_run.toml \
  --override=config/optimization_profiles/mu_uf0.toml \
  --nm1=5 --k=10

# muc + Uf0 + V0
julia --threads=auto --project=. bin/mottjain.jl optimize \
  --config=config/my_run.toml \
  --override=config/optimization_profiles/mu_uf0_v0.toml \
  --nm1=5 --k=10

# muc + Uf + Uf0 + Vf0
julia --threads=auto --project=. bin/mottjain.jl optimize \
  --config=config/my_run.toml \
  --override=config/optimization_profiles/mu_uf_uf0_vf0.toml \
  --nm1=5 --k=10

# muc + Uf + Uf0 + Vf0 + V0
julia --threads=auto --project=. bin/mottjain.jl optimize \
  --config=config/my_run.toml \
  --override=config/optimization_profiles/mu_uf_uf0_vf0_v0.toml \
  --nm1=5 --k=10
```

每种组合的 `free`、初值、固定值、bounds、`score_terms` 和 `score_metric` 都在
相应 profile 中修改。结果目录使用可读的连续编号，例如
`output/optimize/mu_uf_uf0_vf0_nm5_k10_01/`；不会再把 hash 放进目录名。
同一配置中断后重跑会复用原编号，配置发生变化才顺延为 `_02`。真正最优的一行
直接是该目录里的 `best.csv`。

#### FSS：grid、optimize、两者都算和纯后处理

```bash
# 按 [fss].methods；默认 grid 和 optimize 都算，分开保存
julia --threads=auto --project=. bin/mottjain.jl fss \
  --config=config/my_run.toml \
  --override=config/fss_profiles/fss7.toml

# 只硬算 μ 网格
julia --threads=auto --project=. bin/mottjain.jl fss \
  --config=config/my_run.toml --override=config/fss_profiles/fss7.toml \
  --method=grid

# 固定每个外层参数值，沿 size 追踪 muc
julia --threads=auto --project=. bin/mottjain.jl fss \
  --config=config/my_run.toml --override=config/fss_profiles/fss7.toml \
  --method=optimize

# 算数据，再分别画 delta_s 和 delta_o；不做拟合
julia --threads=auto --project=. bin/mottjain.jl fss-all \
  --config=config/my_run.toml --override=config/fss_profiles/fss7.toml \
  --method=both

# 不做 ED，只读取 grid CSV 画 delta_s
julia --project=. bin/mottjain.jl fss-plot \
  --config=config/my_run.toml --override=config/fss_profiles/fss7.toml \
  --method=grid --y=delta_s

# 不做 ED，只读取 optimize CSV 拟合 delta_s
julia --project=. bin/mottjain.jl fss-fit \
  --config=config/my_run.toml --override=config/fss_profiles/fss7.toml \
  --method=optimize --y=delta_s

# 也可跳过默认目录，直接指定任意已有 CSV
julia --project=. bin/mottjain.jl fss-plot \
  --config=config/my_run.toml --source=/完整路径/results.csv --y=delta_s
```

#### scaling、OES 和 RSES

```bash
julia --threads=auto --project=. bin/mottjain.jl scaling --config=config/my_run.toml
julia --threads=auto --project=. bin/mottjain.jl oes     --config=config/my_run.toml
julia --threads=auto --project=. bin/mottjain.jl rses    --config=config/my_run.toml
```

三者都使用 `[model].nm1`、`[hamiltonian]` 和 `[solver].k`。scaling 的显示范围和
可选 factor 在 `[scaling]`；OES/RSES 的切分和 cutoff 在 `[entanglement]`。

#### 从 optimize 到 generator 再到 tower 的完整案例

假设 optimize 已生成一个你认为值得继续看的 `best.csv`：

```bash
# 1. 登记为全局候选点；只写 CSV，不做 ED
julia --project=. bin/mottjain.jl generator-register \
  --config=config/my_run.toml \
  --point=nm5_candidate_01 \
  --from=output/optimize/mu_uf_uf0_vf0_nm5_k10_01/best.csv \
  --notes="nm1=5，准备检查 generator"

# 2. 建立/复用 ED，并拟合一次固定 Lambda
julia --threads=auto --project=. bin/mottjain.jl generator \
  --config=config/my_run.toml --point=nm5_candidate_01

# 3. 修改这个 point 自己的 config/generator/nm5_candidate_01/tower.toml 后反复试
julia --threads=auto --project=. bin/mottjain.jl tower \
  --config=config/my_run.toml --point=nm5_candidate_01
```

登记位置永远是 `config/generator_points.csv`（除非显式用 `--registry=...`）。登记时
还会从模板建立 `config/generator/<point_id>/generator_fit.toml` 和 `tower.toml`。
generator 的所有数据统一放在 `output/generator/<point_id>/`；不再出现 snapshot
hash 子目录。`ed_snapshot.jld2` 中才有完整 basis 和本征向量。

只想先保存 ED、不拟合 Lambda：

```bash
julia --threads=auto --project=. bin/mottjain.jl generator \
  --config=config/my_run.toml --point=nm5_candidate_01 --ed-only
```

之后不加 `--ed-only` 再运行同一命令，会复用 ED 并补做 Lambda 拟合。只有你明确
改变了 `config/generator/nm5_candidate_01/generator_fit.toml` 的 S/dS 训练定义并
希望替换旧拟合时，才运行：

```bash
julia --threads=auto --project=. bin/mottjain.jl generator \
  --config=config/my_run.toml --point=nm5_candidate_01 --refit
```

`generator-register --replace` 会替换同名候选点，`generator --force` 会强制重做
ED，`tower --force` 会覆盖相同 `tower_01` 后处理。这三个选项都可能替换已有
结果，平时不要加。

### generator 三个运行开关的区别

generator 的三个关键开关应特别区分：

| 参数 | 会不会 ED | 会不会拟合 Lambda | 用途 |
|---|---:|---:|---|
| 不加额外开关 | 缺快照才 ED | 缺固定 Lambda 才拟合 | 推荐日常用法 |
| `--ed-only` | 缺快照才 ED | 否 | 先只保存/检查本征态 |
| `--refit` | 复用已有 ED | 是，覆盖固定 Lambda | 明确改变 S/dS 训练定义 |
| `--force` | 是，强制重算 | 是，并清除旧 tower | ED 文件疑似损坏或明确替换该 point 时才用 |

### 谱、scalar/J gap 与密度

```bash
julia --project=. bin/mottjain.jl spectrum --config=config/my_run.toml
julia --project=. bin/mottjain.jl gap     --config=config/my_run.toml
julia --project=. bin/mottjain.jl density --config=config/my_run.toml
```

`spectrum` 只计算 `[hamiltonian].mu` 一个点。它把不同离散 `(Z,R)` sector
中同 `(L²,C₂)`、同能量的副本合并，但不同能量始终保留为不同能级；然后对
`[spectrum].l2_values × c2_values` 的每个组合保留最低
`levels_per_block` 个。`spectrum.csv` 每个 `(L²,C₂)` 组合占一列，默认顺序是
`(0,0),(2,0),(6,0),(0,3),(2,3),(6,3)`，每列 7 个 rescaled energy。
不保存本征向量，也不附加算符身份。`[spectrum].factor` 必须明确给出，程序不会
用暂定算符身份替你拟合它。已有结果默认复用；加 `--force` 才重算。

`gap` 在每个 `(nm1,μ)` 只求一次谱，同时计算原始 `(L²,C₂)=(0,0)`
列表第二项对应的 scalar gap，以及最低 `(2,3)` 态对应的 J gap。输出为
`gap_results.csv`、`scalar_gap.png` 和 `j_gap.png`；两张图均为正方形、
带浅色网格，纵轴范围为 `0–1.0`。`[gap].k` 默认为 5，不使用全局
`[solver].k` 的较大设置。

`density` 只保留每个 sector 的基态候选并计算全局基态的
`⟨Nf⟩`，不再对高能态做 `L²/C₂` 分类。`[density].k` 默认为 3；
`density.png` 为 `650×650` 的正方形图。

### 临界 μ

公共的 μ 范围和 k 留在 `my_run.toml`；本次态标准放在子配置：

```toml
# config/critical_profiles/critical5.toml
[critical]
score_terms = ["ds_s", "j", "curlj", "dj_rank1", "t_rank2"]
score_metric = "q"

[output]
run_name = "critical5"
```

```bash
julia --threads=auto --project=. bin/mottjain.jl critical \
  --config=config/my_run.toml \
  --override=config/critical_profiles/critical5.toml
```

程序只计算 `[critical]` 中 `mu_min`、`mu_max`、`mu_count` 生成的
均匀网格，然后选所配置 score 最小的网格点。不调用 `Optim`，不会在两个
网格点之间额外求谱。提供的 `critical5.toml` 就是旧
`find_critical_point.jl` 的五项规则；`[critical].k` 默认为 10。没有选择包含
`score_terms` 的案例 profile 时，程序会拒绝开始计算，避免无意套用错误标准。

### 单参数或多参数优化

下面整段属于一个 optimization profile，而不是公共 `my_run.toml`：

```toml
[optimization]
score_terms = ["ds_s", "dds_ds", "c2_6", "boxs_s", "j", "curlj", "dj_rank3", "t_rank1"]
score_metric = "cost"
algorithm = "auto"
tie_u0_to_uf = true
u0_over_uf = 9.0
free = ["Uf", "Uf0", "Vf0", "V0", "mu"]
max_iterations = 200

[optimization.values]
Uf = 0.5
Uf0 = 3.0
U0 = 4.5
Vf = 0.0
Vf0 = 0.5
V0 = 1.0
t = 0.5
mu = 0.05

[optimization.bounds]
Uf = [0.0, 10.0]
Uf0 = [0.0, 20.0]
U0 = [0.0, 100.0]
Vf = [-20.0, 20.0]
Vf0 = [-20.0, 20.0]
V0 = [-20.0, 20.0]
t = [-10.0, 10.0]
mu = [-100.0, 100.0]
```

运行时必须选中相应 profile，例如：

```bash
julia --project=. bin/mottjain.jl optimize \
  --config=config/my_run.toml \
  --override=config/optimization_profiles/mu_uf_uf0_vf0_v0.toml
```

这里三部分各司其职：

- `free`：哪些参数允许优化器改变；
- `[optimization.values]`：自由参数的初值，以及非自由参数的固定值；
- `[optimization.bounds]`：自由参数允许的范围，非自由参数的边界不会使用。

因此不需要再进源文件注释代码。比如固定 `Uf=0.5,U0=4.5,Vf=0,Vf0=0.3,t=0.5`，
只优化 `Uf0,V0,mu`：

```toml
[optimization]
free = ["Uf0", "V0", "mu"]
tie_u0_to_uf = true

[optimization.values]
Uf = 0.5
Uf0 = 3.0
Vf = 0.0
Vf0 = 0.3
V0 = 1.0
t = 0.5
mu = 0.05
```

把 `Vf` 也加入优化只需写：

```toml
free = ["Uf", "Uf0", "Vf", "Vf0", "V0", "mu"]
```

并确保 `[optimization.bounds]` 中存在 `Vf`。如果要让 `U0` 独立优化，则写
`tie_u0_to_uf=false`，把 `U0` 加入 `free`；若保持 `true`，程序在所有情况下都
严格使用 `U0=u0_over_uf*Uf`，无论 `Uf` 是自由还是固定参数。

你原 `optimization.jl` 中几种注释切换可直接对应为：

| 原方案 | `free` | 需要在 `values` 固定的关键项 |
|---|---|---|
| “超级大满贯” | `Uf,Uf0,Vf0,V0,mu` | `Vf=0,t=0.5` |
| “大满贯” | `Uf,Uf0,Vf0,mu` | `Vf=0,V0=1,t=0.5` |
| “固定小的V” | `Uf0,V0,mu` | `Uf=0.5,Vf=0,Vf0=0.3,t=0.5` |
| “固定大的V” | `Uf0,Vf0,mu` | `Uf=0.5,Vf=0,V0=1,t=0.5` |
| “加Vf” | 在相应方案的 `free` 中再加入 `Vf` | 给 `Vf` 初值和 bounds |

这些方案在 `tie_u0_to_uf=true` 时都会自动令 `U0=9Uf`。

`algorithm="auto"` 时，一个自由参数使用有界 Brent，两个及以上使用旧
`optimization.jl` 的 Nelder–Mead 和越界 penalty。这是独立 `optimize` 命令的
逻辑，没有改变。FSS 内层寻找 μc 则固定每个外层参数值，沿 size 续接上一尺寸的
μc，并用扩窗和宽区间 Brent 候选检查；它仍走专门的 `H0+μNf` 缓存入口，避免
为每个外层点建立不必要的多参数矩阵。

每次 evaluation 都会立即追加到 `evaluations.csv`，同时记录自由参数列表、算法、
八个实际 Hamiltonian 系数、score 和错误原因。

### 用小模板快速切换 optimization

项目已经提供五个覆盖模板：

| 文件 | 自由参数 |
|---|---|
| `config/optimization_profiles/mu_only.toml` | `mu` |
| `config/optimization_profiles/mu_uf0.toml` | `Uf0,mu` |
| `config/optimization_profiles/mu_uf0_v0.toml` | `Uf0,V0,mu` |
| `config/optimization_profiles/mu_uf_uf0_vf0.toml` | `Uf,Uf0,Vf0,mu` |
| `config/optimization_profiles/mu_uf_uf0_vf0_v0.toml` | `Uf,Uf0,Vf0,V0,mu` |

所有模板都固定 `Vf=0`，并且各自保存 `free/values/bounds/score_terms/score_metric`。
运行时先读取公共 `my_run.toml`，再叠加这个案例子配置。同一个模板可直接换系统
大小和 k：

```bash
julia --threads=auto --project=. bin/mottjain.jl optimize \
  --config=config/my_run.toml \
  --override=config/optimization_profiles/mu_uf0_v0.toml \
  --nm1=7 --k=12
```

把 `optimize` 临时换成 `plan` 可以先检查，不会开始求谱：

```bash
julia --project=. bin/mottjain.jl plan \
  --config=config/my_run.toml \
  --override=config/optimization_profiles/mu_uf0_v0.toml \
  --nm1=7 --k=12
```

optimization 的可见目录名只使用“模板名 + nm1 + k + 连续数字”，例如
`mu_uf0_v0_nm7_k12_01`。完整配置签名只隐藏保存在 `case_identity.toml` 中，用于
判断重跑应该续接 `_01` 还是新建 `_02`，不会再显示成乱码。要自己指定可读前缀时
可加 `--run-name=my_name`，得到 `my_name_01`。

`mu_only.toml` 是独立 `optimize` 命令只找 μ 的模板；FSS 扫描中的 μc 仍直接用
`fss --method=optimize`，因为它还包含 `nm_values × scan_values` 的外层循环。

### Finite-size scaling

公共扫描范围留在 `my_run.toml`：

```toml
[model]
nm_values = [4, 5, 6, 7]

[fss]
methods = ["grid", "optimize"]
k = 15
scan_parameter = "Uf0"
scan_values = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
mu_min = 0.0
mu_max = 0.12
mu_count = 9
optimize_strategy = "size_continuation"
optimize_anchor_nm = 4
optimize_local_half_width = 0.02
optimize_local_count = 9
optimize_max_expansions = 3
optimize_abs_tol = 1.0e-4
optimize_max_iterations = 60
```

本次 score 放在 `config/fss_profiles/fss7.toml`：

```toml
[fss]
# 旧 FSS1.jl 以 k=15 为基线；当前扩展 Hamiltonian 的完整扫描使用 k=30，
# 才能覆盖 Uf0=1.5 等参数点需要的高能级。不要使用测试用的 k=3。
k = 30
score_terms = ["ds_s", "dds_ds", "boxs_s", "j", "curlj", "dj_rank1", "t_rank1"]
score_metric = "cost"

[output]
run_name = "fss7"
```

另有 `config/fss_profiles/fss5.toml`，使用目前统一的五项
`["ds_s", "j", "curlj", "dj_rank1", "t_rank1"]`、`score_metric="cost"`
和 `k=15`。要跑这一套时，只把下面命令中的 `fss7.toml` 换成
`fss5.toml`；它会写到独立的 `output/fss/fss5_XX/`，不会与七项结果混合。

```bash
# 未加 --method：按 TOML 的 methods；默认两种都算
julia --threads=auto --project=. bin/mottjain.jl fss \
  --config=config/my_run.toml --override=config/fss_profiles/fss7.toml

# 本次命令只算固定 μ 网格
julia --threads=auto --project=. bin/mottjain.jl fss \
  --config=config/my_run.toml --override=config/fss_profiles/fss7.toml --method=grid

# 本次命令固定外层参数，并沿 size 追踪 μc
julia --threads=auto --project=. bin/mottjain.jl fss \
  --config=config/my_run.toml --override=config/fss_profiles/fss7.toml --method=optimize
```

两种方法都生成每个 `(nm1, scan_value)` 的结果。区别在 optimize 会把同一个
`scan_value` 在不同 size 上连成一条独立的 μc 轨迹：

- `grid`：准确计算 `mu_count` 个 μ，取 score 最小的网格点；
- `optimize`：`N < optimize_anchor_nm` 独立宽搜但不作为后续 seed；
  `N = optimize_anchor_nm` 宽搜并建立 anchor；更大的 N 围绕前一 size 的 μc 做
  局部细网格和逐谷底 Brent。局部网格最低在窗口边界时自动扩大窗口，同时另跑
  一次全范围 Brent 作为候选；最终比较所有实际算过的有效点，选择最低 score。

程序内部仍按 size 在外层执行，以便同一 size 的不同 `scan_value` 复用昂贵的
Basis/Hamiltonian 缓存；但上一 size 的 μc 按 `scan_value` 分开保存，不会拿
`Uf0=1.5` 的结果给 `Uf0=2.2` 当 seed。

结果不会混在一起，分别保存为 `fss_grid_results.csv` 和
`fss_optimize_results.csv`。`optimize` 的每次 μ 求值还会立即写入
`fss_optimize_evaluations.csv`，可以直接检查完整的 `q(μ)` 轨迹。因此两种结果
可以直接比较，也可以分别画图。结果行会记录 `search_mode`、上一尺寸给出的
`center_muc`、实际局部窗口、扩窗次数、宽 Brent 候选以及最后的 `best_source`：

```bash
julia --project=. bin/mottjain.jl fss-plot \
  --config=config/my_run.toml --override=config/fss_profiles/fss7.toml --method=grid --y=delta_s
julia --project=. bin/mottjain.jl fss-plot \
  --config=config/my_run.toml --override=config/fss_profiles/fss7.toml --method=optimize --y=delta_s
```

若要改扫 `Vf0`，把 `scan_parameter` 改成 `"Vf0"` 并相应修改
`scan_values`；其余 Hamiltonian 系数取 `[hamiltonian]` 中的固定值。建议复制一份
FSS profile 并修改其中的 `run_name`，不要与 Uf0 扫描混写。

当前 no-W 候选已经提供两份这样的独立 profile：
`config/fss_profiles/no_w_uf0.toml` 和 `no_w_vf0.toml`。它们分别复现旧
`FSS1.jl` 先扫 Uf0、再扫 Vf0 的结构；这是两组一维扫描，不是 Uf0×Vf0 二维网格。

`fss-all` 会在数据计算后，分别对 grid/optimize 结果画 `delta_s` 和 `delta_o`，
但不自动拟合。只有明确调用独立的 `fss-fit` 命令时，才会尝试下面的联合拟合：

\[
\Delta(N,g)=\Delta_\infty+a_gN^{-\omega/2}
              =\Delta_\infty+a_gx^\omega,\qquad x=N^{-1/2},
\]

所有扫描参数共享 `Delta_inf` 与 `omega`，每条曲线有独立振幅 `a_g`。至少需要
三个系统大小，且数据点必须多于参数数目；条件不足时会明确报错，不会生成一个
看似正常但实际上欠定的拟合。

### Scaling dimension 和纠缠谱

```bash
julia --project=. bin/mottjain.jl scaling   --config=config/my_run.toml
julia --project=. bin/mottjain.jl oes       --config=config/my_run.toml
julia --project=. bin/mottjain.jl rses      --config=config/my_run.toml
```

这些命令同时保存机器可读 CSV 与 PNG；向量、生成元等对象保存为 JLD2。

### 共形生成元：全局参数表、ED 快照与 tower 后处理

generator 不再直接使用 `[model].nm1` 和 `[hamiltonian]`。它只接受全局参数表
`config/generator_points.csv` 中明确选中的一行，防止把 optimization A 的 `muc`
和 optimization B 的其它耦合混在一起。

参数表的核心列是：

```text
point_id, enabled, nm1,
Uf, Uf0, U0, Vf, Vf0, V0, t, muc,
factor, objective, q, cost, delta_s, delta_o,
score_definition, source, notes
```

其中七个非化学势 Hamiltonian 系数是 `Uf,Uf0,U0,Vf,Vf0,V0,t`；`muc` 单独保存。
`factor` 和各种 score 不是 ED 输入，但会跟随这个候选点进入后处理输出。

#### 方法一：从 optimization 的 best.csv 登记，不手抄数值

先看完某次 optimization，确认它值得研究，然后运行：

```bash
julia --project=. bin/mottjain.jl generator-register \
  --config=config/my_run.toml \
  --point=nm6_candidate_01 \
  --from=output/optimize/mu_uf0_v0_nm6_k70_01/best.csv \
  --notes="nm1=6，目前 tower 看起来最好"
```

这条命令不做 ED。它向全局 CSV 加一行，把 `mu_initial` 改名为更明确的 `muc`，
复制七个耦合、factor、score 和来源路径；同时建立
`config/generator/nm6_candidate_01/` 下该点自己的两个 TOML。同名 `point_id` 默认
拒绝覆盖，避免误写；只有明确加入 `--replace` 才替换。

#### 方法二：手动编辑全局 CSV

也可以直接复制 `example_disabled` 那一行，修改：

- `point_id`：给它一个唯一且稳定的名字；
- `enabled=true`；
- `nm1`、七个耦合和 `muc`；
- 已知时填写 `factor`，不知道可以留空；
- `notes` 写下为什么保留它。

程序会严格检查必需数值、重复 point ID、NaN/Inf 和 disabled 状态，不会缺参数时
偷偷回退到 TOML 的 `[hamiltonian]`。

#### 第一阶段：固定 ED 和 microscopic Lambda

运行：

```bash
julia --threads=auto --project=. bin/mottjain.jl generator \
  --config=config/my_run.toml --point=nm6_candidate_01
```

它先保存普通四个 `(Z,R)` sector，并在 `include_adjoint=true` 时保存旧
`for_generator_special` 使用的 SU(3) adjoint weight sector。每个本征态的完整向量和
basis 都在 `ed_snapshot.jld2`，普通表格另存为：

- `physical_levels_standard.csv`：合并等能副本后的 physical rank/member；
- `physical_levels_adjoint.csv`：same-angular 检验所需额外 sector；
- `snapshot_metadata.toml`：模型尺寸、全部 Hamiltonian 系数、solver 精度、basis 维数、
  Julia/MottJainED/FuzzifiED 版本与源码签名。

逐 sector 原始能量、point 参数、basis 和本征向量已经完整包含在
`ed_snapshot.jld2`；不再额外生成重复的 `spectrum_*.csv` 和 `point.csv`。挑选
tower 态通常只需查看两张 `physical_levels_*.csv`。

随后 generator 单独读取 `config/generator/nm6_candidate_01/generator_fit.toml`，
固定其中指定的训练态 `S→dS`，构造 18 个 microscopic 候选并拟合一次 Lambda。
结果保存在该 point 输出目录下的 `generator/`：

- `generator_fit.jld2`：固定 Lambda Terms、系数、奇异值和数值秩；
- `generator_coefficients.csv`、`generator_summary.csv`；
- `generator_selected_states.csv`：训练时真正使用的 S/dS。

训练配置快照直接存进 `generator_fit.jld2`，不再重复生成两份 generator TOML
元数据。

普通重跑会复用这个 Lambda。如果你后来修改了 `generator_fit.toml`，程序会拒绝
静默替换；只有明确加入 `--refit` 才会覆盖。`--ed-only` 仍可用于只保存 ED，
但推荐的正常流程是不加它，让 generator 同时固定 Lambda。

#### 第二阶段：反复试 tower 态，不重算 ED 或 Lambda

```bash
julia --threads=auto --project=. bin/mottjain.jl tower \
  --config=config/my_run.toml --point=nm6_candidate_01
```

它读取 `config/generator/nm6_candidate_01/tower.toml`。其中：

- `[states.S]` 等表定义 `family,l2,c2,rank,member`；
- `rank` 是 `physical_levels_*.csv` 中合并副本后的 `physical_rank`；
- `member` 选择该能级内部哪个副本，也可改用 `z=...`、`r=...`；
- 每个 `[[overlaps]]` 定义 input、一个或多个 targets、目标角动量和计算模式；
- 多个 targets 的 `total_overlap` 是进入整个候选子空间的权重，比强行认定某一个
  近简并态更稳健；
- `mode="same_angular"` 实现旧代码的 `L- → Lambda_z → L+` 检验；
- `required=false` 的可疑关系缺态时只写 `skipped`，不会丢掉其它结果。

每次修改 tower TOML 会得到可读的 `tower_01`、`tower_02` 连续目录，旧尝试不会
覆盖。内部一致性签名只写进隐藏元数据。输出包括：

- `selected_states.csv`：本次每个名字实际选中了哪个 physical rank/member/sector；
- `tower_overlaps.csv`：逐关系、逐 target overlap 与 total overlap；
- `analysis_metadata.toml`：point、固定 generator、factor 和本次 tower 配置快照；
- `tower_analysis.jld2`：仅当 `save_generated_vectors=true` 时才额外保存生成后向量。

运行结束后，`tower_overlaps.csv` 的主要内容还会按照旧
`conformal_generator.jl` 的固定格子直接打印到终端。每一行显示 input、目标
`l'`、最多若干个 `Target(dE/f) + overlap` 和 total overlap；常用的一到两个
target 会保持原来的两列宽度。`mode`、relation 名和原始能量仍完整保存在 CSV
中。重跑并复用已有 `tower_01` 时也会再次打印，不需要手动打开 CSV。

`tower` 严格不做 ED，也不重新拟合 generator；缺少任意一个固定文件都会明确
报错并要求先运行 `generator --point=...`。

可见目录只使用 `point_id`。程序仍在元数据里保存由 `nm1 + 八个 Hamiltonian 系数
+ k/容差 + adjoint 设置 + ED 源码版本` 生成的隐藏签名。若同一 point 下已有 ED
却与当前设置冲突，程序会拒绝复用；通常应换一个新的 point_id。只有明确使用
`generator --force` 才会替换这个 point 的旧 ED。
由于旧 Lambda 和 tower 不可能再与新本征向量一一对应，`--force` 会同步清除该
point 下面旧的 `generator/` 与 `tower/`，随后重新拟合 Lambda。

## 5. CFT tower score

这一节是一个仍待物理确认的独立分析标准；基础 `spectrum` 和 `scaling` 不会
再自动把它当作已确定的结论。`critical/optimize/fss` 仍以它为目标时，必须把
结果视为依赖当前态认定的试验性输出。

程序内部只有一个候选关系池。凡是会最小化 CFT score 的功能，都能在自己的案例
子 TOML 里用 `score_terms` 任意选择其中若干项；公共 `my_run.toml` 不再保存它。
三个旧名称只是省事的预设：

| 预设名称 | 原始文件 | 默认使用位置 | 默认组合 | 默认 metric |
|---|---|---|---|---|
| `critical5` | `find_critical_point.jl` | `[critical]` | 5 项 | `q` |
| `fss7` | `FSS1.jl` | `[fss]` | 7 项 | `cost` |
| `optimization8` | `optimization.jl` | `[optimization]` | 8 项 | `cost` |

日常复制相应 profile 后编辑 `score_terms`。如果完全删掉这一行，也可以用
`score="critical5"`、`"fss7"` 或 `"optimization8"` 载入整套旧预设。
比如 critical 只选择四项：

```toml
[critical]
score_terms = ["ds_s", "j", "curlj", "t_rank2"]
score_metric = "q"
```

完整候选池如下。`A-B` 表示相应原始态能量之差，`E0` 是全局基态：

| `score_terms` 名称 | 数值关系 | target |
|---|---|---:|
| `ds_s` | `(2,0)[1]-(0,0)[2]` | 1 |
| `dds_ds` | `(6,0)[2]-(2,0)[1]` | 1 |
| `c2_6` | `(6,6)[1]-(2,6)[1]` | 1 |
| `boxs_s` | `(0,0)[3]-(0,0)[2]` | 2 |
| `boxo_o` | `(0,3)[2]-(0,3)[1]` | 2 |
| `j` | `(2,3)[1]-E0` | 2 |
| `curlj` | `(2,3)[3]-E0` | 3 |
| `dj_rank1` | `(6,3)[1]-E0` | 3 |
| `dj_rank3` | `(6,3)[3]-E0` | 3 |
| `t_rank1` | `(6,0)[1]-E0` | 3 |
| `t_rank2` | `(6,0)[2]-E0` | 3 |

`boxo_o` 就是旧 `optimization.jl` 中曾写出但注释掉的 `□O-O`。方括号都是
旧数组的原始 rank，不是合并 multiplet 后的 distinct rank。三套旧预设仍逐项
复现原文件，只是现在可以在任何功能里自由增删。

所有规则都会同时记录两个诊断值。`q` 是先拟合能量因子后，在 scaling-dimension
单位中的 RMS 残差；`cost=||u||²-(u·v)²/||v||²` 是旧 FSS/optimization 使用的
方向夹角目标。`score_metric` 决定程序真正最小化哪一个，另一项仍写入 CSV。
这些 score 都不是“存在 CFT”的充分判据；仍需结合跨尺寸收敛、量子数稳定性、
生成元 overlap 与其它观测量。

不要只选一项：因为 `factor` 也是由同一批关系拟合的，单项时 `q` 和 `cost`
都会恒等于零，不能用于寻找临界点。至少要有两个互相独立的关系，实际优化建议
保留更多约束。减少 `score_terms` 本身不会显著缩短对角化；若所选项不再需要高
rank，可以再谨慎降低该 section 的 `k`，但必须先确认所有所需态都稳定出现。

量子数识别与 rank 规则：

- 能量简并子空间中重新对角化 `L2`，再在同一 `L2` 子空间对角化 `C2`；
- `cft_score` 为了复现旧文件中的数组索引，保留不同 `(Z,R)` 扇区的等能副本，
  再按原始能量顺序取 rank；不使用 `level_catalog` 合并后的 distinct rank。

容差由 `energy_tol`、`quantum_tol`、`degeneracy_tol` 控制。若结论对这些容差
非常敏感，应把它视为诊断信号，而不是通过放宽容差隐藏。

## 6. 输出、恢复与可追溯性

```text
output/
├── spectrum/                  # 未写 run_name 时直接存这里
├── gap/
├── density/
├── critical/<case_name>_01/
├── fss/<case_name>_01/
├── optimize/
│   └── mu_uf_uf0_vf0_nm5_k10_01/
│       ├── evaluations.csv
│       └── best.csv
├── scaling/
├── oes/
├── rses/
└── generator/
    └── <point_id>/
        ├── register.log
        ├── generator.log
        ├── tower.log
        ├── ed_snapshot.jld2
        ├── physical_levels_standard.csv
        ├── physical_levels_adjoint.csv
        ├── snapshot_metadata.toml
        ├── generator/
        │   ├── generator_fit.jld2
        │   ├── generator_coefficients.csv
        │   ├── generator_summary.csv
        │   └── generator_selected_states.csv
        └── tower/
            ├── tower_01/       # overlaps + selected states + 一份 metadata
            └── tower_02/
```

每个目录含 `run_metadata.toml`，记录 Julia 版本、线程数、本项目与 FuzzifiED git
revision。长扫描按 job 追加 checkpoint；失败也会保存错误原因。

数值命令的普通 `@info` 进度默认不再刷终端，而是追加到该任务目录的
`run.log`。终端启动时只显示一次日志路径，之后只保留 warning、error 和未捕获
异常。generator 三步分别使用 `register.log/generator.log/tower.log`；tower 的最终
overlap 表是有意保留的终端短输出。optimization
的每次参数 evaluation 仍会写 `evaluations.csv`；关闭终端 trace 不会丢失优化历史。
`plan` 是给人直接阅读的短输出，仍显示在终端。

普通功能如需区分多个案例，可以在其子 profile 写 `run_name`。同一功能相关配置
不变时重跑会复用 `_01`；相关配置改变后自动顺延 `_02`：

```toml
[output]
root = "/scratch/username/mottjain-output"
run_name = "uf0_scan_v3"
```

## 7. 服务器运行

重计算功能分别使用 `slurm/` 中自己的作业文件。每个文件顶部有独立的
partition、CPU 和 wall time，下面“用户配置区”有该功能自己的 config、
profile、point 和实际 `threads`。修改一次后可以反复提交，不需要重新拼命令：

```bash
sbatch slurm/spectrum.sbatch
sbatch slurm/gap.sbatch
sbatch slurm/density.sbatch
sbatch slurm/critical.sbatch
sbatch slurm/fss.sbatch
sbatch slurm/optimize.sbatch
sbatch slurm/scaling.sbatch
sbatch slurm/generator.sbatch
sbatch slurm/tower.sbatch
sbatch slurm/oes.sbatch
sbatch slurm/rses.sbatch
```

例如 FSS 只需在 `slurm/fss.sbatch` 中保留自己的 `config/profile/method`；优化的
profile、`nm1`、`k` 则只放在 `slurm/optimize.sbatch`，二者互不影响。资源默认值
只是试跑起点，增大系统时仍需按实际峰值修改相应文件。完整说明见
`slurm/README.md`。

当前 `sdicnormal` 的 `DefMemPerCPU=7824 MB`，64 CPU 节点总内存约513024 MB。
专用脚本因此不写 `--mem`：16 CPU约得到122 GiB，32 CPU约245 GiB，56 CPU约
428 GiB。若为了内存申请56 CPU，可以仍写 `threads=16`；没有必要强迫
FuzzifiED 开56线程。更多线程不一定更快，尤其在稀疏矩阵乘法受内存带宽限制时。

`plan`、`generator-register`、`fss-plot`、`fss-fit` 不做昂贵 ED，通常直接在登录
节点运行，不为它们单独申请计算节点。

### 中断以后会不会从头计算

同一份相关配置会回到同一个输出目录；不要加 `--force`。不同功能的断点粒度如下：

| 功能 | 已完成部分是否复用 | 中断点的代价 |
|---|---|---|
| `spectrum` | 完整 `spectrum.csv` 复用 | 单点 ED 整体重算 |
| `gap` | 按每个 `(nm1, μ)` 跳过 | 正在计算的组合重算 |
| `density` | 按每个 μ 跳过 | 正在计算的那个 μ 重算 |
| `fss`/`fss-all` | 按 `(method,nm1,scan_value)` 跳过，并恢复已完成 size 的 muc 供后续续接 | 当前 size 的局部搜索重算 |
| `optimize` | 保留每次 evaluation，并从历史最佳点重新启动 | 不保存 Nelder–Mead simplex/Brent 内部状态，可能重复一些点 |
| `generator` | 完整 `ed_snapshot.jld2` 和 Lambda 会复用 | 若在 ED 快照写完前中断，该 point 的 ED 整体重算 |
| `tower` | 完整且配置签名相同的结果直接复用 | 未完成的 tower 本次整体重算，但不重做 ED/Lambda |
| `critical` | 当前没有逐 μ checkpoint | 中断后整个 μ 网格重算 |
| `scaling`、`oes`、`rses` | 当前是单次整体任务 | 中断后整体重算 |

CSV/JLD2 的正式文件使用原子写入，避免被中断时留下一个看似完整的半文件。
已经成功的离散 job 由参数生成稳定 ID；修改 `k`、系统大小或相关 Hamiltonian
参数后，ID/案例目录会变化，不会把旧结果误当成新结果。

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

### 大系统怎样选择 `k`

`k` 是每个离散 symmetry sector 请求的低能本征态数，不是全体系总态数。
`score_valid=true` 只说明当前 `k` 已经找到了所选 relation 要求的 sector/rank，
是必要条件但不是收敛证明。可靠做法是在最大 `nm1` 的代表性参数点做两到三次：

1. 先用预期值，例如 `k=10`；
2. 再用 `k=15` 或 `20`；
3. 比较被 score 选中的各态能量/rank、`factor`、`objective` 和最终 `muc`；
4. 连续两次增加 `k` 后这些量在所需精度内不变，才把较小的那个 `k` 用于大扫描。

不能因为五项 score 的最大 raw rank 是3，就断言 `k=3` 足够：`k` 是求解器在原始
离散 sector 中保留的态数，之后还要按近似 `(L²,C₂)` 分类，目标表示前面可能夹着
其他态。反过来，`k` 也不必机械地随 `nm1` 成比例增长，应由上述收敛比较决定。
提高 `k` 后，旧的较小 `k` 本征系统不能直接补算成更大的 `k`；它会作为一个新
数值案例重新对角化，但旧数据仍保留，不会被覆盖。

大尺寸仍可能很慢，因为 Hilbert 空间增长是真实的；重构消除重复构造成本，不能
改变 Hilbert 空间的指数增长。

## 9. 常见问题

**`score_valid=false` 或 missing levels**：先提高 `solver.k`，再检查被拒绝态的
`L2/C2` 偏差，不要直接扩大 μ 优化范围。

**最优 μ 落在边界**：结果会写 `at_boundary=true`。扩大区间后重跑，并检查附近
是否发生 level crossing。

**FSS fit 不可识别**：至少准备三个系统大小；扫描曲线数增加时，每条曲线也增加
一个振幅参数，因此总数据点必须同步增加。

**服务器找不到 FuzzifiED**：在项目根目录重新运行
`julia --project=. scripts/setup.jl`。若下载失败，应在允许联网的登录节点完成
setup；网络超时可以重新执行同一命令。不要在源文件里硬编码服务器路径，也不要改回本地
`Pkg.develop`，否则服务器和本地可能使用不同源码。

**想直接使用 Julia API**：

```julia
using MottJainED

c = Couplings(Uf0=2.5, Vf0=0.4, mu=0.05)
s = SolverSettings(k=30)
model = build_model(nm1=5)
cache = prepare_spectrum(model, c, s)
states = solve_spectrum(cache, c.mu)
score = cft_score(
    states; settings=s,
    terms=[:ds_s, :j, :curlj, :dj_rank1, :t_rank2], metric=:q,
)
```

加载包只定义功能，不会自动开始任何扫描。
