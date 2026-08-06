# 从这里开始：MottJainED 文件结构、概念和完整工作流程

这份说明假设你此前没有使用过 Julia 项目、TOML、命令行入口或软件包结构。
阅读完以后，你应该能够回答：

1. 每个文件为什么存在；
2. 哪些文件需要自己修改，哪些不要动；
3. 一条命令怎样一步步变成 Hamiltonian、本征求解和输出文件；
4. 当前有哪些功能，每个功能经过哪些函数；
5. 想检查或修改某一步时应打开哪个文件。

---

## 一、先区分三个完全不同的东西

项目中都有 `.toml` 或 `.jl` 文件，但作用不同。

### 1. Julia 软件环境：`Project.toml` 和 `Manifest.toml`

这两个文件回答的是：

> “运行这套代码需要安装哪些 Julia 软件包，以及具体使用哪个版本？”

它们不保存 (U_f)、(mu)、`nm1`，也不执行物理计算。

### 2. 物理任务配置：`config/default.toml` 或 `config/my_run.toml`

这类文件回答的是：

> “这一次实验用什么 Hamiltonian、系统大小、μ 范围和本征态数？”

它们也不是代码，只是程序读取的数据。

### 3. Julia 代码：所有 `.jl` 文件

`.jl` 文件才定义：

- 怎样构造 fuzzy-sphere model；
- 怎样构造 Hamiltonian；
- 怎样求本征态；
- 怎样识别 (L^2,C_2)；
- 怎样做 FSS、画图和保存数据。

最重要的区分是：

| 类型 | 管什么 | 日常是否修改 |
|---|---|---|
| `Project.toml` | 直接软件依赖及允许版本 | 通常不改 |
| `Manifest.toml` | 所有依赖的精确版本 | 不手动改 |
| `config/*.toml` | 物理和数值参数 | 经常复制后修改 |
| `src/*.jl` | 计算方法和程序逻辑 | 需要改变算法时才改 |

---

## 二、TOML 是什么文件类型

TOML 是一种纯文本配置格式，全名是 **Tom's Obvious Minimal Language**。
它的目标是让人和程序都容易读取。

最简单的 TOML：

```toml
[model]
nm1 = 5
nm_values = [4, 5, 6]

[hamiltonian]
Uf = 0.5
mu = 0.05
```

含义：

- `[model]` 创建一个叫 `model` 的区域；
- `nm1 = 5` 是名字—数值配对；
- `[4, 5, 6]` 是一个数组；
- `#` 后面是注释，程序忽略。

Julia 读入后，大致得到：

```julia
Dict(
    "model" => Dict("nm1" => 5, "nm_values" => [4, 5, 6]),
    "hamiltonian" => Dict("Uf" => 0.5, "mu" => 0.05),
)
```

TOML 本身没有“执行”能力。写：

```toml
mu = 0.05
```

不会开始计算；必须由 Julia 代码读取它。

### 为什么三种文件都使用 TOML

因为 TOML 只是通用格式，文件名和读取者决定用途：

- Julia 的包管理器读取 `Project.toml`、`Manifest.toml`；
- MottJainED 的 `load_config` 读取 `config/my_run.toml`。

它们格式相同，含义不同，就像多个 CSV 可以分别保存能谱和密度。

---

## 三、`Project.toml` 做什么

打开根目录的 `Project.toml`，前四行是项目身份：

```toml
name = "MottJainED"
uuid = "..."
authors = ["..."]
version = "0.1.0"
```

- `name`：因此代码中可以写 `using MottJainED`；
- `uuid`：Julia 用它唯一识别这个项目，避免同名冲突；
- `version`：当前项目版本。

`[deps]` 列出代码直接需要的软件包：

```toml
[deps]
CSV = "..."
CairoMakie = "..."
FuzzifiED = "..."
Optim = "..."
```

右边不是路径，也不是密码，而是每个 Julia 包的唯一 UUID。

`[sources]` 告诉 Julia：

```toml
[sources]
FuzzifiED = {path = "../FuzzifiED.jl"}
```

这里的 FuzzifiED 使用本地相邻文件夹，不从网络下载。

`[compat]` 是允许的版本范围。例如：

```toml
Optim = "1, 2"
```

表示 Optim 1.x 或 2.x 可以使用。

### 什么时候修改 `Project.toml`

只有增加/删除代码依赖时才需要。例如未来代码新增：

```julia
using HDF5
```

才需要通过 Julia `Pkg.add("HDF5")` 更新项目。改变 (U_f) 不应修改它。

### 为什么使用“包”的文件结构

这里的“包”不等于要把代码发布到网上，也不等于把计算藏起来。
它只是 Julia 约定的一种整理方式：

- `src/MottJainED.jl` 统一加载所有函数，避免每个脚本反复 `include`；
- `Project.toml` 统一记录依赖，避免本地和服务器各缺不同的软件包；
- `using MottJainED` 把函数放在独立命名空间，避免全局变量互相覆盖；
- `test/runtests.jl` 可以直接检查同一套代码；
- Julia 可以预编译不常变化的代码，减少重复加载成本。

CLI 只是这套包外面的一个方便入口，不是唯一用法。如果以后想在 Julia REPL
逐步研究，也可以：

```julia
using MottJainED
model = build_model(nm1=5)
cache = prepare_spectrum(model, Couplings(), SolverSettings())
states = solve_spectrum(cache, 0.05)
```

也就是说：包结构负责整理和复用，CLI 负责让一次完整任务可以用一条命令启动。

---

## 四、`Manifest.toml` 做什么

`Manifest.toml` 第一行已经写明：

```toml
# This file is machine-generated - editing it directly is not advised
```

它是 Julia 自动生成的“精确软件清单”。

`Project.toml` 可能只说：

```text
需要 CairoMakie 0.15
```

但 CairoMakie 又依赖 Makie、Colors、GeometryBasics 等许多包。
`Manifest.toml` 会记录整个依赖树的：

- 包名；
- 精确版本；
- UUID；
- 下载内容哈希；
- 包之间的依赖关系；
- 本地 FuzzifiED 路径；
- Julia 版本。

这就是为什么它很长。你不需要读懂它，也不要手工加入注释，因为下一次
`Pkg.resolve()` 可能重写它。

服务器执行：

```bash
julia scripts/setup.jl
```

时，`Pkg.instantiate()` 会根据 Manifest 尽量还原本地相同的软件版本。

可以把二者类比为：

| 文件 | 类比 |
|---|---|
| `Project.toml` | 我需要“面粉、鸡蛋、牛奶” |
| `Manifest.toml` | 面粉品牌/批次、鸡蛋规格、牛奶版本及所有供应链细节 |

---

## 五、整个文件夹分别是什么

```text
MottJainED/
├── Project.toml             Julia 直接依赖与项目身份
├── Manifest.toml            Julia 自动生成的精确依赖锁定
├── README.md                项目首页和文档入口
├── config/
│   └── default.toml         默认物理/数值参数模板
├── bin/
│   └── mottjain.jl          终端命令的短入口
├── scripts/
│   ├── setup.jl             新机器第一次安装依赖
│   └── run_slurm.sh         Slurm 服务器任务模板
├── src/
│   ├── MottJainED.jl        包总入口，按顺序加载其它 src 文件
│   ├── Types.jl             所有核心数据类型
│   ├── Model.jl             物理模型与 Hamiltonian Terms
│   ├── Spectrum.jl          basis/matrix 缓存、本征求解、量子数分类
│   ├── CFT.jl               CFT tower score 与临界 μ 优化
│   ├── Storage.jl           CSV、metadata、job ID、断点续跑工具
│   ├── Workflows.jl         spectrum/gap/density/critical/FSS 等完整任务
│   ├── Conformal.jl         共形生成元候选、拟合与 overlap
│   ├── Entanglement.jl      OES 和 RSES
│   └── CLI.jl               把命令名字分派给 Workflows
├── test/
│   └── runtests.jl          自动测试
└── docs/
    ├── START_HERE.md        当前这份概念和代码地图
    └── USER_GUIDE.md        参数与命令参考
```

---

## 六、`scripts` 中两个文件做什么

### `scripts/setup.jl`：只负责安装环境

它的调用链是：

```text
运行 setup.jl
  → Pkg.activate(MottJainED 根目录)
  → Pkg.develop(../FuzzifiED.jl)
  → Pkg.instantiate()
  → Pkg.precompile()
  → 结束
```

它不会调用 `build_model`，不会对角化，不会产生物理数据。

使用时机：

- 新服务器第一次复制项目后；
- 删除 Julia 环境后；
- Project/Manifest 发生重要变化后。

不是每次 spectrum/FSS 前都运行。

### `scripts/run_slurm.sh`：申请服务器资源后运行一个命令

这是 Bash/Slurm 脚本，不是 Julia 文件。

前面的：

```bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
```

向 Slurm 申请：

- 16 CPU；
- 64 GB 内存；
- 最长 48 小时。

脚本最后：

```bash
julia ... bin/mottjain.jl fss-all --config=config/default.toml
```

才是真正启动 Julia 计算的地方。

如果服务器不用 Slurm，就不使用这个文件。

---

## 七、代码是怎样被加载的

先把常见命令逐块拆开：

```bash
julia --project=. bin/mottjain.jl spectrum --config=config/my_run.toml
```

| 片段 | 作用 |
|---|---|
| `julia` | 启动 Julia 程序 |
| `--project=.` | 把当前文件夹的 Project/Manifest 选作软件环境 |
| `bin/mottjain.jl` | 告诉 Julia 接下来执行哪个 `.jl` 文件 |
| `spectrum` | 传给该文件的第一个普通参数，表示选择求谱功能 |
| `--config=...` | 再传一个选项，指出本次物理参数文件 |

所以 `--project=.` 和 `bin/mottjain.jl` 不重复：前者选择“用哪些软件包”，
后者选择“运行哪段程序”。

`bin/mottjain.jl` 内部还会 `Pkg.activate` 一次项目根目录。在从项目根目录运行
时，这与 `--project=.` 确实选择了同一个环境；保留这一层检查是为了从别的目录
用绝对路径启动入口时仍不会选错环境。因此下面两种方式都可以：

```bash
# 初学阶段建议用这一条，环境写得最明确
julia --project=. bin/mottjain.jl plan

# 入口自身也会激活项目，所以这一条也能运行
julia bin/mottjain.jl plan
```

运行：

```bash
julia bin/mottjain.jl spectrum --config=config/my_run.toml
```

第一段加载链：

```text
bin/mottjain.jl
  → 激活 MottJainED 项目
  → using MottJainED
  → Julia 按包名打开 src/MottJainED.jl
  → src/MottJainED.jl 依次 include：
       Types.jl
       Model.jl
       Spectrum.jl
       CFT.jl
       Storage.jl
       Workflows.jl
       Conformal.jl
       Entanglement.jl
       CLI.jl
  → 所有类型和函数定义完成，但尚未开始物理计算
  → 调用 MottJainED.main(ARGS)
```

为什么 include 顺序重要：

- `Spectrum.jl` 使用 `Types.jl` 定义的 `ModelCache`；
- `CFT.jl` 使用 `Spectrum.jl` 的 `level_catalog`；
- `Workflows.jl` 再组合前面所有底层函数；
- `CLI.jl` 最后调用 Workflows。

---

## 八、最重要的数据类型是什么

这些类型都定义在 `src/Types.jl`。

### `Couplings`

保存一次 Hamiltonian 的八个系数：

```julia
Couplings(Uf, Uf0, U0, Vf, Vf0, V0, t, mu)
```

它只是参数容器，不包含矩阵。

### `SolverSettings`

保存数值方法参数：

- `k`：每个离散扇区求多少个低能态；
- `eig_tol`：ARPACK 容差；
- `energy_tol`：简并能量块容差；
- `quantum_tol`：(L^2,C_2) 接近整数的容差；
- `degeneracy_tol`：跨 `(Z,R)` 等能副本合并容差；
- `dense_cutoff`：稠密/稀疏对角化分界；
- `warm_start`：是否复用上一点的基态初始向量。

### `ModelParameters`

保存一个 `nm1` 下不会随耦合参数变化的物理结构：

- charge-1/charge-3 轨道数；
- QN 对称性；
- configuration space；
- (n_f,n_0,L^2,C_2,L_\pm)；
- Hamiltonian 每个独立项的 FuzzifiED `Terms`。

### `SectorCache`

保存一个 `(Z,R)` 扇区中可复用的数值对象：

- `Basis`；
- 固定部分矩阵 `h0`；
- 粒子数矩阵 `number_f`；
- `l2`、`c2` 矩阵；
- 上一次基态 `warm`。

### `ModelCache`

把同一个系统大小的四个 `SectorCache` 放在一起。

### `SpectrumState`

一个本征态的结果记录：

- `energy`；
- `l2`；
- `c2`；
- 所属 `(Z,R)`；
- 扇区内 rank；
- 可选本征向量和 Basis。

### `PhysicalLevel`

把不同 `(Z,R)` 中属于同一 SU(3) multiplet 的等能副本合并成一个物理能级，
同时记录 `multiplicity` 和原始成员。

### `CFTScore`

保存 tower 诊断结果：

- `valid`：所需能级是否齐全；
- `q`：tower relation 的 RMS 误差；
- `factor`：能量到 scaling dimension 的拟合比例；
- `delta_s`、`delta_o`；
- 每条关系的原始 gap、理论目标和失败原因。

---

## 九、用 spectrum 完整追踪一次计算

命令：

```bash
julia bin/mottjain.jl spectrum --config=config/my_run.toml
```

### 第 1 步：`bin/mottjain.jl`

文件：`bin/mottjain.jl`

作用：激活项目、加载包，把：

```julia
ARGS = ["spectrum", "--config=config/my_run.toml"]
```

交给：

```julia
MottJainED.main(ARGS)
```

### 第 2 步：`main`

文件：`src/CLI.jl`

函数：`main(args)`

内部步骤：

```text
_parse_cli(args)
  → command = "spectrum"
  → options["config"] = "config/my_run.toml"

load_config(config_path)
  → TOML.parsefile
  → 得到配置 Dict

_couplings(config)
  → Couplings(Uf, Uf0, ..., mu)

_solver(config)
  → SolverSettings(k, tolerances, ...)

_nm1(config)
  → 读取单尺寸 nm1

_range(config["spectrum"])
  → 生成需要计算的 mu 数组
```

然后 `main` 调用：

```julia
run_spectrum_scan(nm1, mus, couplings, settings; output=...)
```

### 第 3 步：建立物理模型

文件：`src/Workflows.jl`

函数：`run_spectrum_scan`

首先调用：

```julia
model = build_model(nm1=nm1)
```

`build_model` 位于 `src/Model.jl`，它计算：

```text
s = (nm1-1)/2
charge-1: nf1=3, no1=3*nm1
charge-3: nm0=3*nm1-2, nf0=1, no0=nm0
```

随后构造：

- 总电荷、(2L_z)、SU(3) Cartan 对称量；
- 离散 flavor permutation 和 y-rotation；
- configuration space `Confs`；
- charge-1 算符 (f_1,f_2,f_3)；
- charge-3 算符 (f_0)；
- (n_f,n_0,N_f,L^2,C_2,L_\pm)；
- Hamiltonian 八个独立 Terms：`Uf/Uf0/U0/Vf/Vf0/V0/t/mu`。

这一步没有选择具体 μ，也没有求本征态。

### 第 4 步：构造可复用矩阵缓存

文件：`src/Spectrum.jl`

函数：

```julia
cache = prepare_spectrum(model, couplings, settings)
```

它对四个 `(Z,R)`：

```text
(+1,+1), (+1,-1), (-1,+1), (-1,-1)
```

分别构造：

1. `Basis`；
2. 不含 μ 的 `H0`；
3. `Nf`；
4. `L2`；
5. `C2`。

因此以后每个 μ 只需要：

```text
H(mu) = H0 + mu*Nf
```

不再重建 Basis 和每个观测量。

`lower_sparse` 会把 FuzzifiED 的非排序 CSC-like 存储 canonicalize，避免稀疏加法
漏元素；`hermitian_opmat` 再恢复 Hermitian 标记。

### 第 5 步：对每个 μ 求谱

`run_spectrum_scan` 对 `mus` 循环，调用：

```julia
states = solve_spectrum(cache, mu)
```

文件：`src/Spectrum.jl`

每个 sector 内：

```text
_eigensystem
  → H0 + mu*Nf
  → 小矩阵：LinearAlgebra.eigen
  → 大矩阵：FuzzifiED.GetEigensystem/ARPACK
  → 可选 warm start
```

之后 `_resolve_quantum_numbers!`：

1. 找能量简并块；
2. 在块内重新对角化 (L^2)；
3. 同一个 (L^2) 子块再对角化 (C_2)；
4. 计算每个态的 `energy/l2/c2`。

四个 sector 的结果合并、按能量排序，形成 `Vector{SpectrumState}`。

### 第 6 步：计算 CFT tower score

文件：`src/CFT.jl`

调用：

```julia
score = cft_score(states; settings=settings)
```

先由 `level_catalog`：

- 把接近整数的 `l2/c2` 分类；
- 合并不同 `(Z,R)` 的等能副本；
- 得到按 `(L²,C₂)` 索引的 distinct physical levels。

再检查七条关系所需的能级是否齐全，拟合 `factor` 并计算 `q`、`delta_s`、
`delta_o`。

### 第 7 步：保存结果

文件：`src/Storage.jl` 和 `src/Workflows.jl`

```text
stable_id
  → 根据 nm1/mu/Hamiltonian/k 生成稳定 job_id

spectrum_dataframe
  → SpectrumState 转成长表 DataFrame

atomic_csv
  → 每个 μ 写 spectra/<job_id>.csv

append_csv
  → 完成后向 summary.csv 追加一行
```

如果中断后重跑，`completed_job_ids` 从 summary.csv 找到已成功 job，直接跳过。

### spectrum 最终输出

```text
output/<run_name>/spectrum/
├── run_metadata.toml
├── summary.csv
└── spectra/
    ├── <job_id-1>.csv
    ├── <job_id-2>.csv
    └── ...
```

---

## 十、当前所有功能和完整调用链

下面每条链中，箭头表示“前一个函数调用下一个函数”。公共的模型/谱底层不再
每次逐行重复，但都明确列出。

### A. `plan`：只预览，不计算

```text
bin/mottjain.jl
  → CLI.main
  → load_config
  → _plan
  → _range 统计 spectrum μ 点数
  → 统计 nm_values × scan_values
  → 打印计划
  → 不调用 build_model，不产生输出
```

### B. `spectrum`：固定参数扫描 μ 的完整低能谱

```text
CLI.main
  → load_config / _couplings / _solver / _range
  → Workflows.run_spectrum_scan
  → Model.build_model
  → Spectrum.prepare_spectrum
  → 对每个 μ：
       Spectrum.solve_spectrum
         → _eigensystem
         → _resolve_quantum_numbers!
       CFT.cft_score
         → Spectrum.level_catalog
       spectrum_dataframe
       Storage.atomic_csv / append_csv
  → summary.csv + 每个 μ 的 spectrum CSV
```

### C. `gap`：不同尺寸的 singlet gap–μ

```text
CLI.main
  → Workflows.run_gap_scan
  → 对每个 nm1：
       Model.build_model
       Spectrum.prepare_spectrum
       对每个 μ：
         Spectrum.solve_spectrum
         Spectrum.level_catalog
         取 ground energy
         取 (L²,C₂)=(0,0) 的第二个 distinct level
         gap = E_singlet - E_ground
         scaled_gap = gap*sqrt(nm1)
         Storage.append_csv
  → CairoMakie 画 singlet_gap.png
```

### D. `density`：基态粒子数与密度

```text
CLI.main
  → Workflows.run_density_scan
  → Model.build_model
  → Spectrum.prepare_spectrum
  → 对每个 μ：
       Spectrum.solve_spectrum(keep_vectors=true)
       取全局最低能态及其 Basis
       <Nf> = <ground|Nf|ground>
       N0 = (3*nm1-Nf)/3（固定总电荷约束）
       Storage.append_csv
  → density.csv + density.png
```

### E. `critical`：固定 U/V/t 优化 μ

```text
CLI.main
  → Workflows.run_critical_search
  → Model.build_model
  → Spectrum.prepare_spectrum（只做一次）
  → CFT.optimize_mu
       → coarse μ grid
       → 每个 μ：solve_spectrum → cft_score
       → 选择最佳粗网格相邻区间
       → Optim.Brent 一维精细优化
       → 最优 μ 再算一次完整 score
  → critical_point.csv
  → tower_residuals.csv
```

### F. `optimize`：同时优化多个 Hamiltonian 参数

```text
CLI.main
  → _optimization_bounds
  → Workflows.run_parameter_optimization
  → Model.build_model
  → _prepare_linear_family
       → 每个 sector 构造 fixed matrix
       → 每个自由参数各构造一个 derivative matrix
       → 构造 L2/C2
  → Optim.Fminbox(NelderMead)
       → 每次参数 evaluation：
           _solve_linear
             → fixed + Σ parameter_i*derivative_i
             → Spectrum.solve_spectrum
           CFT.cft_score
           立即 append evaluations.csv
  → best.csv
```

如果已有 `evaluations.csv`，会读取其中最小有效 `q` 的参数作为新起点；不会恢复
整个 Nelder–Mead simplex。

### G. `fss`：跨尺寸、跨耦合扫描并在每点优化 μ

```text
CLI.main
  → 构造 FSSSettings
  → Workflows.run_fss_scan
  → 对每个 nm1：
       Model.build_model（一次）
       对每个 scan_value：
         with_coupling 修改被扫描参数
         第一个值：Spectrum.prepare_spectrum
         后续值：Spectrum.retune_spectrum!
           （复用 Basis/L2/C2/Nf，只重建 H0）
         CFT.optimize_mu
           → 多次 solve_spectrum → cft_score
         Storage.append_csv
  → fss_results.csv
```

一行结果包含：`nm1`、`scan_value`、最优 `mu`、`q`、`factor`、`delta_s`、
`delta_o`、是否在 μ 边界、收敛状态和 Hamiltonian 全部系数。

### H. `fss-plot`：只读数据画图

```text
CLI.main
  → Workflows.plot_fss
  → CSV.read(fss_results.csv)
  → latest_rows 去除同 job 的旧记录
  → _valid_fss 只保留成功且 score_valid 的行
  → 每个 scan_value 画 y 对 nm1^(-1/2)
  → <y>_fss.png
```

不调用 `build_model`，不做本征求解，所以很快。

### I. `fss-fit`：联合有限尺寸拟合

```text
CLI.main
  → Workflows.fit_fss
  → 读取并过滤 fss_results.csv
  → 对给定 omega：
       建立线性 design matrix
       解析最小二乘求 Delta_inf 和每条曲线 amplitude
  → 在 log(omega) 粗网格扫描
  → Brent 精细优化 omega
  → <y>_fit.csv + <y>_fit.png
```

拟合模型：

```text
y(nm1,g) = Delta_inf + amplitude_g*(nm1^(-1/2))^omega
```

### J. `fss-all`

```text
run_fss_scan
  → plot_fss(y=delta_s)
  → fit_fss(y=delta_s)
```

只是把 G、H、I 顺序执行。

### K. `scaling`：一个参数点的 scaling-dimension spectrum

```text
CLI.main
  → Workflows.plot_scaling_dimensions
  → Model.build_model
  → Spectrum.prepare_spectrum
  → Spectrum.solve_spectrum(mu=hamiltonian.mu)
  → CFT.cft_score 得到 factor（或配置手动给 factor）
  → Delta_i = (E_i-E_ground)/factor
  → ell = (sqrt(1+4L2)-1)/2
  → scaling_nm<n>.csv + scaling_nm<n>.png
```

### L. `generator`：拟合共形生成元

```text
CLI.main
  → Workflows.run_generator_analysis（定义在 Conformal.jl）
  → Model.build_model
  → Spectrum.prepare_spectrum
  → Spectrum.solve_spectrum(keep_vectors=true, k=generator.k)
  → Conformal.build_conformal_store
       → Spectrum.level_catalog
       → 按 default_conformal_specs 选择 G/S/dS/T/J/... 波函数
  → Conformal.generator_candidates
       → 构造 18 个 microscopic L=1,m=0 Terms
  → Conformal.fit_generator(S,dS,candidates)
       → 每个 candidate 作用在 |S>
       → 组成 design matrix
       → 截断 SVD 最小二乘拟合 |dS>
       → 合成 Lambda Terms 并计算 fidelity
  → generator_coefficients.csv
  → generator_summary.csv
  → generator.jld2（保存向量/Terms 等复杂对象）
```

`project_angular_momentum` 和 `generator_overlap` 是后续手动研究其它 parent/descendant
overlap 的 API；默认命令当前主要完成 (S\to\partial S) 的 generator fit。

### M. `oes`：orbital entanglement spectrum

```text
CLI.main
  → Entanglement.run_orbital_entanglement
  → Model.build_model
  → _ground_state
       → prepare_spectrum
       → solve_spectrum(keep_vectors=true)
  → _orbital_cut 决定 A/B 轨道
  → _diagonal_sectors 枚举 QA/Lz/F3/F8 可行扇区
  → FuzzifiED.GetEntSpec
  → _entanglement_dataframe：lambda → xi=-log(lambda)
  → _save_entanglement：选择权重最大的 QA、移动 xi 最小值、统计低支
  → spectrum.csv / qa_weights.csv / counting.csv / summary.csv / spectrum.png
```

### N. `rses`：real-space entanglement spectrum

```text
CLI.main
  → Entanglement.run_realspace_entanglement
  → Model.build_model
  → _ground_state
  → _hemisphere_amplitudes
       → SpecialFunctions.beta_inc 计算每个轨道落在区域 A 的振幅
  → 构造总电荷/2Lz 的 A/B 扇区
  → FuzzifiED.GetEntSpec
  → _entanglement_dataframe
  → _save_entanglement
  → 与 OES 同类的 CSV/PNG 输出
```

---

## 十一、每个 `src` 文件的函数索引

### `Types.jl`

| 名字 | 意义 |
|---|---|
| `Couplings` | Hamiltonian 八个系数容器 |
| `with_coupling` | 复制 Couplings，只修改一个系数 |
| `SolverSettings` | 本征求解和分类容差 |
| `ModelParameters` | 一个尺寸下不随耦合变化的模型结构/Terms |
| `SectorCache` | 一个 `(Z,R)` 扇区的 Basis 和矩阵缓存 |
| `ModelCache` | 四个离散扇区缓存 |
| `SpectrumState` | 一个原始本征态记录 |
| `PhysicalLevel` | 合并对称副本后的 distinct level |
| `CFTScore` | tower score 结果 |
| `FSSSettings` | FSS 外层扫描设置 |

### `Model.jl`

| 函数 | 意义 |
|---|---|
| `pad_qn_diag/offdiag` | 把 charge-1 QN 扩展到 charge-1+charge-3 轨道空间 |
| `pad_term/pad_sphere_observable` | 给 charge-3 左侧偏移轨道编号 |
| `build_model` | 构造模型、对称性、观测量和八个 Hamiltonian Terms |
| `hamiltonian_terms` | 用 Couplings 线性组合 Terms |
| `lower_sparse` | canonicalize FuzzifiED 下三角存储 |
| `float_opmat` | 强制构造 Float64 OpMat |
| `hermitian_opmat` | 将 canonical 下三角矩阵恢复为 Hermitian OpMat |

### `Spectrum.jl`

| 函数 | 意义 |
|---|---|
| `prepare_spectrum` | 建四 sector 的 H0/Nf/L2/C2 缓存 |
| `retune_spectrum!` | 改 U/V/t 后只重建 H0 |
| `_eigensystem` | 稠密或 ARPACK 本征求解 |
| `_resolve_quantum_numbers!` | 简并子空间重新对角化 L2/C2 |
| `solve_spectrum` | 合并四 sector 的低能态 |
| `level_catalog` | 量子数分类并合并等能对称副本 |
| `spectrum_dataframe` | 转为 CSV 表格 |

### `CFT.jl`

| 函数 | 意义 |
|---|---|
| `cft_score` | 七条 tower relation 的 factor、q、Delta |
| `optimize_mu` | coarse grid + Brent 优化 μ |
| `score_dataframe` | 每条 tower residual 表格 |

### `Storage.jl`

| 函数 | 意义 |
|---|---|
| `ensure_output` | 创建输出目录 |
| `atomic_csv` | 临时文件写完后原子替换，避免半个 CSV |
| `append_csv` | 长任务每完成一点就追加 checkpoint |
| `completed_job_ids` | 找到已经成功的任务 |
| `stable_id` | 根据任务参数生成稳定短哈希 |
| `write_run_metadata` | 记录 Julia/FuzzifiED/git/线程信息 |
| `latest_rows` | 同 job 多条记录只保留最后一条 |

### `Workflows.jl`

| 函数 | 对应命令 |
|---|---|
| `run_spectrum_scan` | `spectrum` |
| `run_gap_scan` | `gap` |
| `run_density_scan` | `density` |
| `run_critical_search` | `critical` |
| `run_parameter_optimization` | `optimize` |
| `run_fss_scan` | `fss` |
| `plot_fss` | `fss-plot` |
| `fit_fss` | `fss-fit` |
| `plot_scaling_dimensions` | `scaling` |

### `Conformal.jl`

| 函数/类型 | 意义 |
|---|---|
| `ConformalState/Store` | 带名字的 CFT 候选态数据库 |
| `default_conformal_specs` | 默认 G/S/dS/T/J 等量子数和 rank |
| `generator_candidates` | 18 个 microscopic generator Terms |
| `fit_generator` | 截断 SVD 拟合 generator 系数 |
| `project_angular_momentum` | 多项式投影到指定 ell |
| `generator_overlap` | 计算生成向量落入指定目标态的权重 |
| `run_generator_analysis` | `generator` 完整任务 |

### `Entanglement.jl`

| 函数 | 意义 |
|---|---|
| `_ground_state` | 得到最低能波函数及 Basis |
| `_orbital_cut` | 定义 orbital A/B |
| `_diagonal_sectors` | 枚举纠缠扇区 |
| `_hemisphere_amplitudes` | RSES 每轨道区域 A 振幅 |
| `_save_entanglement` | 选 dominant QA、统计、保存和画图 |
| `run_orbital_entanglement` | `oes` |
| `run_realspace_entanglement` | `rses` |

### `CLI.jl`

CLI 是 **Command Line Interface（命令行接口）**。它不做 Hamiltonian 数学，只做：

1. 解析终端文字；
2. 读取 TOML；
3. 把 TOML 转成 `Couplings/SolverSettings`；
4. 选择并调用上表某个 workflow。

| 函数 | 意义 |
|---|---|
| `load_config` | TOML 文件 → Dict |
| `_couplings` | `[hamiltonian]` → Couplings |
| `_solver` | `[solver]` → SolverSettings |
| `_range` | `min/max/count` → 均匀数值数组 |
| `_parse_cli` | `ARGS` → command 和 options |
| `_output_root` | 确定输出目录 |
| `_plan` | 只打印任务数量 |
| `main` | command 的总分派器 |

---

## 十二、日常到底应该改什么

### 每次换物理参数

复制并修改：

```text
config/my_run.toml
```

不要改 `Project.toml`，不要改 `Manifest.toml`。

### 想改变 CFT score 定义

查看：

```text
src/CFT.jl → cft_score
```

### 想增加/删除 Hamiltonian 项

查看：

```text
src/Types.jl → Couplings
src/Model.jl → build_model / hamiltonian_terms
```

### 想改变 FSS 拟合公式

查看：

```text
src/Workflows.jl → fit_fss
```

### 想改变共形态的 rank 选择

查看：

```text
src/Conformal.jl → default_conformal_specs
```

### 想改变纠缠谱切分或统计

查看：

```text
src/Entanglement.jl
```

---

## 十三、最小使用流程

第一次在新机器：

```text
确认 MottJainED 与 FuzzifiED.jl 同级
  → julia scripts/setup.jl
```

一次新的物理实验：

```text
复制 config/default.toml 为 config/<实验名>.toml
  → 只修改该配置
  → julia bin/mottjain.jl plan --config=...
  → 确认任务数量和 k
  → 本地短任务直接运行，服务器长任务用 Slurm
  → 查看 output/<run_name>/<command>/
  → 画图/拟合命令直接读取已有 CSV，不重新对角化
```

如果只记住一句话：

> `config` 决定“算什么”，CLI 决定“选择哪个功能”，Workflows 组织步骤，
> Model/Spectrum/CFT 才做核心数学，Storage 保存结果。
