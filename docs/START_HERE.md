# 从这里开始：MottJainED 文件结构、概念和完整工作流程

这份说明假设你此前没有使用过 Julia 项目、TOML、命令行入口或软件包结构。
阅读完以后，你应该能够回答：

1. 每个文件为什么存在；
2. 哪些文件需要自己修改，哪些不要动；
3. 一条命令怎样一步步变成 Hamiltonian、本征求解和输出文件；
4. 当前有哪些功能，每个功能经过哪些函数；
5. 想检查或修改某一步时应打开哪个文件。

需要直接复制运行命令时，请看 `docs/USER_GUIDE.md` 第 4 节“命令”：那里集中列出
目前所有功能的案例、可改参数和输出位置；程序本身不再增加逐命令的 `--help` 层。

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
| `Manifest-v1.11.toml` 等 | 对应 Julia minor 版本的精确依赖 | 不手动改 |
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
FuzzifiED = {url = "https://github.com/FuzzifiED/FuzzifiED.jl.git", rev = "29a0cc9e06bcb5b30d3cf9f6db6416917f8a573f"}
```

这里把 FuzzifiED 固定到已经在本地与服务器核对过的 Git commit。
`Pkg.instantiate()` 会自动下载这份源码，因此服务器不需要另放一个相邻的
`FuzzifiED.jl` 文件夹，也不会意外跟随 `main` 的后续变化。

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

## 四、版本专用 Manifest 做什么

`Manifest-v1.11.toml` 第一行已经写明：

```toml
# This file is machine-generated - editing it directly is not advised
```

它是 Julia 自动生成的“精确软件清单”。

`Project.toml` 可能只说：

```text
需要 CairoMakie 0.15
```

但 CairoMakie 又依赖 Makie、Colors、GeometryBasics 等许多包。
Manifest 会记录整个依赖树的：

- 包名；
- 精确版本；
- UUID；
- 下载内容哈希；
- 包之间的依赖关系；
- FuzzifiED 的精确源码哈希和 Git revision；
- Julia 版本。

这就是为什么它很长。你不需要读懂它，也不要手工加入注释，因为下一次
`Pkg.resolve()` 可能重写它。

服务器执行：

```bash
julia --project=. scripts/setup.jl
```

时，Julia会优先读取与自己版本匹配的文件。例如本地Julia 1.11读取
`Manifest-v1.11.toml`；服务器Julia 1.12第一次运行时按`Project.toml`解析，
`setup.jl`随后保存为`Manifest-v1.12.toml`。不同minor版本不能共用一个通用
Manifest，因为stdlib及其JLL依赖也可能改变。

可以把二者类比为：

| 文件 | 类比 |
|---|---|
| `Project.toml` | 我需要“面粉、鸡蛋、牛奶” |
| `Manifest-v1.11.toml` | Julia 1.11所用面粉品牌/批次及所有供应链细节 |
| `Manifest-v1.12.toml` | Julia 1.12自己的对应供应链细节 |

---

## 五、整个文件夹分别是什么

```text
MottJainED/
├── Project.toml             Julia 直接依赖与项目身份
├── Manifest-v1.11.toml      本地 Julia 1.11 的精确依赖锁定
├── README.md                项目首页和文档入口
├── config/
│   ├── default.toml         默认物理/数值参数模板
│   ├── generator_points.csv 全局 generator 候选参数点表
│   ├── critical_profiles/   每个 critical 案例自己的 score
│   ├── fss_profiles/        每个 FSS 案例自己的 score
│   ├── optimization_profiles/  free/values/bounds/score 都在案例模板内
│   └── generator/
│       ├── templates/       新 point 的 fit/tower 模板
│       └── <point_id>/      该 Hamiltonian 独有的 generator_fit.toml/tower.toml
├── bin/
│   └── mottjain.jl          终端命令的短入口
├── scripts/
│   └── setup.jl             新机器第一次安装依赖
├── slurm/                   每个重计算功能自己的服务器作业文件
│   ├── fss.sbatch           FSS 的资源和固定命令
│   ├── optimize.sbatch      参数优化的资源和固定命令
│   ├── generator.sbatch     generator ED/Lambda 的资源和固定命令
│   └── ...                  spectrum/gap/density/critical/tower/ES 等
├── src/
│   ├── MottJainED.jl        包总入口，按顺序加载其它 src 文件
│   ├── Types.jl             所有核心数据类型
│   ├── Model.jl             物理模型与 Hamiltonian Terms
│   ├── Spectrum.jl          basis/matrix 缓存、本征求解、量子数分类
│   ├── CFT.jl               可自由组合的 CFT relation 池与固定 μ 网格选择
│   ├── Storage.jl           CSV、metadata、job ID、断点续跑工具
│   ├── Workflows.jl         spectrum/gap/density/critical/FSS 等完整任务
│   ├── Conformal.jl         共形生成元候选、拟合与 overlap
│   ├── GeneratorWorkflow.jl 全局参数表、ED 快照与 tower 后处理
│   ├── Entanglement.jl      OES 和 RSES
│   └── CLI.jl               把命令名字分派给 Workflows
├── test/
│   └── runtests.jl          自动测试
└── docs/
    ├── START_HERE.md        当前这份概念和代码地图
    └── USER_GUIDE.md        参数与命令参考
```

---

## 六、`scripts` 与 `slurm` 做什么

### `scripts/setup.jl`：只负责安装环境

它的调用链是：

```text
运行 setup.jl
  → Pkg.activate(MottJainED 根目录)
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

### `slurm/*.sbatch`：一个计算功能一个服务器作业

这些是 Bash/Slurm 脚本，不是 Julia 文件。旧的单一 `run_slurm.sh` 已被拆开，
所以切换功能不再需要注释/取消注释 Julia 语句。例如：

```bash
sbatch slurm/fss.sbatch
sbatch slurm/optimize.sbatch
sbatch slurm/generator.sbatch
```

每个文件前面的：

```bash
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
```

向 Slurm 申请：

- 16 CPU；
- `sdicnormal` 自动随 CPU 配给约 122 GiB 内存；
- 最长 48 小时。

当前 `sdicnormal` 是每 CPU 约 7824 MB 内存，节点为64 CPU、约501 GiB；因此专用
脚本不再写 `--mem`。例如申请56 CPU大约得到428 GiB，但可以在用户配置区保留
`threads=16`，只让 FuzzifiED 使用16线程。申请资源数和实际计算线程不是一回事。

同一个文件的“用户配置区”保存它自己的 profile、point 等少量选择。例如
`fss.sbatch` 中：

```bash
config="config/my_run.toml"
profile="config/fss_profiles/fss5.toml"
method="both"
threads=16
```

脚本末尾已经固定调用对应的 Julia 功能，一般不需要改。不同 `nm1` 所需内存可能
差很多，因此各文件顶部的 CPU/时间只是起点；调整 `fss.sbatch` 不会改变
`generator.sbatch`。完整文件表见 `slurm/README.md`。

如果服务器不用 Slurm，就不使用这个目录。

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
| `--override=...` | 读取基础配置后，只覆盖小模板中明确写出的字段 |
| `--nm1=7 --k=12` | 临时覆盖常改的系统大小和当前任务的 k |

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
- `CFT.jl` 使用 `Spectrum.jl` 生成的原始 `SpectrumState` 和量子数分类；
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
- `definition`：`critical5`、`fss7` 或 `optimization8`；
- `terms`：这次实际选择的 relation 名称列表；
- `metric` / `objective`：当前工作流实际最小化 `q` 还是 `cost`；
- `q`：tower relation 的 RMS 误差；
- `cost`：旧 FSS/optimization 的方向夹角目标；
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

config["hamiltonian"]["mu"]
  → spectrum 只计算这个单点

config["spectrum"]
  → 读取 k、l2_values、c2_values、levels_per_block、factor
```

然后 `main` 调用：

```julia
run_spectrum(nm1, couplings, settings; output=...)
```

### 第 3 步：建立物理模型

文件：`src/Workflows.jl`

函数：`run_spectrum`

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

这一步尚未求本征态。

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

在配置的单个 μ 上只需要：

```text
H(mu) = H0 + mu*Nf
```

不再重建 Basis 和每个观测量。

`lower_sparse` 会把 FuzzifiED 的非排序 CSC-like 存储 canonicalize，避免稀疏加法
漏元素；`hermitian_opmat` 再恢复 Hermitian 标记。

### 第 5 步：在单个 μ 求谱

`run_spectrum` 调用：

```julia
states = solve_spectrum(cache, couplings.mu)
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

### 第 6 步：合并物理能级、筛选并保存结果

文件：`src/Storage.jl` 和 `src/Workflows.jl`

```text
level_catalog
  → 在容差内把 L²/C₂ 认成整数标签
  → 合并不同 (Z,R) sector 中同量子数、同能量的副本
  → 不合并能量不同的能级

spectrum_level_table
  → 对配置的每个 (L²,C₂) 组合取最低若干能级，每个组合成为一列
  → 计算 rescaled_energy=(E-E0)/factor
  → 默认列序为 (0,0),(2,0),(6,0),(0,3),(2,3),(6,3)

atomic_csv
  → 写 spectrum.csv
```

`spectrum` 不套用 CFT tower 标准，也不输出任何算符身份。本征向量只在内存中
用于识别量子数，不写 JLD2。若某个组合不足配置数量，程序会提示增大
`[spectrum].k`。

如果 `spectrum.csv` 已完整存在，重跑会直接复用；显式加入 `--force` 才覆盖。

### spectrum 最终输出

```text
output/spectrum/<可选案例名_01>/
├── case_identity.toml
├── resolved_config.toml
├── run_metadata.toml
└── spectrum.csv
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
  → 展示 spectrum 单点 μ、筛选量子数、数量和 k
  → 统计 nm_values × scan_values
  → 打印计划
  → 不调用 build_model，不产生输出
```

### B. `spectrum`：单个 μ 的筛选 rescaled 物理能级

```text
CLI.main
  → load_config / _couplings / _solver
  → 从 [spectrum] 读取 k、L²/C₂ 列表、每组数量和 factor
  → Workflows.run_spectrum
  → Model.build_model
  → Spectrum.prepare_spectrum
  → Spectrum.solve_spectrum(hamiltonian.mu)
      → _eigensystem
      → _resolve_quantum_numbers!
  → level_catalog 合并严格同能的离散-sector副本
  → spectrum_level_table 筛选并 rescale
  → Storage.atomic_csv
  → spectrum.csv
```

### C. `gap`：不同尺寸的 scalar gap 与 J gap–μ

```text
CLI.main
  → 读取 [gap] 的 k（默认 5，与旧 ES_mu.jl 一致）
  → Workflows.run_gap_scan
  → 对每个 nm1：
       Model.build_model
       Spectrum.prepare_spectrum
       对每个 μ：
         Spectrum.solve_spectrum
         取 ground energy
         scalar_gap = 原始 (L²,C₂)=(0,0) 列表第二项 - E_ground
         j_gap = 最低 (L²,C₂)=(2,3) 态 - E_ground
         分别乘 sqrt(nm1)
         Storage.append_csv
  → gap_results.csv
  → CairoMakie 画 scalar_gap.png 和 j_gap.png
```

两种 gap 来自同一次求谱，所以增加 J gap 不会让每个 `(nm1,μ)`
重复对角化。两张图都是 `650×650`、正方形坐标区、浅色网格，纵轴固定为
`0–1.0`。

### D. `density`：基态粒子数与密度

```text
CLI.main
  → 读取 [density].k（默认 3，与旧 particle_density_mu.jl 一致）
  → Workflows.run_density_scan
  → Model.build_model
  → Spectrum.prepare_spectrum
  → 对每个 μ：
       Spectrum.solve_ground_state
       只比较四个 (Z,R) sector 的最低态，不做 L²/C₂ 分类
       取全局最低能态及其 Basis
       <Nf> = <ground|Nf|ground>
       N0 = (3*nm1-Nf)/3（固定总电荷约束）
       Storage.append_csv
  → density.csv + density.png
```

`density.png` 沿用旧图的 `650×650` 正方形布局，纵轴范围为 `0–3`，
并保留 `0,1,2` 的灰色虚线参考线。

### E. `critical`：在固定 μ 网格上选最小 score

```text
CLI.main
  → Workflows.run_critical_search
  → Model.build_model
  → Spectrum.prepare_spectrum（只做一次）
  → CFT.scan_mu
       → 只枚举配置的 mu_min:mu_max，mu_count 个点
       → 每个 μ：solve_spectrum → TOML 选定的 cft_score
       → 选 objective 最小的网格点，不调用 Optim
  → critical_scan.csv
  → critical_point.csv
  → tower_residuals.csv
```

`score_terms` 不放在公共 my_run 中；复制 `config/critical_profiles/critical5.toml`
并在子配置里增删 relation，就能让 critical 使用任意组合。

### F. `optimize`：优化一个或多个 Hamiltonian 参数

```text
CLI.main
  → _optimization_couplings
       → [optimization.values] 覆盖本次优化的初值/固定值
  → _optimization_bounds
  → Workflows.run_parameter_optimization
  → Model.build_model
  → _prepare_linear_family
       → 每个 sector 构造 fixed matrix
       → 每个自由参数各构造一个 derivative matrix
       → 构造 L2/C2
  → 一个 free：Optim.Brent
    多个 free：Optim.NelderMead（越界点返回 penalty）
       → 每次参数 evaluation：
           _solve_linear
             → fixed + Σ parameter_i*derivative_i
             → Spectrum.solve_spectrum
           CFT.cft_score
           立即 append evaluations.csv
  → best.csv
```

`free` 决定哪些参数变化；`[optimization.values]` 对自由参数表示初值、对其余
参数表示固定值。默认严格沿用旧 `optimization.jl` 当前启用的“超级大满贯”方案。
如果已有兼容的 `evaluations.csv`，会读取其中最小有效 `objective` 的参数作为
新起点；不会恢复整个 Nelder–Mead simplex。

### G. `fss`：跨尺寸、跨耦合，并用两种方法找 μc

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
         method=grid：
           CFT.scan_mu
             → 只计算 mu_count 个网格点 → cft_score
             → 选网格中 objective 最小点
         method=optimize：
           optimize_fss_mu_with_score
             → N=3 独立宽搜，不约束后续 size
             → N=4 宽搜，作为固定 scan_value 的 anchor
             → N=5/6 围绕前一 size 的 muc 做局部细网格
             → 窗口边界最低时自动扩窗
             → 与全范围 Brent 候选比较实际 score
         两种方法分别 Storage.append_csv
  → fss_grid_results.csv
  → fss_optimize_results.csv
  → fss_optimize_evaluations.csv（完整 q(μ) 求值轨迹）
```

两种方法默认都运行。也可以在 TOML 写 `methods=["grid"]` /
`methods=["optimize"]`，或在命令末尾临时加 `--method=grid`、
`--method=optimize`。一行结果包含：`nm1`、`scan_value`、`muc`、`objective`、
`q`、`cost`、`factor`、`delta_s`、`delta_o`、是否在 μ 边界和 Hamiltonian
全部系数；optimize 结果还记录 continuation 中心、实际窗口、扩窗次数、宽 Brent
候选和最终最低点的来源。

### H. `fss-plot`：只读数据画图

```text
CLI.main
  → Workflows.plot_fss
  → CSV.read(fss_grid_results.csv 或 fss_optimize_results.csv)
  → latest_rows 去除同 job 的旧记录
  → _valid_fss 只保留成功且 score_valid 的行
  → 每个 scan_value 画 y 对 nm1^(-1/2)
  → <y>_<method>_fss.png
```

不调用 `build_model`，不做本征求解，所以很快。

### I. `fss-fit`：联合有限尺寸拟合

```text
CLI.main
  → Workflows.fit_fss
  → 读取并过滤某一种 method 的 FSS CSV
  → 对给定 omega：
       建立线性 design matrix
       解析最小二乘求 Delta_inf 和每条曲线 amplitude
  → 只在 81 个固定 log(omega) 网格点中选残差最小点
  → 不调用 Optim
  → <y>_fit.csv + <y>_fit.png
```

拟合模型：

```text
y(nm1,g) = Delta_inf + amplitude_g*(nm1^(-1/2))^omega
```

### J. `fss-all`

```text
run_fss_scan
  → 对每个启用的 method：plot_fss(y=delta_s)
  → 对每个启用的 method：fit_fss(y=delta_s)
```

只是把 G、H、I 顺序执行。

### K. `scaling`：一个参数点的 scaling-dimension spectrum

```text
CLI.main
  → Workflows.plot_scaling_dimensions
  → Model.build_model
  → Spectrum.prepare_spectrum
  → Spectrum.solve_spectrum(mu=hamiltonian.mu)
  → 按旧 scaling_dimension.jl 的五条能隙关系得到 factor（或配置手动给 factor）
  → Delta_i = (E_i-E_ground)/factor
  → ell = (sqrt(1+4L2)-1)/2
  → scaling_nm<n>.csv + scaling_nm<n>.png
```

### L. `generator` 与 `tower`：先固定 ED，再反复试 tower

```text
CLI.main
  → GeneratorWorkflow.load_generator_point
       → 只按 --point 从 config/generator_points.csv 取一行
       → 得到 nm1、Uf/Uf0/U0/Vf/Vf0/V0/t、muc、factor
  → GeneratorWorkflow.ensure_generator_snapshot
  → Model.build_model
  → Spectrum.prepare_spectrum
  → Spectrum.solve_spectrum(keep_vectors=true, k=generator.k)
  → 可选求旧 for_generator_special 的 adjoint weight sector
  → 原子保存 ed_snapshot.jld2
  → 保存两张 physical_levels_*.csv 和一份可读元数据（原始向量只存一份）
  → GeneratorWorkflow.run_generator_fit
       → 读取 config/generator/<point_id>/generator_fit.toml 中固定的 S/dS
  → Conformal.generator_candidates
       → 构造 18 个 microscopic L=1,m=0 Terms
  → Conformal.fit_generator(S,dS,candidates)
       → 每个 candidate 作用在 |S>
       → 组成 design matrix
       → 截断 SVD 最小二乘拟合 |dS>
       → 合成 Lambda Terms 并计算 fidelity
  → output/generator/<point_id>/generator/generator_fit.jld2（固定保存 Lambda）
  → tower 命令调用 GeneratorWorkflow.run_tower_analysis
       → 只读取固定 Lambda 和 config/generator/<point_id>/tower.toml
       → 按 family/L2/C2/physical rank/member 选择可变 tower 态
  → 对每个 [[overlaps]] 做角动量投影和目标子空间 overlap
  → same_angular 模式做 L- → Lambda_z → L+
  → selected_states.csv + tower_overlaps.csv + analysis_metadata.toml
  → 终端用旧脚本的固定格子打印 Input、l'、Target(dE/f)、overlap 与 Total
  → 仅 save_generated_vectors=true 时额外保存 tower_analysis.jld2
```

推荐用 `generator --point=ID` 一次固定 ED 和 Lambda，再用 `tower --point=ID` 反复
修改其它态认定。tower 配置变化只产生 `tower_01/tower_02` 连续目录，不重复
ED/拟合。可见数据目录始终是 point_id；隐藏签名负责检查 Hamiltonian、k/容差和
代码版本是否真的一致。

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

### `CFT.jl`

| 函数 | 意义 |
|---|---|
| `cft_score` | 从候选 relation 池任意组合，计算 factor、q、cost、Delta |
| `scan_mu` | 只计算用户指定的 μ 网格并选 objective 最小点 |
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
| `start_task_logging` | Info 追加到 `run.log`，终端只保留 Warn/Error |
| `latest_rows` | 同 job 多条记录只保留最后一条 |

### `Workflows.jl`

| 函数 | 对应命令 |
|---|---|
| `run_spectrum` | `spectrum` |
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
| `run_generator_analysis` | 旧的单点直接拟合 API；CLI 已改用快照工作流 |

### `GeneratorWorkflow.jl`

| 函数/类型 | 意义 |
|---|---|
| `GeneratorPoint` | 全局 CSV 中一行经过校验的参数点 |
| `GeneratorEDSnapshot` | 按 sector 打包的 basis、能量、量子数和本征向量 |
| `register_optimization_point` | 把选中的 optimization `best.csv` 登记到全局表 |
| `ensure_generator_snapshot` | 建立或按精确身份复用 ED 快照 |
| `locate_generator_snapshot` | tower 严格定位已有快照，不允许偷偷重算 |
| `run_generator_fit` | 在 ED 快照上拟合一次并固定保存 Lambda |
| `load_generator_fit` | 读取与快照绑定的固定 Lambda |
| `run_tower_analysis` | 只用固定 Lambda，按独立 TOML 选态和计算多层 overlap |

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

不要改 `Project.toml`，不要手工改任何 `Manifest-v*.toml`。

### 想切换 optimization 固定/自由参数

只改三处：

```toml
[optimization]
free = ["Uf0", "V0", "mu"]

[optimization.values]
Uf = 0.5       # 不在 free：固定值
Uf0 = 3.0      # 在 free：初值
V0 = 1.0       # 在 free：初值
mu = 0.05      # 在 free：初值

[optimization.bounds]
Uf0 = [0.0, 20.0]
V0 = [-20.0, 20.0]
mu = [-100.0, 100.0]
```

不再通过注释 `make_tms_hmt` 调用来切换方案。完整旧方案对照见
`docs/USER_GUIDE.md` 的“单参数或多参数优化”。

### 想改变 CFT score 定义

日常选择不需要改代码，复制相应的 critical/FSS/optimization 子 profile，再在
子 TOML 中增删列表元素：

```toml
score_terms = ["ds_s", "j", "curlj", "dj_rank1", "t_rank2"]
score_metric = "q"        # 或 cost
```

完整候选名称和对应 `(L²,C₂,rank)` 见 `docs/USER_GUIDE.md`。只有要新增候选关系
或改变已有候选的态定义时才查看：

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
clone MottJainED 并进入项目目录
  → julia --project=. scripts/setup.jl
  → 首次解析/预编译可能持续十几分钟
```

一次新的物理实验：

```text
复制 config/default.toml 为 config/<实验名>.toml
  → 只修改该配置
  → julia --project=. bin/mottjain.jl plan --config=...
  → 确认任务数量和 k
  → 本地短任务直接运行，服务器长任务用 Slurm
  → 查看 output/<command>/<可选案例名_01>/
  → 画图/拟合命令直接读取已有 CSV，不重新对角化
```

如果只记住一句话：

> `config` 决定“算什么”，CLI 决定“选择哪个功能”，Workflows 组织步骤，
> Model/Spectrum/CFT 才做核心数学，Storage 保存结果。
