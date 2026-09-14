# 独立的 N=7 ED 实验层

**当前入口（2026-09-14）：** Vf0=0.45 已接受实测 `muc=0.152375` 并完成本地归档。
下一个 Vf0=0.65 点使用相同的有硬上限自动 scout/refine 流程，详见
[N7_AUTO_PIPELINE.md](N7_AUTO_PIPELINE.md)。下面手工 scout/refine 命令保留为已验证方法和
历史记录，不是当前要重复提交的步骤。

这套代码只存在于 `MottJainED-fast-ed` 工作树中。它没有修改
`FuzzifiED.jl`，也没有替换 `MottJainED` 原有的 `prepare_spectrum`、
`solve_spectrum` 或 FSS/optimization 工作流。

## 它解决什么问题

原流程会在一个进程中同时保留四个 `(Z,R)` sector，并串行求解。新流程把工作
拆成三步：

1. 固定非 `mu` 的 Hamiltonian 参数，逐 sector 构造 `H0`、`Nf`、`L2`、`C2`，
   每完成一个便原子写入独立 JLD2 文件并释放内存。
2. 每个 `(mu, sector)` 由一个独立进程加载和求解。四个 sector 可以作为 Slurm
   array 并行；中断后只需重跑缺失任务。
3. 四个 CSV 全部存在且身份一致后才合并，并调用原来的 `cft_score` 计算 `q`、
   `DeltaS`、`DeltaO` 和五条 relation residual。

因此它主要降低 N=7 峰值内存、把单点四个 sector 的 wall time 变成可并行任务，
并让不同 `mu` 或 `k` 复用矩阵。它不会改变底层稀疏本征值算法，所以单个 sector
自身仍然是昂贵计算。

缓存 ID 包含非 `mu` 参数、矩阵相关的 MottJainED 源码哈希、FuzzifiED 版本与
源码哈希。配置或源码不一致时不会误用旧矩阵。所有结果采用临时文件加原子改名，
避免中断留下伪装成完整结果的半个文件。

## 先做本地验证

```bash
julia --project=. scripts/fast_ed.jl prepare \
  --config=config/fast_ed/n5_retained_validation.toml

for sector in 1 2 3 4; do
  julia --project=. scripts/fast_ed.jl solve \
    --config=config/fast_ed/n5_retained_validation.toml \
    --mu-index=1 --sector-index="$sector"
done

julia --project=. scripts/fast_ed.jl collect \
  --config=config/fast_ed/n5_retained_validation.toml

julia --project=. scripts/fast_ed.jl compare \
  --config=config/fast_ed/n5_retained_validation.toml --mu-index=1
```

`compare` 只允许 N<=6，防止为了验证而意外重复一次昂贵的 N=7 计算。
把上面 profile 换成 `config/fast_ed/n6_retained_validation.toml` 可做同样的
N=6 检查。

## 已完成的 retained N=7 流程（历史记录）

当前新点的全新重跑入口见后文及 [N7_MU_SEARCH.md](N7_MU_SEARCH.md)；本节保留原有实测流程。

N=7 首次服务器实测的四个缓存合计为 25 GB，构造用时 14 分 40 秒，峰值内存
22.90 GB。建议至少留 40 GB 可用工作盘。最大 sector 的 k=20 pilot 用时
7 分 59 秒，峰值内存约 13.3 GB。`sdicnormal` 的 8 CPU 额度约为 61 GB，已有
充足余量；prepare 与 sector array 模板现均使用 8 CPU/8 线程。首次 prepare
曾保守申请 16 CPU，不能把这个历史额度误用作当前默认。

第一次建立共享缓存时不要立刻连锁提交全部20个求解任务：

```bash
mkdir -p slurm-logs
cache_job=$(sbatch --parsable slurm/fast_ed_prepare.sbatch)
cache_job=${cache_job%%;*}
```

缓存完成后先运行 `seff $cache_job`。然后只提交 task 8；它对应中心 mu 的
`(+,+)` sector，也是 N=6 中维数最大的 sector：

```bash
pilot_job=$(sbatch --parsable --array=8 slurm/fast_ed_sector_array.sbatch)
pilot_job=${pilot_job%%;*}
```

pilot 完成后运行 `seff ${pilot_job}_8`。两次 MaxRSS 都有安全余量后，再提交完整
scout；task 8 会识别已有 CSV 并直接复用：

```bash
array_job=$(sbatch --parsable --cpus-per-task=8 --array=0-19%4 \
  slurm/fast_ed_sector_array.sbatch)
array_job=${array_job%%;*}
```

确认 array 全部成功（失败的 index 重提后）再提交收集，避免原 array 中一次失败
使 `afterok` 永久阻塞：

```bash
collect_job=$(sbatch --parsable slurm/fast_ed_collect.sbatch)
collect_job=${collect_job%%;*}
```

如果 cache 或 pilot 因内存失败，可以用 `--cpus-per-task=20` 重提同一命令；
计算线程仍保持8，已完成的 sector cache/result 会自动复用。

默认 profile 是 `config/fast_ed/n7_retained_k20.toml`。中心
`mu=0.14625732779985` 是从已审计的 N=5、N=6 临界点做的一步线性外推；scout
实际计算中心左右共五点。`%4` 把同时运行的任务限制为四个，避免20个进程同时
读取几十 GB 的共享缓存。

收集任务会在 solver 结果目录
`output/fast_ed/runs/<run-name>/nmN_<cache-id>/kK_<settings-id>/` 生成：

- `scan_summary.csv`：五个 mu 的 q、factor、DeltaS、DeltaO；
- `best_summary.csv`：当前 q 最低的一行，是画 N=5,6,7 FSS 图需要的 N=7 输入；
- `best_relations.csv`：最佳点的五条 CFT relation；
- `best_spectrum.csv`：最佳点的完整低能态表；
- `collection_manifest.toml`：完整性、最佳 mu 和 cache/solver 身份。

把这五个小文件和 prepare/sector 作业的 `seff` 输出交回来即可；不需要传递
20--30 GB 的 JLD2 矩阵缓存。N=5、6 的基线结果已经在本地。

### N=7 精细 mu 扫描

首次五点 scout 的最低采样点为 `mu=0.14625732779985`、`q=0.0793799`，但
步长 `0.005` 对 DeltaS 太粗。中心三点的二次估计把最低点放在约
`mu=0.1451--0.1453`。独立 profile
`config/fast_ed/n7_retained_k20_refine.toml` 因此固定原 Hamiltonian，只计算
`0.14450:0.00025:0.14600` 的七个 mu。它复用已有 25 GB cache：

```bash
refine_job=$(sbatch --parsable --cpus-per-task=8 --array=0-27%4 \
  --export=ALL,CONFIG=config/fast_ed/n7_retained_k20_refine.toml \
  slurm/fast_ed_sector_array.sbatch)
refine_job=${refine_job%%;*}
```

确认28个 array task 全部成功后收集：

```bash
collect_job=$(sbatch --parsable \
  --export=ALL,CONFIG=config/fast_ed/n7_retained_k20_refine.toml \
  slurm/fast_ed_collect.sbatch)
collect_job=${collect_job%%;*}
```

交接文件位于上述 solver 结果目录，而不是 `<run-name>` 的直接下一级。可以用
`find output/fast_ed/runs/n7_retained_refine -name collection_manifest.toml` 定位。
2026-09-11 已完成精扫：实际最佳点为 `mu=0.14475`，当前五条 relation 所需最大
sector rank 仅为 9，因此当前 N=7 FSS 使用 `k=20`，不提交额外 `k=30` 检查。
如果需要增加 `mu`，只编辑 `mus=[...]`
并把 array 范围设为 `0:(4*mu点数-1)%4`；已有 `mu`/sector 会自动复用。

## Uf0/Vf0 局部 FSS：恢复 stage12 引导搜索，N=7 逐点串行

这组 FSS 只使用可信的 N=5、6、7，并分别做两条一维扫描：

- `Uf0 = [1.65, 1.834, 2.00]`，固定 `Vf0=0.55`；
- `Vf0 = [0.45, 0.55, 0.65]`，固定 `Uf0=1.834`。

其余参数保持 retained point。2026-09-11 用户要求恢复之前 stage 的正确、高效做法：
撤回 `205fbff` 中 N=5 直接承担 41 点宽网格及两组串行的配置。第一阶段直接运行
`scripts/two_size_tuning.jl`，搜索设置与 stage12 一致：N=3、4 宽网格及谷底精修，
N=5 从同参数 N=4 muc 续接，N=6 从 N=5 续接；局部窗口半宽 0.02、9 点、边界扩窗，
N=5 保留宽 Brent 对照，N=6 每点宽审计，`k=30` 和五项 q score 不变。
小尺寸只引导 muc，不进入 matching 或最终 FSS。

Uf0 三个点与 Vf0 两个端点作为两个独立 Slurm array task，各用 8 CPU/线程。
中心 `(Uf0,Vf0)=(1.834,0.55)` 只在 Uf0 任务计算一次，随后给两条曲线共用。
合计 5 个唯一 Hamiltonian 点、10 行 N=5/6 结果；含 guide 时为 20 行。

以下是 N=5/6 阶段的提交入口（2026-09-12 已下载完整结果并通过审计，无需重算）：

```bash
mkdir -p slurm-logs
sbatch slurm/fss_retained_local_n56.sbatch
```

首次运行的结果目录为（相同配置重提会续用同一目录，已完成参数/尺寸自动跳过）：

```text
output/two_size_tuning/n56_retained_local_uf0_01/
output/two_size_tuning/n56_retained_local_vf0_01/
```

每个任务计算后自动审计 `search/fss_optimize_results.csv`：Uf0 共 12 行，Vf0 共 8 行，
检查完整性、收敛、边界、有限结果、N=6 local/wide 一致性，并验证 N=5/6 确实使用
同参数的前一尺寸 muc 作为 continuation seed。审计失败以非零状态结束，保留输出。
下载两个完整小目录，检查 `two_size_matching.csv`、`search/fss_optimize_results.csv`
及 `search/fss_optimize_evaluations.csv` 后，才安排第一个新 N=7 点。
旧 `output/fast_ed/fss_n56_*` 的部分结果保留，但不会混入新流程。

**2026-09-13 已完成步骤：** Uf0=1.65、Vf0=0.55、V0=0.34 的七点精扫 548498
全部成功，28 CSV/TOML、560 态通过身份/选态检查，两种独立评分重算完全一致。
最低采样点 μ=0.14425、q=0.09688635 落在精扫左端；随后两点补扫已完成。
前三个精扫点估计谷底约0.144344，仅据此选择两个实测补点：
`config/fast_ed/n7_uf0_165_k20_followup.toml`，μ=`[0.14400,0.144375]`。
前者检查左侧回升，后者位于已有0.14425/0.14450之间，靠近估计谷底。

只提交 `--array=0-7%4`，每任务8CPU/线程；沿用548421的cache，不再prepare或删除。
结果保存到 `output/fast_ed/runs/n7_uf0_165_followup/`，不改底层求解器或五项q。
两个补点必须和已有scout/精扫一起在本地分析，不能对两点目录单独要求内部最小值。
本轮各sector耗时5m55s–8m44s，整个精扫从首个开始至最后结束约53m35s。
合并14个实测mu后接受 `mu=0.144375, q=0.09683425844993852`；左右近邻回升，
选态/身份/连续性审计通过。其小结果已在本地归档，下一点为Uf0=2.00自动pilot。
用户网络不稳定，命令必须拆成 pull 前、单独 pull、成功后提交三段。
完整精扫审计、资源记录与两点补算提交方式见 [N7_MU_SEARCH.md](N7_MU_SEARCH.md)。

N=7 阶段每次只处理一个新参数点。同一参数点收集后，可先运行只读审计：

```bash
julia --project=. scripts/fast_ed.jl audit \
  --config=POINT.toml --require-interior
```

确认小型 CSV/TOML 结果已经保存在 run 目录后，只有显式设置
`fast_ed.allow_cache_release=true` 的临时参数点 profile 才允许释放自己的缓存：

```bash
julia --project=. scripts/fast_ed.jl release-cache \
  --config=POINT.toml --require-interior --confirm-release
```

该命令再次验证 collection 完整、best μ 不在扫描边界、cache/result 身份一致，且
只删除 profile 精确指向的 `nmN_<cache-id>` 目录。retained 中心 profile 没有删除
许可，因此不会被这个串行清理流程误删。若审计失败，缓存保留，后续点不应开始。

## 目录

- 矩阵缓存：`output/fast_ed/cache/nmN_<cache-id>/`
- 分 sector 与合并结果：`output/fast_ed/runs/<run-name>/.../`
- 每个 `mu` 的最终文件：`merged_spectrum.csv`、`score_summary.csv`、
  `score_relations.csv`
- solver 结果目录的画图交接文件：`scan_summary.csv`、`best_summary.csv`、
  `best_relations.csv`、`best_spectrum.csv`、`collection_manifest.toml`

目前没有更改 FuzzifiED 底层。如果以后 profiling 证明需要改 JLL/Fortran 的
稀疏矩阵乘法或 eigensolver 接口，必须另建 FuzzifiED fork/构建目录并先取得许可。

## 当前自动N=7点（2026-09-14）

`config/fast_ed/n7_original_audit_auto.toml`计算stage1--12之前最终locator audit建议的原始
Hamiltonian：`(Uf,Uf0,U0,Vf,Vf0,V0,t)=(0.46,1.834,4.14,0,0.41,0.525,0.5)`。
N5/N6五项q定位结果分别为`mu=0.1220681019828/0.1218143040052`、
`q=0.11755568086489525/0.06517470973203797`。阻尼延拓中心为`0.1216874050164`，首次
五点scout为`0.1116874050164–0.1316874050164`，共20个独立sector task。

旧N5/N6宽扫描覆盖`mu=0.05778640450004–0.21057280900008`；两个尺寸均只有一个离散局部
q极小值，local refine与wide Brent的最优mu相差约`4e-5`，所以该N7点不设人为guard。
流程仍是一点一cache、每个sector任务8CPU/8线程、最多四个并发，完成后只下载最终tar.gz
和sha256。启动脚本在登录节点只检查空队列并提交一个1 CPU bootstrap；Julia加载/预编译、
旧cache核验与受控清理、流水线初始化都在bootstrap allocation内完成。

上一点Vf0=.65已经完成本地raw与归档SHA审计。其服务器cache只有在
`n7_vf0_065_cache_retirement.toml`实际计算、manifest与手输ID都严格对应
`acae800d795e294b`，且本地归档SHA严格等于
`e2f3e2f6545e8a7026ff402e8fb06271d3a8650e990751b5e6a27f08244f1d77`时才可受控释放。
retained中心cache继续保护。新原始点的`auto_release=false`，必须等其最终包下载并独立审计。

## 原始 audit 基线的 Uf0/Vf0 局部 N=5/6 与首个 N=7

原中心 `(Uf0,Vf0)=(1.834,0.41)` 的 N=5/6 已完成并通过审计，不重复计算。围绕它只新增
四个 Hamiltonian：`(1.65,0.41)`、`(2.00,0.41)`、`(1.834,0.30)`、
`(1.834,0.52)`；其余耦合固定为原始 audit 基线
`(Uf,U0,Vf,V0,t)=(0.46,4.14,0,0.525,0.5)`。

`slurm/fss_original_local_n56.sbatch` 把四点拆成四个独立 array task，最多四个并发。
每个 task 只运行 `scripts/two_size_tuning.jl`：N=3、4 做宽范围 muc 引导，N=5、6
续接并进入 matching；`k=30`、五项 q score 和 stage12 的 adaptive wide audit 保持不变。
任务在 N=6 审计和输出后停止，不会提交、生成或清理任何 N=7 cache。

```bash
mkdir -p slurm-logs
sbatch slurm/fss_original_local_n56.sbatch
```

首次运行对应四个 `output/two_size_tuning/n56_original_*_01/` 目录。四点均已下载并通过
repository CSV audit 与基于 raw gaps/targets 的独立复算；N=5/6 均未落在边界，wide challenger
收敛且与 local optimum 一致。联合原中心后的主要结果为：

| 点 | q(N=5) | q(N=6) | DeltaS drift | DeltaO drift |
|---|---:|---:|---:|---:|
| `Uf0=1.65` | 0.0999568 | 0.0625464 | -0.0402149 | -0.124394 |
| `Vf0=0.52` | 0.1084034 | 0.0725442 | -0.0350447 | -0.125767 |
| 原中心 | 0.1175557 | 0.0651747 | -0.0493959 | -0.138465 |
| `Uf0=2.00` | 0.1449348 | 0.0832512 | -0.0548342 | -0.146052 |
| `Vf0=0.30` | 0.1467821 | 0.0832203 | -0.0539804 | -0.142118 |

因此首个 N=7 只计算 `config/fast_ed/n7_original_uf0_165_auto.toml`，即
`(Uf0,Vf0,V0)=(1.65,0.41,0.525)`。它在两个尺寸的 q 都是五点中最低，且两种 drift
相对原中心也同时缩小；`Vf0=0.52` 保留为第二候选，等待首点 N=7 审计后再决定。

首点的阻尼延拓中心为 `mu=0.12139298275993`，初始五点 scout 为
`0.11139298275993–0.13139298275993`，`k=20`、cold start、20 个独立 sector task；后续
refine/followup 保持最多 16 个 mu 和 3 轮的有界规则。N=3/4 各只有一个引导 valley，且
N=5/6 local/wide 一致，所以不另加 guard。该 profile 的 `retire_before_start=[]` 且
`auto_release=false`：不会删除旧 stage8 中心、原始 audit 中心或任何其他 cache。

```bash
bash scripts/submit_fast_ed_pipeline.sh config/fast_ed/n7_original_uf0_165_auto.toml
```

正常完成后只需下载 `output/fast_ed/archives/n7_original_uf0_165_auto_final.tar.gz`
及对应 `.sha256`；不要在服务器运行 collect，也不必例行下载日志或 sacct。

首个局部点已完成并通过独立 raw-spectrum 审计。第二个 N=7 点使用
`config/fast_ed/n7_original_vf0_052_auto.toml`，Hamiltonian 为
`(Uf0,Vf0,V0)=(1.834,0.52,0.525)`。N=5/6 的 muc 为
`0.11307786788428/0.11270441406241`，阻尼延拓中心为 `0.112517687151475`；首轮仍为
中心及其 `±0.005,±0.01` 五点，共 20 个 sector task。N=3/4 各只有一个引导 valley，且
N=5/6 local/wide 一致，因此不另设 guard。

启动时只会通过 `n7_original_uf0_165_cache_retirement.toml` 受控释放刚刚完成、已经下载并
审计的 `nm7_3bfe378c1c294c11`。retirement profile 同时核对该 cache ID、manifest身份与
归档SHA `0a0479121496c7e05e4011f3e1835a273a4495af420e713b34c112396a56a185`；两个中心
cache及其他Hamiltonian不在清理列表中。新点保持 `auto_release=false`。

```bash
bash scripts/submit_fast_ed_pipeline.sh config/fast_ed/n7_original_vf0_052_auto.toml
```

正常完成后下载 `output/fast_ed/archives/n7_original_vf0_052_auto_final.tar.gz` 及对应
`.sha256`，再在本地独立审计并合并 N=5/6/7。

## 已完成的本地等价性检查

- N=5、k=20：80 态键完全相同，最大能量误差 `3.24e-14`，最大 L2/C2
  误差 `1.38e-11` / `1.54e-12`。
- N=6、k=20：80 态键完全相同，最大能量误差 `9.68e-14`，最大 L2/C2
  误差 `1.47e-11` / `8.52e-13`。
- 两个尺寸的 `q`、`DeltaS`、`DeltaO` 与旧路径差异都在约 `1e-13`。
- N=6 单个最大 sector 用 8 线程约 5.4 秒；旧路径四 sector 串行求解约
  21.4 秒。这说明四个各占 8 核的 Slurm 任务有接近四路 wall-time 并行的空间。
  不要在只有 8 核的本机同时跑四个 8 核进程。
