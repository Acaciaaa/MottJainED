# 专用 Slurm 作业

这里采用“一个重计算功能一个 `.sbatch` 文件”的结构。平时不需要拼很长的
`sbatch ... Julia ...` 命令，只需打开对应文件修改两处：

1. 顶部 `#SBATCH` 中的 CPU 和最长时间；
2. `用户配置区`中的 config、profile、point、`threads` 或少量命令参数。

然后在 MottJainED 根目录提交，例如：

```bash
sbatch slurm/fss.sbatch
sbatch slurm/optimize.sbatch
sbatch slurm/generator.sbatch
```

文件对应关系：

| 功能 | 文件 |
|---|---|
| spectrum | `spectrum.sbatch` |
| gap | `gap.sbatch` |
| density | `density.sbatch` |
| critical | `critical.sbatch` |
| FSS 数据和 DeltaS/DeltaO 图 | `fss.sbatch`（数组任务分别扫描 Uf0/Vf0） |
| Hamiltonian 参数优化 | `optimize.sbatch` |
| SO(3)lver N=6 六项嵌套参数优化 | `so3lver_n6_nested_optimize.sbatch` |
| SO(3)lver N=7 同协议六项嵌套比较 | `so3lver_n7_nested_optimize.sbatch` |
| scaling dimension 图 | `scaling.sbatch` |
| generator ED 与 Lambda 拟合 | `generator.sbatch` |
| tower overlap | `tower.sbatch` |
| orbital entanglement spectrum | `oes.sbatch` |
| real-space entanglement spectrum | `rses.sbatch` |
| N=7 FastED k=40 generator/tower | `fast_ed_vector_sector_array.sbatch` + `fast_ed_tower_collect.sbatch`（由 `scripts/submit_n7_original_tower_k40.sh` 提交） |

## sdicnormal 的 CPU 与内存

当前服务器查询结果是 `DefMemPerCPU=7824 MB`：不写 `--mem` 时，Slurm 按
`--cpus-per-task` 自动分配约 7.824 GB/CPU。当前节点有 64 CPU 和约 513024 MB
总内存。常用换算约为：

| 申请 CPU | Slurm 默认内存（约） |
|---:|---:|
| 8 | 61 GiB |
| 16 | 122 GiB |
| 32 | 245 GiB |
| 56 | 428 GiB |
| 64 | 489 GiB |

这个表只说明未设置内存时的隐式额度，不是推荐申请量。已有同类任务的`seff`数据时，
脚本应按实测`MaxRSS`加合理余量显式设置总内存；CPU数则按实际并行工作量决定，不能为了
取得内存而申请不会使用的CPU。只有尚无可信测量、或任务确实需要默认比例时才暂用分区默认。

若某个高内存任务的内存需求远高于CPU需求，应显式申请总内存，而不是借CPU取得内存。例如：

```bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=428G          # 仅在实测确实需要时
threads=16
```

线程过多可能受内存带宽、同步开销影响，不一定更快；当前工作流也不会并行启动多次ED。
CPU和内存都应分别由同类任务的计时与`MaxRSS`决定，且`threads`不可大于申请的CPU数。

增大 `nm1`、`k` 或保存更多本征向量时，应修改对应文件自己的 CPU/时间，
不会影响其他功能。当前集群规则以后若改变，以服务器的 `scontrol show partition
sdicnormal` 为准。

`plan`、`generator-register`、`fss-plot` 和 `fss-fit` 不做昂贵 ED，通常直接在
登录节点运行，不需要为它们申请专用作业。

每个作业会生成 `slurm-功能名-jobid.out`。程序自己的详细进度和结果仍写入
`output/<功能>/...`；Slurm 文件主要保留启动信息、warning 和错误。
