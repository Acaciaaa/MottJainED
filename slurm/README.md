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
| FSS 数据、图和拟合 | `fss.sbatch`（数组任务分别扫描 Uf0/Vf0） |
| Hamiltonian 参数优化 | `optimize.sbatch` |
| scaling dimension 图 | `scaling.sbatch` |
| generator ED 与 Lambda 拟合 | `generator.sbatch` |
| tower overlap | `tower.sbatch` |
| orbital entanglement spectrum | `oes.sbatch` |
| real-space entanglement spectrum | `rses.sbatch` |

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

因此这些脚本故意不写 `#SBATCH --mem`。需要更多内存时增加 CPU 申请数；设置
显式 `--mem` 反而会绕开这个分区的默认比例。

申请 CPU 的数量和实际计算线程已经分开。例如：

```bash
#SBATCH --cpus-per-task=56  # 为作业取得约 428 GiB 内存
threads=16                  # FuzzifiED 实际先使用 16 线程
```

申请56个 CPU 不代表一定要开56线程。线程过多可能受内存带宽、同步开销影响，
不一定更快；当前工作流也不会并行启动多次 ED。建议先用 8–16 个计算线程，
根据计时再尝试增加，但 `threads` 不可大于申请的 CPU 数。

增大 `nm1`、`k` 或保存更多本征向量时，应修改对应文件自己的 CPU/时间，
不会影响其他功能。当前集群规则以后若改变，以服务器的 `scontrol show partition
sdicnormal` 为准。

`plan`、`generator-register`、`fss-plot` 和 `fss-fit` 不做昂贵 ED，通常直接在
登录节点运行，不需要为它们申请专用作业。

每个作业会生成 `slurm-功能名-jobid.out`。程序自己的详细进度和结果仍写入
`output/<功能>/...`；Slurm 文件主要保留启动信息、warning 和错误。
