# 专用 Slurm 作业

这里采用“一个重计算功能一个 `.sbatch` 文件”的结构。平时不需要拼很长的
`sbatch ... Julia ...` 命令，只需打开对应文件修改两处：

1. 顶部 `#SBATCH` 中的 CPU、内存和最长时间；
2. `用户配置区`中的 config、profile、point 或少量命令参数。

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
| FSS 数据、图和拟合 | `fss.sbatch` |
| Hamiltonian 参数优化 | `optimize.sbatch` |
| scaling dimension 图 | `scaling.sbatch` |
| generator ED 与 Lambda 拟合 | `generator.sbatch` |
| tower overlap | `tower.sbatch` |
| orbital entanglement spectrum | `oes.sbatch` |
| real-space entanglement spectrum | `rses.sbatch` |

这些文件中的内存和时间只是适合开始试跑的默认申请量，不代表任意 `nm1` 都足够。
增大 `nm1`、`k` 或保存更多本征向量时，应修改对应文件自己的资源行，不会影响
其他功能。

`plan`、`generator-register`、`fss-plot` 和 `fss-fit` 不做昂贵 ED，通常直接在
登录节点运行，不需要为它们申请专用作业。

每个作业会生成 `slurm-功能名-jobid.out`。程序自己的详细进度和结果仍写入
`output/<功能>/...`；Slurm 文件主要保留启动信息、warning 和错误。
