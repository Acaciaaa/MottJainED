# MottJainED

这是 `mott_jain` 的独立、配置驱动重构版。旧目录和上游 `FuzzifiED.jl`
都不会被本项目修改；项目把 FuzzifiED 固定到已经核对过的 Git commit，
新机器可通过 Julia 的包管理器自动取得同一份源码和对应平台的二进制依赖。

如果你不熟悉 Julia 项目、TOML、命令行或本项目的文件结构，请先读：

**[`docs/START_HERE.md`](docs/START_HERE.md) —— 从零解释每个文件以及所有功能的完整调用流程。**

第一次使用：

```bash
cd MottJainED
julia --project=. scripts/setup.jl
julia --project=. bin/mottjain.jl plan
```

复制并编辑 `config/default.toml` 后运行，例如：

```bash
julia --threads=auto --project=. bin/mottjain.jl critical \
  --config=config/my_run.toml --override=config/critical_profiles/critical5.toml
julia --threads=auto --project=. bin/mottjain.jl fss-all \
  --config=config/my_run.toml --override=config/fss_profiles/fss7.toml
```

optimization 可叠加小模板并临时指定尺寸，例如：

```bash
julia --threads=auto --project=. bin/mottjain.jl optimize \
  --config=config/my_run.toml \
  --override=config/optimization_profiles/mu_uf0_v0.toml \
  --nm1=7 --k=12
```

共形生成元采用“全局候选点 → ED 快照 → 可反复 tower 后处理”的流程：

```bash
# 把看中的 optimization best.csv 加入 config/generator_points.csv
julia --project=. bin/mottjain.jl generator-register \
  --config=config/my_run.toml --point=nm6_candidate_01 \
  --from=output/optimize/mu_uf0_v0_nm6_k70_01/best.csv

# 做/复用一次昂贵 ED，并固定保存 S→dS 拟合出的 Lambda
julia --threads=auto --project=. bin/mottjain.jl generator \
  --config=config/my_run.toml --point=nm6_candidate_01

# 修改 config/generator/nm6_candidate_01/tower.toml 后反复运行，不重算 ED/Lambda
julia --threads=auto --project=. bin/mottjain.jl tower \
  --config=config/my_run.toml --point=nm6_candidate_01
```

查命令和配置参数可再看 [`docs/USER_GUIDE.md`](docs/USER_GUIDE.md)。
