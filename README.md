# MottJainED

这是 `mott_jain` 的独立、配置驱动重构版。旧目录和本地 `FuzzifiED.jl`
不会被本项目修改；本项目通过相对路径依赖后者。

如果你不熟悉 Julia 项目、TOML、命令行或本项目的文件结构，请先读：

**[`docs/START_HERE.md`](docs/START_HERE.md) —— 从零解释每个文件以及所有功能的完整调用流程。**

第一次使用：

```bash
cd /Users/ruiqi/Documents/hkust/research/fuzzysphere/MottJainED
julia --project=. scripts/setup.jl
julia --project=. bin/mottjain.jl plan
```

复制并编辑 `config/default.toml` 后运行，例如：

```bash
julia --project=. bin/mottjain.jl critical --config=config/my_run.toml
julia --project=. bin/mottjain.jl fss-all --config=config/my_run.toml
```

查命令和配置参数可再看 [`docs/USER_GUIDE.md`](docs/USER_GUIDE.md)。
