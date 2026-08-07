#!/usr/bin/env julia

# 这个文件只是“命令行入口”，不包含具体物理计算。
# 真正的模型、谱计算、FSS 等代码位于 ../src/。

# Pkg 是 Julia 自带的项目/依赖管理工具。
import Pkg

# @__DIR__ 是本文件所在的 MottJainED/bin；".." 回到 MottJainED 根目录。
# activate 后，Julia 会使用根目录里的 Project.toml 和 Manifest.toml。
Pkg.activate(normpath(joinpath(@__DIR__, "..")); io=devnull)

# 加载 src/MottJainED.jl；该文件继续 include src 中的各功能文件。
using MottJainED

# ARGS 是终端中写在 mottjain.jl 后面的字符串。
# 例如 `spectrum --config=config/my_run.toml` 会成为 ARGS 的两个元素。
# main 解析这些参数、选择功能；返回 0 表示任务成功，exit 将状态交给操作系统/Slurm。
exit(MottJainED.main(ARGS))
