# 只在一台新机器/服务器第一次安装项目时运行本脚本。
# 它不做任何物理计算，只准备 Julia 软件环境。
import Pkg

# setup.jl 位于 MottJainED/scripts，所以返回上一层得到项目根目录。
root = normpath(joinpath(@__DIR__, ".."))

# 选择本项目的 Project.toml/Manifest.toml。
Pkg.activate(root)

# 根据 Project.toml 和 Manifest.toml 安装/核对全部依赖。
# FuzzifiED 已固定到核对过的 GitHub commit；Pkg 会自动下载同一份源码，
# 不要求服务器上存在相邻的 FuzzifiED.jl 文件夹。
Pkg.instantiate()

# 预编译依赖，减少第一次正式计算的启动时间。
Pkg.precompile()
println("MottJainED environment is ready: $root")
