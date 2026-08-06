# 只在一台新机器/服务器第一次安装项目时运行本脚本。
# 它不做任何物理计算，只准备 Julia 软件环境。
import Pkg

# setup.jl 位于 MottJainED/scripts，所以返回上一层得到项目根目录。
root = normpath(joinpath(@__DIR__, ".."))

# 选择本项目的 Project.toml/Manifest.toml。
Pkg.activate(root)

# 告诉 Julia：FuzzifiED 不是从网络安装，而是使用相邻的 ../FuzzifiED.jl。
Pkg.develop(path=normpath(joinpath(root, "..", "FuzzifiED.jl")))

# 根据 Project.toml 和 Manifest.toml 安装/核对全部依赖。
Pkg.instantiate()

# 预编译依赖，减少第一次正式计算的启动时间。
Pkg.precompile()
println("MottJainED environment is ready: $root")
