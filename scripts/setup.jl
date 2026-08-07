# 只在一台新机器/服务器第一次安装项目时运行本脚本。
# 它不做任何物理计算，只准备 Julia 软件环境。
import Pkg

# setup.jl 位于 MottJainED/scripts，所以返回上一层得到项目根目录。
root = normpath(joinpath(@__DIR__, ".."))

# 选择本项目；Julia会自动优先读取与自身minor版本匹配的Manifest。
Pkg.activate(root)

# 根据Project和当前版本可用的Manifest安装/核对全部依赖。
# FuzzifiED 已固定到核对过的 GitHub commit；Pkg 会自动下载同一份源码，
# 不要求服务器上存在相邻的 FuzzifiED.jl 文件夹。
Pkg.instantiate()

# Julia 1.11 与 1.12 的 stdlib/JLL 依赖图不同，不能共用同一份通用 Manifest。
# 首次在一个新 Julia minor 版本上解析后，把 Pkg 生成的通用文件改成版本专用名；
# 以后该版本会优先读取它，其他 Julia minor 版本则不会误读。
generic_manifest = joinpath(root, "Manifest.toml")
versioned_manifest = joinpath(root, "Manifest-v$(VERSION.major).$(VERSION.minor).toml")
if !isfile(versioned_manifest) && isfile(generic_manifest)
    mv(generic_manifest, versioned_manifest)
    println("Created Julia $(VERSION.major).$(VERSION.minor) manifest: $versioned_manifest")
end

# 预编译依赖，减少第一次正式计算的启动时间。
Pkg.precompile()
println("MottJainED environment is ready: $root")
