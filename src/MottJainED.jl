module MottJainED

# 这个文件是整个项目的“总目录”。它不做具体物理计算，只负责：
# 1. 载入第三方依赖；2. 依次载入本项目各功能文件；3. 导出给用户调用的名字。
using CSV
using CairoMakie
using DataFrames
using Dates
using FuzzifiED
using JLD2
using LinearAlgebra
using Logging
using Optim
using Printf
using SHA
using SparseArrays
using SpecialFunctions
using Statistics
using TOML
using WignerSymbols

# 本项目的 Hamiltonian 是实矩阵，统一使用 Float64 可以节省内存和时间。
# SilentStd=true 关闭 FuzzifiED 内部较冗长的标准输出，由本项目统一报告进度。
FuzzifiED.ElementType = Float64
FuzzifiED.SilentStd = true

function __init__()
    # Julia 每次载入 MottJainED 时都会执行 __init__；再次设置可避免其他代码
    # 在预编译之后修改 FuzzifiED 全局选项，导致运行行为不一致。
    FuzzifiED.ElementType = Float64
    FuzzifiED.SilentStd = true
end

# include 的顺序就是源码依赖顺序：先定义数据类型，再建立模型和求谱，
# 最后组合成工作流并接到命令行入口。
include("Types.jl")
include("Model.jl")
include("Spectrum.jl")
include("CFT.jl")
include("Storage.jl")
include("Workflows.jl")
include("Conformal.jl")
include("GeneratorWorkflow.jl")
include("Entanglement.jl")
include("CLI.jl")

# export 之后，用户写 `using MottJainED` 就能直接使用这些公开接口；
# 没有 export 的函数通常是内部辅助函数，仍可用 MottJainED.函数名 调用。
export Couplings, SolverSettings, SectorKey, ModelParameters, ModelCache,
       SpectrumState, CFTScore, ScanJob, FSSSettings,
       build_model, hamiltonian_terms, prepare_spectrum,
       solve_spectrum, level_catalog, cft_score, scan_mu,
       run_spectrum_scan, run_gap_scan, run_density_scan,
       run_critical_search, run_parameter_optimization,
       run_fss_scan, plot_fss, fit_fss,
       plot_scaling_dimensions,
       ConformalState, ConformalStateStore, build_conformal_store,
       default_conformal_specs, generator_candidates, run_generator_analysis,
       fit_generator, generator_overlap, project_angular_momentum,
       GeneratorPoint, StoredSectorSpectrum, GeneratorEDSnapshot,
       ensure_generator_registry, load_generator_point,
       register_optimization_point, ensure_generator_case_config,
       ensure_generator_snapshot,
       locate_generator_snapshot, load_generator_snapshot,
       snapshot_states, run_generator_fit, load_generator_fit,
       run_tower_analysis,
       run_orbital_entanglement, run_realspace_entanglement,
       load_config, main

end
