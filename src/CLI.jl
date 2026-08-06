# CLI = command-line interface（命令行接口）。这个文件不实现新的物理公式；
# 它只把 `ARGS` 中的命令和 TOML 配置翻译成 Workflows.jl 的 Julia 函数调用。

"""读取一个 TOML 配置文件，得到嵌套的 `Dict{String,Any}`。"""
function load_config(path::AbstractString=joinpath(PACKAGE_ROOT, "config", "default.toml"))
    isfile(path) || throw(ArgumentError("Configuration file not found: $path"))
    return TOML.parsefile(path)
end

# TOML 中 `[solver]` 等 section 读出后是 Dict；这两个 helper 统一处理缺省值。
_section(config, name) = get(config, String(name), Dict{String,Any}())
_get(section, name, default) = get(section, String(name), default)

function _couplings(config)
    # [hamiltonian] -> Couplings
    return Couplings(_section(config, :hamiltonian))
end

function _solver(config)
    # [solver] -> SolverSettings，并把 TOML 的通用 Number 类型转为明确 Julia 类型。
    section = _section(config, :solver)
    return SolverSettings(
        k=Int(_get(section, :k, 20)),
        eig_tol=Float64(_get(section, :eig_tol, 1e-8)),
        energy_tol=Float64(_get(section, :energy_tol, 1e-7)),
        quantum_tol=Float64(_get(section, :quantum_tol, 2e-3)),
        degeneracy_tol=Float64(_get(section, :degeneracy_tol, 2e-6)),
        dense_cutoff=Int(_get(section, :dense_cutoff, 128)),
        ncv_extra=Int(_get(section, :ncv_extra, 12)),
        warm_start=Bool(_get(section, :warm_start, true)),
    )
end

function _with_k(settings::SolverSettings, k::Int)
    # 某些任务（generator）需要更多低能态：仅替换 k，其他容差保持不变。
    return SolverSettings(
        k=k, eig_tol=settings.eig_tol, energy_tol=settings.energy_tol,
        quantum_tol=settings.quantum_tol, degeneracy_tol=settings.degeneracy_tol,
        dense_cutoff=settings.dense_cutoff, ncv_extra=settings.ncv_extra,
        warm_start=settings.warm_start,
    )
end

function _nm1(config)
    return Int(_get(_section(config, :model), :nm1, 5))
end

function _nm_values(config)
    return Int.(_get(_section(config, :model), :nm_values, [4, 5]))
end

function _range(section; prefix="mu")
    # 支持两种写法：直接给 mus=[...]，或给 mu_min/mu_max/mu_count 生成等距网格。
    direct = _get(section, Symbol(prefix * "s"), nothing)
    direct === nothing || return Float64.(direct)
    lower = Float64(_get(section, Symbol(prefix * "_min"), 0.0))
    upper = Float64(_get(section, Symbol(prefix * "_max"), 0.12))
    count = Int(_get(section, Symbol(prefix * "_count"), 13))
    return collect(range(lower, upper; length=count))
end

function _output_root(config, command, override=nothing)
    # 默认结果层级：<root>/<run_name>/<command>/；相对路径以项目根目录为准。
    root = override === nothing ? String(_get(_section(config, :output), :root, "output")) : String(override)
    root = isabspath(root) ? root : joinpath(PACKAGE_ROOT, root)
    run_name = String(_get(_section(config, :output), :run_name, "default"))
    return joinpath(root, run_name, String(command))
end

function _parse_cli(args)
    # 例如 ["spectrum","--config=my.toml","--force"] 被拆成：
    # command="spectrum"，options=Dict("config"=>"my.toml", "force"=>"true")。
    positionals = String[]
    options = Dict{String,String}()
    for argument in args
        if startswith(argument, "--")
            parts = split(argument[3:end], "="; limit=2)
            options[parts[1]] = length(parts) == 2 ? parts[2] : "true"
        else
            push!(positionals, argument)
        end
    end
    command = isempty(positionals) ? "help" : first(positionals)
    return command, options
end

_option_bool(options, key, default=false) =
    haskey(options, key) ? lowercase(options[key]) in ("1", "true", "yes", "on") : default

function _print_help()
    println("""
MottJainED — 可复现的 SU(3) fuzzy-sphere 计算流程

用法：
  julia --project=. bin/mottjain.jl COMMAND [--config=FILE] [--output=DIR] [--force]

命令：
  plan        只显示任务大小，不进行数值计算
  spectrum    在一组 mu 上计算完整低能谱
  gap         对多个系统大小扫描 singlet gap
  density     扫描 charge-1/charge-3 基态密度
  critical    固定其余 Hamiltonian 参数，只优化 mu
  optimize    同时优化选定的 Hamiltonian 参数
  fss         生成支持断点续算的 finite-size-scaling 数据
  fss-plot    读取已有 FSS CSV 画图，不重新计算
  fss-fit     联合拟合共享的 Delta_inf 和 omega
  fss-all     依次计算 FSS、画 delta_s、尝试联合拟合
  scaling     画一个参数点的 scaling-dimension spectrum
  generator   拟合并保存共形生成元候选系数
  oes         orbital entanglement spectrum
  rses        real-space entanglement spectrum
  help        显示这段帮助

物理和数值参数都在 config/default.toml（或其副本）中。结果按稳定 job ID
断点保存；只有显式加入 --force 才重新计算已经成功的任务。
""")
end

function _plan(config)
    # 只解析并展示任务规模，不建 basis、不造 Hamiltonian、不做对角化。
    spectrum = _range(_section(config, :spectrum))
    fss = _section(config, :fss)
    scan_values = Float64.(_get(fss, :scan_values, collect(1.5:0.5:4.0)))
    println("配置预览（这里还没有开始数值计算）")
    println("  单尺寸 nm1                 = $(_nm1(config))")
    println("  多尺寸 nm_values           = $(_nm_values(config))")
    println("  spectrum 的 mu 点数        = $(length(spectrum))")
    println("  FSS 总任务点数              = $(length(_nm_values(config)) * length(scan_values))")
    println("  每个 sector 的本征态数 k    = $(_solver(config).k)")
    println("  Julia 线程数                = $(Threads.nthreads())")
    return 0
end

function _optimization_bounds(section, free)
    # [optimization.bounds] -> Dict(:参数 => (下界,上界))
    raw = _get(section, :bounds, Dict{String,Any}())
    bounds = Dict{Symbol,Tuple{Float64,Float64}}()
    for name in free
        entry = get(raw, String(name), nothing)
        entry === nothing && throw(ArgumentError("Missing optimization bound for $name"))
        bounds[name] = (Float64(entry[1]), Float64(entry[2]))
    end
    return bounds
end

"""
命令行的总路由函数。

`bin/mottjain.jl` 把 Julia 自动生成的 `ARGS` 原样传进来。本函数依次：解析命令、
读配置、设置线程、构造参数对象，再根据 command 调用一个高级 workflow。
成功返回 0；异常会向上传递，让 shell/Slurm 看到非零退出状态。
"""
function main(args=ARGS)
    command, options = _parse_cli(args)
    command in ("help", "-h", "--help") && (_print_help(); return 0)
    config_path = get(options, "config", joinpath(PACKAGE_ROOT, "config", "default.toml"))
    config = load_config(config_path)
    command == "plan" && return _plan(config)
    output_override = get(options, "output", nothing)
    force = _option_bool(options, "force", false)
    couplings = _couplings(config)
    settings = _solver(config)
    solver_section = _section(config, :solver)
    BLAS.set_num_threads(Int(_get(solver_section, :blas_threads, 1)))
    configured_threads = Int(_get(solver_section, :fuzzified_threads, 0))
    FuzzifiED.NumThreads = configured_threads > 0 ? configured_threads : Threads.nthreads()
    nm1 = _nm1(config)

    # 从这里开始，每个 elseif 就对应用户手册中的一个可运行功能。
    if command == "spectrum"
        # Workflows.run_spectrum_scan：逐 μ 求谱与 tower score。
        run_spectrum_scan(
            nm1, _range(_section(config, :spectrum)), couplings, settings;
            output=_output_root(config, :spectrum, output_override), force=force,
            keep_vectors=_option_bool(options, "keep-vectors", false),
        )
    elseif command == "gap"
        # Workflows.run_gap_scan：多个 nm1 的 singlet gap 扫描及交叉图。
        run_gap_scan(
            _nm_values(config), _range(_section(config, :gap)), couplings, settings;
            output=_output_root(config, :gap, output_override), force=force,
        )
    elseif command == "density"
        # Workflows.run_density_scan：基态 Nf/N0 随 μ 变化。
        run_density_scan(
            nm1, _range(_section(config, :density)), couplings, settings;
            output=_output_root(config, :density, output_override), force=force,
        )
    elseif command == "critical"
        # Workflows.run_critical_search：固定其他系数，仅优化 μ。
        section = _section(config, :critical)
        run_critical_search(
            nm1, couplings, settings;
            mu_min=Float64(_get(section, :mu_min, 0.0)),
            mu_max=Float64(_get(section, :mu_max, 0.12)),
            coarse_points=Int(_get(section, :coarse_points, 9)),
            output=_output_root(config, :critical, output_override),
        )
    elseif command == "optimize"
        # Workflows.run_parameter_optimization：带边界的多耦合参数优化。
        section = _section(config, :optimization)
        free = Symbol.(_get(section, :free, ["Uf", "Uf0", "Vf0", "V0", "mu"]))
        run_parameter_optimization(
            nm1, couplings, free, _optimization_bounds(section, free), settings;
            max_iterations=Int(_get(section, :max_iterations, 200)),
            output=_output_root(config, :optimization, output_override),
        )
    elseif command in ("fss", "fss-all")
        # fss 只算数据；fss-all 随后还画 ΔS 并尝试联合外推。
        section = _section(config, :fss)
        fss = FSSSettings(
            nm_values=_nm_values(config),
            scan_parameter=Symbol(_get(section, :scan_parameter, "Uf0")),
            scan_values=Float64.(_get(section, :scan_values, collect(1.5:0.5:4.0))),
            mu_min=Float64(_get(section, :mu_min, 0.0)),
            mu_max=Float64(_get(section, :mu_max, 0.12)),
            coarse_points=Int(_get(section, :coarse_points, 9)),
            mu_abs_tol=Float64(_get(section, :mu_abs_tol, 1e-5)),
            max_iterations=Int(_get(section, :max_iterations, 60)),
        )
        directory = _output_root(config, :fss, output_override)
        run_fss_scan(couplings, fss, settings; output=directory, force=force)
        if command == "fss-all"
            source = joinpath(directory, "fss_results.csv")
            plot_fss(source; y=:delta_s)
            try
                fit_fss(source; y=:delta_s)
            catch err
                @warn "FSS data were saved, but the joint fit is not yet identifiable" error=sanitize_error(err)
            end
        end
    elseif command == "fss-plot"
        # 只读取已有 CSV 画图，不重新对角化。
        source = get(options, "source", joinpath(_output_root(config, :fss, output_override), "fss_results.csv"))
        y = Symbol(get(options, "y", "delta_s"))
        plot_fss(source; y=y)
    elseif command == "fss-fit"
        # 只读取已有 CSV 拟合 Δ∞ 和 ω，不重新对角化。
        source = get(options, "source", joinpath(_output_root(config, :fss, output_override), "fss_results.csv"))
        y = Symbol(get(options, "y", "delta_s"))
        fit_fss(source; y=y)
    elseif command == "scaling"
        # Workflows.plot_scaling_dimensions：一个参数点的低能 tower 图。
        section = _section(config, :scaling)
        factor = _get(section, :factor, nothing)
        plot_scaling_dimensions(
            nm1, couplings, settings;
            factor=factor === nothing ? nothing : Float64(factor),
            l2_max=Float64(_get(section, :l2_max, 20.0)),
            c2_max=Float64(_get(section, :c2_max, 12.0)),
            output=_output_root(config, :scaling, output_override),
        )
    elseif command == "generator"
        # Conformal.run_generator_analysis：求更多态并拟合 microscopic Λ。
        generator_settings = _with_k(
            settings, Int(_get(_section(config, :generator), :k, 80)),
        )
        run_generator_analysis(
            nm1, couplings, generator_settings;
            output=_output_root(config, :generator, output_override),
        )
    elseif command == "oes"
        # Entanglement.run_orbital_entanglement：轨道硬切分纠缠谱。
        section = _section(config, :entanglement)
        run_orbital_entanglement(
            nm1, couplings, settings;
            xi_cut=Float64(_get(section, :xi_cut, 10.0)),
            output=_output_root(config, :oes, output_override),
        )
    elseif command == "rses"
        # Entanglement.run_realspace_entanglement：球冠实空间切分纠缠谱。
        section = _section(config, :entanglement)
        run_realspace_entanglement(
            nm1, couplings, settings;
            x=Float64(_get(section, :realspace_x, 0.5)),
            qa_half_window=Int(_get(section, :qa_half_window, 2)),
            lz2_cap=Int(_get(section, :lz2_cap, 14)),
            xi_cut=Float64(_get(section, :xi_cut, 10.0)),
            output=_output_root(config, :rses, output_override),
        )
    else
        throw(ArgumentError("Unknown command '$command'. Run the help command."))
    end
    return 0
end
