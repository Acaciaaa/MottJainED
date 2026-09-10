# CLI = command-line interface（命令行接口）。这个文件不实现新的物理公式；
# 它只把 `ARGS` 中的命令和 TOML 配置翻译成 Workflows.jl 的 Julia 函数调用。

function _deep_merge_config!(base::Dict{String,Any}, overlay::AbstractDict)
    for (key, value) in pairs(overlay)
        name = String(key)
        if value isa AbstractDict && get(base, name, nothing) isa AbstractDict
            nested = Dict{String,Any}(String(k) => deepcopy(v) for (k, v) in pairs(base[name]))
            base[name] = _deep_merge_config!(nested, value)
        else
            base[name] = deepcopy(value)
        end
    end
    return base
end

"""读取基础 TOML；若给出 override，再递归覆盖其中明确写出的字段。"""
function load_config(
    path::AbstractString=joinpath(PACKAGE_ROOT, "config", "default.toml");
    override=nothing,
)
    isfile(path) || throw(ArgumentError("Configuration file not found: $path"))
    config = TOML.parsefile(path)
    if override !== nothing
        override_path = String(override)
        isfile(override_path) || throw(ArgumentError(
            "Override configuration file not found: $override_path",
        ))
        _deep_merge_config!(config, TOML.parsefile(override_path))
    end
    return config
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

function _output_base(config, override=nothing)
    root = override === nothing ?
           String(_get(_section(config, :output), :root, "output")) : String(override)
    return isabspath(root) ? normpath(root) : normpath(joinpath(PACKAGE_ROOT, root))
end

function _output_root(config, command, override=nothing)
    # 新层级始终“功能在前”：<root>/<command>/<可选 run_name>/。
    # 没有 run_name 时直接写入功能目录，不增加中间默认层。
    feature = joinpath(_output_base(config, override), String(command))
    run_name = strip(String(_get(_section(config, :output), :run_name, "")))
    return isempty(run_name) ? feature : joinpath(feature, run_name)
end

function _safe_case_label(label::AbstractString)
    cleaned = replace(strip(String(label)), r"[^A-Za-z0-9_.-]+" => "_")
    cleaned = strip(cleaned, ['_', '.', '-'])
    isempty(cleaned) && return "run"
    return cleaned
end

function _config_identity(config, feature=:all)
    section_map = Dict(
        :spectrum => ("model", "hamiltonian", "solver", "spectrum"),
        :gap => ("model", "hamiltonian", "solver", "gap"),
        :density => ("model", "hamiltonian", "solver", "density"),
        :critical => ("model", "hamiltonian", "solver", "critical"),
        :optimize => ("model", "solver", "optimization"),
        :fss => ("model", "hamiltonian", "solver", "fss"),
        :scaling => ("model", "hamiltonian", "solver", "scaling"),
        :oes => ("model", "hamiltonian", "solver", "entanglement"),
        :rses => ("model", "hamiltonian", "solver", "entanglement"),
    )
    if feature == :all
        snapshot = deepcopy(config)
        output = get!(snapshot, "output", Dict{String,Any}())
        pop!(output, "root", nothing)
    else
        wanted = get(section_map, Symbol(feature), ())
        snapshot = Dict{String,Any}(
            name => deepcopy(config[name]) for name in wanted if haskey(config, name)
        )
    end
    io = IOBuffer()
    TOML.print(io, snapshot; sorted=true)
    return stable_id("resolved-config-v2", String(feature), String(take!(io)))
end

"""给可变配置分配可读的 `_01/_02/...` 目录；同一配置重跑会复用原编号。"""
function _numbered_case_output(feature_root, label, config; feature=:all)
    ensure_output(feature_root)
    prefix = _safe_case_label(label) * "_"
    identity = _config_identity(config, feature)
    maximum_index = 0
    for name in readdir(feature_root)
        startswith(name, prefix) || continue
        suffix = name[(length(prefix) + 1):end]
        isempty(suffix) && continue
        all(isdigit, suffix) || continue
        index = parse(Int, suffix)
        maximum_index = max(maximum_index, index)
        directory = joinpath(feature_root, name)
        identity_path = joinpath(directory, "case_identity.toml")
        if isfile(identity_path)
            saved = TOML.parsefile(identity_path)
            get(saved, "config_identity", "") == identity && return directory
        end
    end
    directory = joinpath(feature_root, prefix * lpad(maximum_index + 1, 2, '0'))
    ensure_output(directory)
    atomic_toml(joinpath(directory, "case_identity.toml"), Dict(
        "config_identity" => identity,
        "created_at" => string(now()),
        "label" => chop(prefix; tail=1),
    ))
    return directory
end

function _optimization_output(config, override_path, output_override)
    feature_root = joinpath(_output_base(config, output_override), "optimize")
    explicit = strip(String(_get(_section(config, :output), :run_name, "")))
    profile = override_path === nothing ? "optimize" :
              splitext(basename(String(override_path)))[1]
    k = Int(_get(_section(config, :optimization), :k, 70))
    label = isempty(explicit) ? "$(profile)_nm$(_nm1(config))_k$(k)" : explicit
    return _numbered_case_output(feature_root, label, config; feature=:optimize)
end

function _ordinary_task_output(config, feature, output_override)
    feature_root = joinpath(_output_base(config, output_override), String(feature))
    run_name = strip(String(_get(_section(config, :output), :run_name, "")))
    isempty(run_name) && return feature_root
    return _numbered_case_output(
        feature_root, run_name, config; feature=Symbol(feature),
    )
end

"""配置中的相对路径统一相对 MottJainED 根目录，而不是当前 shell 目录。"""
_project_path(path::AbstractString) = isabspath(path) ? normpath(path) : normpath(joinpath(PACKAGE_ROOT, path))

function _generator_request(config, options, settings::SolverSettings; output_override=nothing)
    section = _section(config, :generator)
    point_id = get(options, "point", nothing)
    point_id === nothing && throw(ArgumentError(
        "generator/tower requires --point=POINT_ID from the global generator registry",
    ))
    registry = _project_path(get(
        options, "registry", String(_get(section, :registry, "config/generator_points.csv")),
    ))
    data_root = haskey(options, "data-root") ? _project_path(options["data-root"]) :
                joinpath(_output_base(config, output_override), "generator")
    config_root = _project_path(String(_get(section, :config_root, "config/generator")))
    template_root = joinpath(config_root, "templates")
    case_directory = ensure_generator_case_config(
        point_id, config_root, template_root,
    )
    tower_config = _project_path(get(
        options, "tower-config", joinpath(case_directory, "tower.toml"),
    ))
    fit_config = _project_path(get(
        options, "fit-config", joinpath(case_directory, "generator_fit.toml"),
    ))
    generator_settings = _with_k(settings, Int(_get(section, :k, 80)))
    include_adjoint = Bool(_get(section, :include_adjoint, true))
    adjoint_k = Int(_get(section, :adjoint_k, generator_settings.k))
    point = load_generator_point(registry, point_id)
    return (
        point=point, registry=registry, data_root=data_root,
        case_directory=case_directory,
        tower_config=tower_config, fit_config=fit_config, settings=generator_settings,
        include_adjoint=include_adjoint, adjoint_k=adjoint_k,
    )
end

function _command_output_section(command)
    command == "optimize" && return :optimize
    command in ("fss", "fss-all", "fss-plot", "fss-fit") && return :fss
    command in ("generator-register", "tower") && return :generator
    return Symbol(command)
end

function _apply_cli_config_overrides!(config, command, options, override_path)
    if haskey(options, "nm1")
        model = get!(config, "model", Dict{String,Any}())
        model["nm1"] = parse(Int, options["nm1"])
    end
    if haskey(options, "k")
        section_name = command in ("optimize", "plan") ? "optimization" :
                       command in ("fss", "fss-all", "fss-plot", "fss-fit") ? "fss" :
                       command in ("generator", "tower") ? "generator" :
                       command
        section = get!(config, section_name, Dict{String,Any}())
        section["k"] = parse(Int, options["k"])
    end
    output = get!(config, "output", Dict{String,Any}())
    if haskey(options, "run-name")
        output["run_name"] = options["run-name"]
    end
    return config
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
  julia --project=. bin/mottjain.jl COMMAND [--config=FILE] [--override=FILE]
      [--nm1=N] [--k=N] [--run-name=NAME] [--output=DIR] [--force]

命令：
  plan        只显示任务大小，不进行数值计算
  spectrum    在 hamiltonian.mu 计算筛选后的 rescaled 低能物理能级
  gap         对多个系统大小扫描 scalar gap 和 J gap
  density     扫描 charge-1/charge-3 基态密度
  critical    在给定 mu 网格上选所配置 score 的最小点
  optimize    按 free/values/bounds 优化任意 Hamiltonian 参数组合
  fss         FSS 固定 grid / 沿 size 追踪 μc；--method=grid|optimize|both
  fss-plot    读取已有 FSS CSV 画图，不重新计算
  fss-fit     联合拟合共享的 Delta_inf 和 omega
  fss-all     依次计算 FSS，并分别画 delta_s 和 delta_o
  scaling     画一个参数点的 scaling-dimension spectrum
  generator-register  把选中的 optimization best.csv 加入全局参数点表
  generator   建立/复用 ED 快照，拟合并固定保存 microscopic Lambda
  tower       读取已有 ED 和固定 Lambda，按 tower TOML 选态并计算 overlap
  oes         orbital entanglement spectrum
  rses        real-space entanglement spectrum
  help        显示这段帮助

generator 额外使用 --point=ID；加 --ed-only 只保存 ED，加 --refit 才覆盖已有
Lambda 拟合。tower 用 --tower-config=FILE 反复尝试选态，不重新 ED/拟合 Lambda。

`--override` 只覆盖小配置文件中明确写出的字段；`--nm1`、`--k` 可临时覆盖
常用标量。物理和数值参数都在 config/default.toml（或其副本）中。结果按稳定 job ID
断点保存；只有显式加入 --force 才重新计算已经成功的任务。

""")
end

function _plan(config)
    # 只解析并展示任务规模，不建 basis、不造 Hamiltonian、不做对角化。
    spectrum = _section(config, :spectrum)
    critical = _range(_section(config, :critical))
    fss = _section(config, :fss)
    scan_values = Float64.(_get(fss, :scan_values, collect(1.5:0.5:4.0)))
    fss_mu_count = Int(_get(fss, :mu_count, 9))
    fss_methods = _fss_methods(fss, Dict{String,String}())
    fss_parameter = String(_get(fss, :scan_parameter, "Uf0"))
    fss_mu_min = Float64(_get(fss, :mu_min, 0.0))
    fss_mu_max = Float64(_get(fss, :mu_max, 0.12))
    _, critical_terms, critical_metric = _score_options(
        _section(config, :critical);
        default_definition="critical5", default_metric="q",
    )
    _, fss_terms, fss_metric = _score_options(
        fss; default_definition="fss7", default_metric="cost",
    )
    _, optimization_terms, optimization_metric = _score_options(
        _section(config, :optimization);
        default_definition="optimization8", default_metric="cost",
    )
    optimization_section = _section(config, :optimization)
    optimization_free = Symbol.(_get(
        optimization_section, :free, ["Uf", "Uf0", "Vf0", "V0", "mu"],
    ))
    optimization_algorithm = Symbol(_get(optimization_section, :algorithm, "auto"))
    optimization_algorithm == :auto &&
        (optimization_algorithm = length(optimization_free) == 1 ? :brent : :nelder_mead)
    optimization_values = _optimization_couplings(config, optimization_section)
    if Bool(_get(optimization_section, :tie_u0_to_uf, true))
        ratio = Float64(_get(optimization_section, :u0_over_uf, 9.0))
        optimization_values = with_coupling(
            optimization_values, :U0, ratio * optimization_values.Uf,
        )
    end
    println("配置预览（这里还没有开始数值计算）")
    println("  单尺寸 nm1                 = $(_nm1(config))")
    println("  多尺寸 nm_values           = $(_nm_values(config))")
    println("  spectrum 的 mu             = $(_couplings(config).mu)")
    println("  spectrum 的 L2             = $(Int.(_get(spectrum, :l2_values, [0, 2, 6])))")
    println("  spectrum 的 C2             = $(Int.(_get(spectrum, :c2_values, [0, 3])))")
    println("  spectrum 每个 (L2,C2) 数量 = $(Int(_get(spectrum, :levels_per_block, 7)))")
    println("  spectrum 每 sector 的 k    = $(Int(_get(spectrum, :k, _solver(config).k)))")
    println("  critical 的 mu 点数        = $(length(critical))")
    println("  FSS 外层任务点数          = $(length(_nm_values(config)) * length(scan_values))")
    println("  FSS 外层参数               = $fss_parameter")
    println("  FSS 参数点                 = $scan_values")
    println("  FSS mu 范围                = [$fss_mu_min, $fss_mu_max]")
    println("  FSS 方法                    = $fss_methods")
    :grid in fss_methods && println("  FSS grid 求谱点数         = $(length(_nm_values(config)) * length(scan_values) * fss_mu_count)")
    if :optimize in fss_methods
        println("  FSS optimize 搜索策略     = $(_get(fss, :optimize_strategy, "size_continuation"))")
        anchor_nm = Int(_get(fss, :optimize_anchor_nm, 4))
        anchor_sizes = count(nm -> nm <= anchor_nm, _nm_values(config))
        continuation_sizes = length(_nm_values(config)) - anchor_sizes
        anchor_count = anchor_sizes * length(scan_values) * fss_mu_count
        local_count = continuation_sizes * length(scan_values) *
                      Int(_get(fss, :optimize_local_count, 9))
        println("  FSS anchor size            = N <= $anchor_nm")
        println("  FSS anchor 发现网格 ED    = $anchor_count")
        println("  FSS continuation 初始网格 = $local_count（另加精修、宽 Brent 和必要扩窗）")
    end
    println("  通用每 sector 本征态数 k = $(_solver(config).k)")
    println("  critical 每 sector 的 k       = $(Int(_get(_section(config, :critical), :k, 10)))")
    println("  FSS 每 sector 的 k        = $(Int(_get(fss, :k, 15)))")
    println("  critical score ($(critical_metric)) = $critical_terms")
    println("  FSS score ($(fss_metric))      = $fss_terms")
    println("  optimize score ($(optimization_metric)) = $optimization_terms")
    println("  optimize 每 sector 的 k   = $(Int(_get(optimization_section, :k, 70)))")
    println("  optimize 自由参数          = $optimization_free")
    println("  optimize 算法              = $optimization_algorithm")
    println("  optimize 初值/固定值       = $(coupling_namedtuple(optimization_values))")
    run_name = strip(String(_get(_section(config, :output), :run_name, "")))
    println("  输出标签                    = $(isempty(run_name) ? "（功能目录）" : run_name)")
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

function _optimization_couplings(config, section)
    # [optimization.values] 只服务 optimize：自由参数取这里作为初值，
    # 不在 free 中的参数取这里作为固定值；未写的项才回退到 [hamiltonian]。
    couplings = _couplings(config)
    values = _get(section, :values, Dict{String,Any}())
    for (raw_name, raw_value) in pairs(values)
        name = Symbol(raw_name)
        name in HAMILTONIAN_FIELDS || throw(ArgumentError(
            "Unknown parameter '$raw_name' in [optimization.values]",
        ))
        couplings = with_coupling(couplings, name, Float64(raw_value))
    end
    return validate(couplings)
end

function _score_options(section; default_definition, default_metric=nothing)
    definition = normalize_score_definition(_get(section, :score, default_definition))
    terms = resolve_score_terms(definition, _get(section, :score_terms, nothing))
    metric = normalize_score_metric(
        _get(section, :score_metric, default_metric), definition,
    )
    return definition, terms, metric
end

function _require_case_score(section, command)
    (haskey(section, "score_terms") || haskey(section, "score")) && return nothing
    throw(ArgumentError(
        "$command requires a case profile containing score_terms. " *
        "Pass the corresponding --override=config/..._profiles/CASE.toml file.",
    ))
end

function _fss_methods(section, options)
    raw = haskey(options, "method") ? options["method"] :
          _get(section, :methods, ["grid", "optimize"])
    names = raw isa AbstractVector ? String.(raw) : split(String(raw), ',')
    methods = Symbol.(lowercase.(strip.(names)))
    :both in methods && (methods = [:grid, :optimize])
    methods = unique(methods)
    all(method -> method in (:grid, :optimize), methods) || throw(ArgumentError(
        "FSS method must be grid, optimize, or both",
    ))
    return methods
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
    override_path = get(options, "override", nothing)
    config = load_config(config_path; override=override_path)
    _apply_cli_config_overrides!(config, command, options, override_path)
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
    if command == "critical"
        _require_case_score(_section(config, :critical), command)
    elseif command == "optimize"
        optimization_case = _section(config, :optimization)
        _require_case_score(optimization_case, command)
        for required_section in ("free", "values", "bounds")
            haskey(optimization_case, required_section) || throw(ArgumentError(
                "optimize case profile is missing [optimization].$required_section",
            ))
        end
    elseif command in ("fss", "fss-all")
        _require_case_score(_section(config, :fss), command)
    end
    feature = _command_output_section(command)
    task_output = if command == "optimize"
        _optimization_output(config, override_path, output_override)
    elseif feature == :generator
        point_id = get(options, "point", nothing)
        point_id === nothing && throw(ArgumentError("$command requires --point=POINT_ID"))
        joinpath(_output_base(config, output_override), "generator", point_id)
    else
        _ordinary_task_output(config, feature, output_override)
    end
    log_filename = command == "generator-register" ? "register.log" :
                   feature == :generator ? "$(command).log" : "run.log"
    logging_state = start_task_logging(task_output; filename=log_filename)
    println("详细进度写入：$(logging_state.path)")

    # 从这里开始，每个 elseif 就对应用户手册中的一个可运行功能。
    try
    if command == "spectrum"
        # 单个 mu 的精简物理能级表；不保存本征向量，也不猜测算符身份。
        section = _section(config, :spectrum)
        spectrum_settings = _with_k(settings, Int(_get(section, :k, settings.k)))
        haskey(section, "factor") || throw(ArgumentError(
            "spectrum requires an explicit positive [spectrum].factor; " *
            "it will not infer a scale from tentative operator identities",
        ))
        factor = Float64(_get(section, :factor, NaN))
        write_resolved_config(
            task_output, config;
            base_config=config_path, override_config=override_path,
        )
        run_spectrum(
            nm1, couplings, spectrum_settings;
            l2_values=Int.(_get(section, :l2_values, [0, 2, 6])),
            c2_values=Int.(_get(section, :c2_values, [0, 3])),
            levels_per_block=Int(_get(section, :levels_per_block, 7)),
            factor=factor,
            output=task_output, force=force,
        )
    elseif command == "gap"
        # Workflows.run_gap_scan：多个 nm1 的 scalar/J gap 扫描及两张图。
        section = _section(config, :gap)
        gap_settings = _with_k(settings, Int(_get(section, :k, 5)))
        run_gap_scan(
            _nm_values(config), _range(section), couplings, gap_settings;
            output=task_output, force=force,
        )
    elseif command == "density"
        # density 只需基态，单独使用小 k，不继承全谱的较大 k。
        section = _section(config, :density)
        density_settings = _with_k(settings, Int(_get(section, :k, 3)))
        run_density_scan(
            nm1, _range(section), couplings, density_settings;
            output=task_output, force=force,
        )
    elseif command == "critical"
        # critical 只计算配置中列出的 μ 网格，不调用 Optim。
        section = _section(config, :critical)
        definition, terms, metric = _score_options(
            section; default_definition="critical5", default_metric="q",
        )
        critical_settings = _with_k(settings, Int(_get(section, :k, 10)))
        run_critical_search(
            nm1, _range(section), couplings, critical_settings;
            definition=definition, terms=terms, metric=metric,
            output=task_output,
        )
    elseif command == "optimize"
        # Workflows.run_parameter_optimization：带边界的多耦合参数优化。
        section = _section(config, :optimization)
        free = Symbol.(_get(section, :free, ["Uf", "Uf0", "Vf0", "V0", "mu"]))
        definition, terms, metric = _score_options(
            section; default_definition="optimization8", default_metric="cost",
        )
        optimization_settings = _with_k(settings, Int(_get(section, :k, 70)))
        tie_u0 = Bool(_get(section, :tie_u0_to_uf, true))
        optimization_couplings = _optimization_couplings(config, section)
        optimization_output = task_output
        write_resolved_config(
            optimization_output, config;
            base_config=config_path, override_config=override_path,
        )
        run_parameter_optimization(
            nm1, optimization_couplings, free,
            _optimization_bounds(section, free), optimization_settings;
            max_iterations=Int(_get(section, :max_iterations, 200)),
            algorithm=Symbol(_get(section, :algorithm, "auto")),
            abs_tol=Float64(_get(section, :abs_tol, 1e-4)),
            definition=definition, terms=terms, metric=metric,
            u0_over_uf=tie_u0 ? Float64(_get(section, :u0_over_uf, 9.0)) : nothing,
            penalty=Float64(_get(section, :penalty, 1.0e6)),
            output=optimization_output,
        )
    elseif command in ("fss", "fss-all")
        # fss 只算数据；fss-all 随后分别画 ΔS 和 ΔO，不做尺寸拟合。
        section = _section(config, :fss)
        definition, terms, metric = _score_options(
            section; default_definition="fss7", default_metric="cost",
        )
        methods = _fss_methods(section, options)
        fss_settings = _with_k(settings, Int(_get(section, :k, 15)))
        fss = FSSSettings(
            nm_values=_nm_values(config),
            scan_parameter=Symbol(_get(section, :scan_parameter, "Uf0")),
            scan_values=Float64.(_get(section, :scan_values, collect(1.5:0.5:4.0))),
            mu_min=Float64(_get(section, :mu_min, 0.0)),
            mu_max=Float64(_get(section, :mu_max, 0.12)),
            mu_count=Int(_get(section, :mu_count, 9)),
            methods=methods, score_definition=definition, score_terms=terms,
            score_metric=metric,
            optimize_strategy=Symbol(lowercase(replace(
                String(_get(section, :optimize_strategy, "size_continuation")), '-' => '_',
            ))),
            optimize_anchor_nm=Int(_get(section, :optimize_anchor_nm, 4)),
            optimize_local_half_width=Float64(_get(section, :optimize_local_half_width, 0.02)),
            optimize_local_count=Int(_get(section, :optimize_local_count, 9)),
            optimize_max_expansions=Int(_get(section, :optimize_max_expansions, 3)),
            optimize_abs_tol=Float64(_get(section, :optimize_abs_tol, 1e-4)),
            optimize_max_iterations=Int(_get(section, :optimize_max_iterations, 60)),
            optimize_wide_mode=Symbol(lowercase(replace(
                String(_get(section, :optimize_wide_mode, "always")), '-' => '_',
            ))),
            optimize_wide_adaptive_nm=Int(_get(section, :optimize_wide_adaptive_nm, 6)),
            optimize_wide_audit_first=Bool(_get(section, :optimize_wide_audit_first, false)),
            optimize_wide_audit_all=Bool(_get(section, :optimize_wide_audit_all, false)),
            optimize_wide_jump_tol=Float64(_get(section, :optimize_wide_jump_tol, 0.03)),
            optimize_wide_mu_tol=Float64(_get(section, :optimize_wide_mu_tol, 5e-3)),
            optimize_wide_objective_tol=Float64(_get(
                section, :optimize_wide_objective_tol, 1e-4,
            )),
        )
        directory = task_output
        run_fss_scan(couplings, fss, fss_settings; output=directory, force=force)
        if command == "fss-all"
            for method in methods
                source = joinpath(directory, "fss_$(method)_results.csv")
                data = CSV.read(source, DataFrame)
                for observable in (:delta_s, :delta_o)
                    valid = _valid_fss(data, observable)
                    if nrow(valid) == 0
                        @warn "FSS data were saved, but no valid rows are available; plot skipped" method observable source
                        continue
                    end
                    plot_fss(
                        source; y=observable,
                        output=joinpath(directory, "$(observable)_$(method)_fss.png"),
                    )
                end
            end
        end
    elseif command == "fss-plot"
        # 只读取已有 CSV 画图，不重新对角化。
        section = _section(config, :fss)
        methods = _fss_methods(section, options)
        directory = task_output
        y = Symbol(get(options, "y", "delta_s"))
        if haskey(options, "source")
            plot_fss(options["source"]; y=y)
        else
            for method in methods
                plot_fss(
                    joinpath(directory, "fss_$(method)_results.csv"); y=y,
                    output=joinpath(directory, "$(y)_$(method)_fss.png"),
                )
            end
        end
    elseif command == "fss-fit"
        # 只读取已有 CSV 拟合 Δ∞ 和 ω，不重新对角化。
        section = _section(config, :fss)
        methods = _fss_methods(section, options)
        directory = task_output
        y = Symbol(get(options, "y", "delta_s"))
        if haskey(options, "source")
            fit_fss(options["source"]; y=y)
        else
            for method in methods
                fit_fss(
                    joinpath(directory, "fss_$(method)_results.csv"); y=y,
                    output=directory, label=String(method),
                )
            end
        end
    elseif command == "scaling"
        # Workflows.plot_scaling_dimensions：一个参数点的低能 tower 图。
        section = _section(config, :scaling)
        factor = _get(section, :factor, nothing)
        plot_scaling_dimensions(
            nm1, couplings, settings;
            factor=factor === nothing ? nothing : Float64(factor),
            l2_max=Float64(_get(section, :l2_max, 20.0)),
            c2_max=Float64(_get(section, :c2_max, 8.0)),
            output=task_output,
        )
    elseif command == "generator-register"
        # 把人工确认值得研究的一次 optimization 最优点复制到全局参数表。
        section = _section(config, :generator)
        registry = _project_path(get(
            options, "registry", String(_get(section, :registry, "config/generator_points.csv")),
        ))
        point_id = get(options, "point", nothing)
        point_id === nothing && throw(ArgumentError("generator-register requires --point=NEW_ID"))
        source = get(options, "from", nothing)
        source === nothing && throw(ArgumentError("generator-register requires --from=PATH/TO/best.csv"))
        point = register_optimization_point(
            registry, point_id, _project_path(source);
            notes=get(options, "notes", ""),
            replace=_option_bool(options, "replace", false),
        )
        config_root = _project_path(String(_get(section, :config_root, "config/generator")))
        case_directory = ensure_generator_case_config(
            point.point_id, config_root, joinpath(config_root, "templates"),
        )
        @info "registered generator point" point_id=point.point_id registry case_directory output=task_output
    elseif command in ("generator", "tower")
        # generator 固定 ED 和 Lambda；tower 严格只读二者，便于反复试其它态。
        request = _generator_request(
            config, options, settings; output_override=output_override,
        )
        snapshot_result = if command == "generator"
            ensure_generator_snapshot(
                request.point, request.settings; data_root=request.data_root,
                include_adjoint=request.include_adjoint,
                adjoint_k=request.adjoint_k, force=force,
            )
        else
            locate_generator_snapshot(
                request.point, request.settings; data_root=request.data_root,
                include_adjoint=request.include_adjoint,
                adjoint_k=request.adjoint_k,
            )
        end
        if command == "generator" && _option_bool(options, "ed-only", false)
            @info "ED-only generator task complete" point_id=request.point.point_id path=snapshot_result.path
        elseif command == "generator"
            run_generator_fit(
                snapshot_result.snapshot, request.point, request.fit_config;
                snapshot_directory=snapshot_result.directory,
                force=_option_bool(options, "refit", false) || force,
            )
        else
            fixed_generator = load_generator_fit(snapshot_result.directory)
            run_tower_analysis(
                snapshot_result.snapshot, request.point, request.tower_config;
                snapshot_directory=snapshot_result.directory,
                generator_fit=fixed_generator, force=force,
            )
        end
    elseif command == "oes"
        # Entanglement.run_orbital_entanglement：轨道硬切分纠缠谱。
        section = _section(config, :entanglement)
        run_orbital_entanglement(
            nm1, couplings, settings;
            xi_cut=Float64(_get(section, :xi_cut, 10.0)),
            output=task_output,
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
            output=task_output,
        )
    else
        throw(ArgumentError("Unknown command '$command'. Run the help command."))
    end
    finally
        stop_task_logging(logging_state)
    end
    return 0
end
