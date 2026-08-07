# 本文件把共形生成元计算拆成三个可复用的数据层：
# 1. 全局 CSV 参数点注册表；2. 含完整本征向量的 ED 快照；3. 可反复修改的 tower 后处理。
# 后处理配置改变时只生成新的 analysis 子目录，不重复最昂贵的 ED。

const GENERATOR_REGISTRY_COLUMNS = [
    "point_id", "created_at", "enabled", "nm1",
    "Uf", "Uf0", "U0", "Vf", "Vf0", "V0", "t", "muc",
    "factor", "objective", "q", "cost", "delta_s", "delta_o",
    "score_definition", "source", "notes",
]

"""全局参数表中一行已经过校验、可直接用于 ED 的物理参数点。"""
Base.@kwdef struct GeneratorPoint
    point_id::String
    nm1::Int
    couplings::Couplings
    factor::Union{Nothing,Float64}=nothing
    metadata::Dict{String,Any}=Dict{String,Any}()
end

"""一个 basis 内的一批本征态；basis 只存一次，避免每个态重复序列化。"""
Base.@kwdef struct StoredSectorSpectrum
    family::Symbol                 # :standard 或 :adjoint
    key::SectorKey                 # adjoint family 使用 (0,0) 占位
    basis::Any
    energies::Vector{Float64}
    l2values::Vector{Float64}
    c2values::Vector{Float64}
    ranks::Vector{Int}
    vectors::Vector{Vector{Float64}}
end

"""一次昂贵 ED 的完整、可复用快照。"""
Base.@kwdef struct GeneratorEDSnapshot
    schema_version::Int=2
    snapshot_id::String
    created_at::String
    point::GeneratorPoint
    settings::SolverSettings
    include_adjoint::Bool
    adjoint_k::Int
    sectors::Vector{StoredSectorSpectrum}
    model_summary::Dict{String,Any}
    provenance::Dict{String,Any}
end

_registry_has(row, name::AbstractString) = Symbol(name) in propertynames(row)

function _registry_value(row, name::AbstractString; aliases=String[], required::Bool=true)
    for candidate in (String(name), aliases...)
        _registry_has(row, candidate) || continue
        value = row[Symbol(candidate)]
        if !ismissing(value) && !(value isa AbstractString && isempty(strip(value)))
            return value
        end
    end
    required && throw(ArgumentError("Generator registry row is missing '$name'"))
    return nothing
end

function _registry_float(row, name::AbstractString; aliases=String[], required::Bool=true)
    value = _registry_value(row, name; aliases=aliases, required=required)
    value === nothing && return nothing
    number = value isa Real ? Float64(value) : parse(Float64, strip(String(value)))
    isfinite(number) || throw(ArgumentError("Generator registry '$name' must be finite"))
    return number
end

function _registry_bool(value, default::Bool=true)
    (value === nothing || ismissing(value)) && return default
    value isa Bool && return value
    value isa Integer && return value != 0
    lowercase(strip(String(value))) in ("1", "true", "yes", "on") && return true
    lowercase(strip(String(value))) in ("0", "false", "no", "off") && return false
    throw(ArgumentError("enabled must be true or false, got '$value'"))
end

"""若全局参数表尚不存在，建立只有表头的可手工编辑 CSV。"""
function ensure_generator_registry(path::AbstractString)
    if !isfile(path)
        table = DataFrame((Symbol(name) => Any[] for name in GENERATOR_REGISTRY_COLUMNS)...)
        atomic_csv(path, table)
    end
    return abspath(path)
end

"""为一个 point 建立独立的 generator fit/tower 配置；已有文件永远不覆盖。"""
function ensure_generator_case_config(
    point_id::AbstractString,
    config_root::AbstractString,
    template_root::AbstractString,
)
    occursin(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$", point_id) || throw(ArgumentError(
        "point_id may contain only letters, digits, '.', '_' and '-'",
    ))
    case_directory = joinpath(abspath(config_root), String(point_id))
    mkpath(case_directory)
    for name in ("generator_fit.toml", "tower.toml")
        target = joinpath(case_directory, name)
        isfile(target) && continue
        template = joinpath(abspath(template_root), name)
        isfile(template) || throw(ArgumentError(
            "Generator template not found: $template",
        ))
        cp(template, target)
    end
    return case_directory
end

"""严格按唯一 `point_id` 读取一行；Hamiltonian 参数不会回退到 TOML。"""
function load_generator_point(path::AbstractString, point_id::AbstractString)
    isfile(path) || throw(ArgumentError(
        "Generator registry not found: $path. Create/edit config/generator_points.csv first.",
    ))
    data = CSV.read(path, DataFrame)
    "point_id" in names(data) || throw(ArgumentError("Generator registry has no point_id column"))
    indices = findall(value -> !ismissing(value) && String(value) == String(point_id), data.point_id)
    isempty(indices) && throw(ArgumentError("No generator point '$point_id' in $path"))
    length(indices) == 1 || throw(ArgumentError(
        "Generator point_id '$point_id' appears $(length(indices)) times; point_id must be unique",
    ))
    row = data[only(indices), :]
    enabled = _registry_bool(
        _registry_has(row, "enabled") ? row[:enabled] : nothing,
    )
    enabled || throw(ArgumentError("Generator point '$point_id' is disabled in the registry"))
    nm1_value = _registry_float(row, "nm1")
    nm1 = round(Int, nm1_value)
    nm1_value == nm1 || throw(ArgumentError("nm1 must be an integer, got $nm1_value"))
    nm1 >= 2 || throw(ArgumentError("nm1 must be at least 2"))
    couplings = Couplings(
        Uf=_registry_float(row, "Uf"),
        Uf0=_registry_float(row, "Uf0"),
        U0=_registry_float(row, "U0"),
        Vf=_registry_float(row, "Vf"),
        Vf0=_registry_float(row, "Vf0"),
        V0=_registry_float(row, "V0"),
        t=_registry_float(row, "t"),
        mu=_registry_float(row, "muc"),
    )
    validate(couplings)
    factor = _registry_float(row, "factor"; required=false)
    factor === nothing || factor > 0 || throw(ArgumentError("factor must be positive when provided"))
    metadata = Dict{String,Any}()
    for name in names(data)
        value = row[Symbol(name)]
        metadata[String(name)] = ismissing(value) ? nothing : value
    end
    return GeneratorPoint(
        point_id=String(point_id), nm1=nm1, couplings=couplings,
        factor=factor, metadata=metadata,
    )
end

function _source_row_value(row, name::AbstractString; aliases=String[], required=true)
    return _registry_value(row, name; aliases=aliases, required=required)
end

"""
把一次 optimization 的 `best.csv` 无抄写地加入全局 generator 参数表。

默认拒绝覆盖同名 point；只有用户明确传 `replace=true` 才替换。
"""
function register_optimization_point(
    registry_path::AbstractString,
    point_id::AbstractString,
    best_path::AbstractString;
    notes::AbstractString="",
    replace::Bool=false,
)
    occursin(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$", point_id) || throw(ArgumentError(
        "point_id may contain only letters, digits, '.', '_' and '-'",
    ))
    isfile(best_path) || throw(ArgumentError("Optimization best.csv not found: $best_path"))
    best = CSV.read(best_path, DataFrame)
    nrow(best) == 1 || throw(ArgumentError("Expected exactly one row in $best_path"))
    source = best[1, :]
    absolute_best_path = abspath(best_path)
    relative_best_path = relpath(absolute_best_path, PACKAGE_ROOT)
    recorded_source = startswith(relative_best_path, "..") ?
                      absolute_best_path : relative_best_path
    getnum(name; aliases=String[], required=true) = begin
        value = _source_row_value(source, name; aliases=aliases, required=required)
        value === nothing ? missing : (value isa Real ? Float64(value) : parse(Float64, String(value)))
    end
    row = (
        point_id=String(point_id), created_at=string(now()), enabled=true,
        nm1=round(Int, getnum("nm1")),
        Uf=getnum("Uf"; aliases=["effective_Uf"]),
        Uf0=getnum("Uf0"; aliases=["effective_Uf0"]),
        U0=getnum("U0"; aliases=["effective_U0"]),
        Vf=getnum("Vf"; aliases=["effective_Vf"]),
        Vf0=getnum("Vf0"; aliases=["effective_Vf0"]),
        V0=getnum("V0"; aliases=["effective_V0"]),
        t=getnum("t"; aliases=["effective_t"]),
        muc=getnum("muc"; aliases=["mu", "mu_initial", "effective_mu"]),
        factor=getnum("factor"; required=false),
        objective=getnum("objective"; required=false),
        q=getnum("q"; required=false), cost=getnum("cost"; required=false),
        delta_s=getnum("delta_s"; required=false),
        delta_o=getnum("delta_o"; required=false),
        score_definition=String(_source_row_value(
            source, "score_definition"; required=false,
        ) === nothing ? "" : _source_row_value(source, "score_definition"; required=false)),
        source=recorded_source, notes=String(notes),
    )
    ensure_generator_registry(registry_path)
    existing = CSV.read(registry_path, DataFrame)
    if "point_id" in names(existing)
        same = findall(value -> !ismissing(value) && String(value) == String(point_id), existing.point_id)
        if !isempty(same) && !replace
            throw(ArgumentError("point_id '$point_id' already exists; choose another name or pass replace=true"))
        end
        isempty(same) || deleteat!(existing, same)
    end
    combined = vcat(existing, DataFrame([row]); cols=:union)
    canonical = [name for name in GENERATOR_REGISTRY_COLUMNS if name in names(combined)]
    extras = [name for name in names(combined) if !(name in canonical)]
    select!(combined, vcat(canonical, extras))
    atomic_csv(registry_path, combined)
    return load_generator_point(registry_path, point_id)
end

function _pack_spectrum(states::Vector{SpectrumState}, family::Symbol)
    isempty(states) && return StoredSectorSpectrum[]
    keys = sort!(unique((state.sector.z, state.sector.r) for state in states))
    packed = StoredSectorSpectrum[]
    for (z, r) in keys
        selected = sort(
            filter(state -> state.sector.z == z && state.sector.r == r, states);
            by=state -> state.rank,
        )
        all(state -> state.vector !== nothing && state.basis !== nothing, selected) ||
            error("ED snapshot requires keep_vectors=true")
        push!(packed, StoredSectorSpectrum(
            family=family, key=SectorKey(z, r), basis=first(selected).basis,
            energies=[state.energy for state in selected],
            l2values=[state.l2 for state in selected],
            c2values=[state.c2 for state in selected],
            ranks=[state.rank for state in selected],
            vectors=[state.vector::Vector{Float64} for state in selected],
        ))
    end
    return packed
end

function snapshot_states(snapshot::GeneratorEDSnapshot, family::Symbol=:standard)
    states = SpectrumState[]
    for sector in snapshot.sectors
        sector.family == family || continue
        for i in eachindex(sector.energies)
            push!(states, SpectrumState(
                sector.energies[i], sector.l2values[i], sector.c2values[i],
                sector.key, sector.ranks[i], sector.vectors[i], sector.basis,
            ))
        end
    end
    sort!(states; by=state -> state.energy)
    return states
end

function _jl_tree_signature(root::AbstractString)
    isdir(root) || return "missing"
    files = String[]
    source = joinpath(root, "src")
    if isdir(source)
        for (directory, _, names_in_directory) in walkdir(source)
            for name in names_in_directory
                endswith(name, ".jl") && push!(files, joinpath(directory, name))
            end
        end
    end
    isfile(joinpath(root, "Project.toml")) && push!(files, joinpath(root, "Project.toml"))
    sort!(files)
    io = IOBuffer()
    for path in files
        write(io, relpath(path, root), UInt8(0), read(path), UInt8(0))
    end
    return bytes2hex(sha1(take!(io)))[1:16]
end

function _file_set_signature(paths::Vector{String})
    io = IOBuffer()
    for path in sort(paths)
        isfile(path) || continue
        write(io, basename(path), UInt8(0), read(path), UInt8(0))
    end
    return bytes2hex(sha1(take!(io)))[1:16]
end

function _generator_provenance()
    fuzzified_root = try
        normpath(joinpath(dirname(pathof(FuzzifiED)), ".."))
    catch
        normpath(joinpath(PACKAGE_ROOT, "..", "FuzzifiED.jl"))
    end
    return Dict{String,Any}(
        "package_git_revision" => git_revision(PACKAGE_ROOT),
        "package_source_signature" => _jl_tree_signature(PACKAGE_ROOT),
        # 只有真正决定 basis/Hamiltonian/eigensystem 的文件进入 ED 身份。
        # Project.toml 的本地依赖路径、CSV/overlap/画图代码都不改变已存本征态，
        # 因此不能让这些纯工程改动无意义地要求重算昂贵 ED。
        "ed_source_signature" => _file_set_signature([
            joinpath(PACKAGE_ROOT, "src", "Types.jl"),
            joinpath(PACKAGE_ROOT, "src", "Model.jl"),
            joinpath(PACKAGE_ROOT, "src", "Spectrum.jl"),
        ]),
        "fuzzified_path" => fuzzified_root,
        "fuzzified_git_revision" => git_revision(fuzzified_root),
        "fuzzified_source_signature" => _jl_tree_signature(fuzzified_root),
        "julia_version" => string(VERSION),
        "julia_threads" => Threads.nthreads(),
    )
end

function generator_snapshot_id(
    point::GeneratorPoint,
    settings::SolverSettings;
    include_adjoint::Bool=true,
    adjoint_k::Int=settings.k,
    provenance::AbstractDict=_generator_provenance(),
)
    return stable_id(
        "generator-ed-v1", point.nm1, coupling_vector(point.couplings), settings,
        include_adjoint, adjoint_k,
        get(
            provenance, "ed_source_signature",
            get(provenance, "package_source_signature", "unknown"),
        ),
        get(provenance, "fuzzified_source_signature", "unknown"),
    )
end

function generator_snapshot_directory(
    data_root::AbstractString,
    point::GeneratorPoint,
    snapshot_id::AbstractString,
)
    occursin(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$", point.point_id) || throw(ArgumentError(
        "point_id may contain only letters, digits, '.', '_' and '-'",
    ))
    # snapshot_id 仍用于内部一致性校验，但不再暴露成乱码目录。
    # point_id 必须唯一代表一组 Hamiltonian 参数。
    return joinpath(abspath(data_root), point.point_id)
end

function _model_summary(model::ModelParameters)
    return Dict{String,Any}(
        "name" => String(model.name), "nm1" => model.nm1, "s" => model.s,
        "nf1" => model.nf1, "no1" => model.no1, "nm0" => model.nm0,
        "nf0" => model.nf0, "no0" => model.no0, "no" => model.no,
        "total_charge" => model.no1, "radius_squared" => Float64(model.nm1),
    )
end

"""求原代码 `for_generator_special` 使用的 SU(3) adjoint weight sector。"""
function _solve_generator_adjoint(
    model::ModelParameters,
    couplings::Couplings,
    settings::SolverSettings,
)
    basis = Basis(Confs(model.no, [model.no1, 0, 2, 0], model.qnd))
    basis.dim == 0 && error("Generator adjoint sector is empty")
    h0terms = hamiltonian_terms(
        model, with_coupling(couplings, :mu, 0.0); include_mu=false,
    )
    sector = SectorCache(
        SectorKey(0, 0), basis,
        lower_sparse(float_opmat(Operator(basis, h0terms))),
        lower_sparse(float_opmat(Operator(basis, model.number_f))),
        float_opmat(Operator(basis, model.l2)),
        float_opmat(Operator(basis, model.c2)),
        Float64[],
    )
    energies, vectors = _eigensystem(sector, couplings.mu, settings)
    order = sortperm(energies)
    energies, vectors = energies[order], vectors[:, order]
    _resolve_quantum_numbers!(energies, vectors, sector.l2, sector.c2, settings)
    order = sortperm(energies)
    energies, vectors = energies[order], vectors[:, order]
    states = SpectrumState[]
    for rank in eachindex(energies)
        vector = copy(vectors[:, rank])
        push!(states, SpectrumState(
            energies[rank], real(dot(vector, sector.l2 * vector)),
            real(dot(vector, sector.c2 * vector)), SectorKey(0, 0), rank,
            vector, basis,
        ))
    end
    return states
end

function _generator_with_k(settings::SolverSettings, k::Int)
    return SolverSettings(
        k=k, eig_tol=settings.eig_tol, energy_tol=settings.energy_tol,
        quantum_tol=settings.quantum_tol,
        degeneracy_tol=settings.degeneracy_tol,
        dense_cutoff=settings.dense_cutoff, ncv_extra=settings.ncv_extra,
        warm_start=settings.warm_start,
    )
end

function _snapshot_spectrum_table(
    snapshot::GeneratorEDSnapshot,
    family::Symbol,
    point::GeneratorPoint=snapshot.point,
)
    states = snapshot_states(snapshot, family)
    isempty(states) && return DataFrame()
    ground = minimum(state.energy for state in snapshot_states(snapshot, :standard))
    rows = NamedTuple[]
    for (global_rank, state) in enumerate(states)
        labels = _state_quantum_labels(state, snapshot.settings.quantum_tol)
        gap = state.energy - ground
        push!(rows, (
            family=String(family), global_rank=global_rank,
            energy=state.energy, energy_gap=gap,
            scaled_dimension=point.factor === nothing ? missing : gap / point.factor,
            l2=labels === nothing ? missing : labels[1],
            c2=labels === nothing ? missing : labels[2],
            l2_raw=state.l2, c2_raw=state.c2,
            z=state.sector.z, r=state.sector.r, sector_rank=state.rank,
            basis_dimension=length(state.vector::Vector{Float64}),
        ))
    end
    return DataFrame(rows)
end

function _physical_level_table(snapshot::GeneratorEDSnapshot, family::Symbol)
    states = snapshot_states(snapshot, family)
    isempty(states) && return DataFrame()
    catalog, rejected = level_catalog(
        states; quantum_tol=snapshot.settings.quantum_tol,
        degeneracy_tol=snapshot.settings.degeneracy_tol,
    )
    rows = NamedTuple[]
    for ((l2, c2), levels) in sort!(collect(catalog); by=first)
        for (physical_rank, level) in enumerate(levels)
            for (member_index, member) in enumerate(level.members)
                push!(rows, (
                    family=String(family), l2=l2, c2=c2,
                    physical_rank=physical_rank, member_index=member_index,
                    energy=level.energy, multiplicity=level.multiplicity,
                    z=member.sector.z, r=member.sector.r,
                    sector_rank=member.rank, member_energy=member.energy,
                ))
            end
        end
    end
    for member in rejected
        push!(rows, (
            family=String(family), l2=missing, c2=missing,
            physical_rank=missing, member_index=missing,
            energy=member.energy, multiplicity=0,
            z=member.sector.z, r=member.sector.r,
            sector_rank=member.rank, member_energy=member.energy,
        ))
    end
    return DataFrame(rows)
end

function _point_table(point::GeneratorPoint)
    c = point.couplings
    values = Dict{String,Any}(
        name => (value === nothing ? missing : value)
        for (name, value) in point.metadata
    )
    merge!(values, Dict{String,Any}(
        "point_id" => point.point_id, "nm1" => point.nm1,
        "Uf" => c.Uf, "Uf0" => c.Uf0, "U0" => c.U0,
        "Vf" => c.Vf, "Vf0" => c.Vf0, "V0" => c.V0,
        "t" => c.t, "muc" => c.mu,
        "factor" => point.factor === nothing ? missing : point.factor,
    ))
    order = vcat(
        [name for name in GENERATOR_REGISTRY_COLUMNS if haskey(values, name)],
        sort([name for name in keys(values) if !(name in GENERATOR_REGISTRY_COLUMNS)]),
    )
    table = DataFrame()
    for name in order
        table[!, Symbol(name)] = [values[name]]
    end
    return table
end

function _snapshot_metadata(snapshot::GeneratorEDSnapshot)
    c = snapshot.point.couplings
    sectors = [Dict{String,Any}(
        "family" => String(sector.family), "z" => sector.key.z,
        "r" => sector.key.r, "basis_dimension" => sector.basis.dim,
        "saved_state_count" => length(sector.energies),
    ) for sector in snapshot.sectors]
    registry_metadata = Dict{String,Any}()
    for (name, value) in snapshot.point.metadata
        (value === nothing || ismissing(value)) && continue
        registry_metadata[name] = value isa Union{AbstractString,Real,Bool} ? value : string(value)
    end
    snapshot.point.factor === nothing || (registry_metadata["factor"] = snapshot.point.factor)
    return Dict{String,Any}(
        "schema_version" => snapshot.schema_version,
        "snapshot_id" => snapshot.snapshot_id,
        "created_at" => snapshot.created_at,
        "point_id" => snapshot.point.point_id,
        "model" => snapshot.model_summary,
        "hamiltonian" => Dict(
            "Uf" => c.Uf, "Uf0" => c.Uf0, "U0" => c.U0,
            "Vf" => c.Vf, "Vf0" => c.Vf0, "V0" => c.V0,
            "t" => c.t, "muc" => c.mu,
        ),
        "solver" => Dict(
            "k" => snapshot.settings.k, "adjoint_k" => snapshot.adjoint_k,
            "eig_tol" => snapshot.settings.eig_tol,
            "energy_tol" => snapshot.settings.energy_tol,
            "quantum_tol" => snapshot.settings.quantum_tol,
            "degeneracy_tol" => snapshot.settings.degeneracy_tol,
            "dense_cutoff" => snapshot.settings.dense_cutoff,
            "ncv_extra" => snapshot.settings.ncv_extra,
            "include_adjoint" => snapshot.include_adjoint,
        ),
        "provenance" => snapshot.provenance,
        "registry_row_at_ed_time" => registry_metadata,
        "sectors" => sectors,
    )
end

function materialize_generator_snapshot(
    directory::AbstractString,
    snapshot::GeneratorEDSnapshot;
    point::GeneratorPoint=snapshot.point,
)
    ensure_output(directory)
    # ed_snapshot.jld2 已经保存 point、逐 sector basis 和所有本征向量。
    # 人工挑 tower 态时只需要合并后的 physical-level 表；不再重复导出 point.csv
    # 与逐 sector 的 spectrum_*.csv，避免同一数据出现三四份。
    for family in (:standard, :adjoint)
        states = snapshot_states(snapshot, family)
        isempty(states) && continue
        atomic_csv(
            joinpath(directory, "physical_levels_$(family).csv"),
            _physical_level_table(snapshot, family),
        )
    end
    atomic_toml(joinpath(directory, "snapshot_metadata.toml"), _snapshot_metadata(snapshot))
    return directory
end

function load_generator_snapshot(path::AbstractString)
    isfile(path) || throw(ArgumentError("Generator ED snapshot not found: $path"))
    data = JLD2.load(path)
    haskey(data, "snapshot") || error("Invalid ED snapshot: missing 'snapshot'")
    snapshot = data["snapshot"]
    snapshot isa GeneratorEDSnapshot || error("Unsupported ED snapshot type")
    return snapshot
end

_ed_identity_path(directory::AbstractString) = joinpath(directory, ".ed_identity.toml")

function _same_solver_settings(left::SolverSettings, right::SolverSettings)
    return all(
        field -> getfield(left, field) == getfield(right, field),
        fieldnames(SolverSettings),
    )
end

function _same_ed_request(
    snapshot::GeneratorEDSnapshot,
    point::GeneratorPoint,
    settings::SolverSettings;
    include_adjoint::Bool,
    adjoint_k::Int,
)
    return snapshot.point.point_id == point.point_id &&
           snapshot.point.nm1 == point.nm1 &&
           coupling_vector(snapshot.point.couplings) == coupling_vector(point.couplings) &&
           _same_solver_settings(snapshot.settings, settings) &&
           snapshot.include_adjoint == include_adjoint &&
           snapshot.adjoint_k == adjoint_k
end

"""
验证 ED 身份，并只为旧版过宽身份规则做一次兼容迁移。

旧快照把 `Project.toml` 的本地路径写法算进 source hash；相对路径改成绝对路径也会
被误判成 Hamiltonian 变化。迁移必须同时满足：快照自身旧签名完整、所有物理参数
与 solver 设置逐字段相同、FuzzifiED 源码相同。迁移结果放在隐藏 sidecar；一旦建立，
以后核心源码或设置再变化就不会重复放宽检查。
"""
function _validate_snapshot_identity(
    snapshot::GeneratorEDSnapshot,
    point::GeneratorPoint,
    settings::SolverSettings,
    provenance::AbstractDict,
    identity_signature::AbstractString,
    directory::AbstractString;
    include_adjoint::Bool,
    adjoint_k::Int,
)
    saved_signature = String(get(
        snapshot.provenance, "ed_identity_signature", snapshot.snapshot_id,
    ))
    saved_signature == identity_signature && return true

    identity_path = _ed_identity_path(directory)
    if isfile(identity_path)
        accepted = TOML.parsefile(identity_path)
        return get(accepted, "accepted_identity_signature", "") == identity_signature &&
               get(accepted, "snapshot_saved_identity_signature", "") == saved_signature
    end

    _same_ed_request(
        snapshot, point, settings;
        include_adjoint=include_adjoint, adjoint_k=adjoint_k,
    ) || return false
    snapshot.schema_version == 2 || return false

    saved_fuzzified = String(get(
        snapshot.provenance, "fuzzified_source_signature", "missing",
    ))
    current_fuzzified = String(get(provenance, "fuzzified_source_signature", "missing"))
    saved_fuzzified == current_fuzzified || return false

    # 用快照当时保存的 provenance 重新生成旧签名，确认 JLD2 中的输入与旧身份一致。
    reconstructed_saved = generator_snapshot_id(
        snapshot.point, snapshot.settings;
        include_adjoint=snapshot.include_adjoint,
        adjoint_k=snapshot.adjoint_k,
        provenance=snapshot.provenance,
    )
    reconstructed_saved == saved_signature || return false

    atomic_toml(identity_path, Dict(
        "accepted_identity_signature" => String(identity_signature),
        "snapshot_saved_identity_signature" => saved_signature,
        "migrated_at" => string(now()),
        "reason" => "legacy ED identity included non-physical Project.toml path text",
    ))
    @info "accepted compatible legacy generator ED snapshot" point_id=point.point_id
    return true
end

function _write_ed_identity(
    directory::AbstractString,
    identity_signature::AbstractString,
    saved_signature::AbstractString=identity_signature,
)
    atomic_toml(_ed_identity_path(directory), Dict(
        "accepted_identity_signature" => String(identity_signature),
        "snapshot_saved_identity_signature" => String(saved_signature),
        "created_at" => string(now()),
    ))
    return nothing
end

"""
建立或复用精确匹配当前 point、求解设置和代码版本的 ED 快照。

JLD2 先原子保存；即使后续选态或 overlap 失败，昂贵本征态仍可复用。
"""
function ensure_generator_snapshot(
    point::GeneratorPoint,
    settings::SolverSettings;
    data_root::AbstractString=joinpath(PACKAGE_ROOT, "output", "generator"),
    include_adjoint::Bool=true,
    adjoint_k::Int=settings.k,
    force::Bool=false,
)
    adjoint_k > 0 || throw(ArgumentError("generator.adjoint_k must be positive"))
    provenance = _generator_provenance()
    identity_signature = generator_snapshot_id(
        point, settings; include_adjoint=include_adjoint,
        adjoint_k=adjoint_k, provenance=provenance,
    )
    directory = generator_snapshot_directory(data_root, point, identity_signature)
    path = joinpath(directory, "ed_snapshot.jld2")
    if isfile(path) && !force
        snapshot = load_generator_snapshot(path)
        _validate_snapshot_identity(
            snapshot, point, settings, provenance, identity_signature, directory;
            include_adjoint=include_adjoint, adjoint_k=adjoint_k,
        ) || throw(ArgumentError(
            "Existing ED under point '$(point.point_id)' was made with different Hamiltonian/solver/code settings. " *
            "Use a new point_id, or pass --force only if you intentionally want to replace that point's ED.",
        ))
        materialize_generator_snapshot(directory, snapshot; point=point)
        @info "reusing generator ED snapshot" point_id=point.point_id path
        return (snapshot=snapshot, path=path, directory=directory, reused=true)
    end

    if force && isfile(path)
        # ED 被明确替换后，旧 Lambda/tower 都不再与新本征向量一致。
        for dependent in (
            joinpath(directory, "generator"),
            joinpath(directory, "tower"),
            joinpath(directory, "latest_analysis.toml"),
        )
            if isdir(dependent)
                rm(dependent; recursive=true, force=true)
            elseif isfile(dependent)
                rm(dependent; force=true)
            end
        end
    end

    @info "building generator ED snapshot" point_id=point.point_id nm1=point.nm1
    model = build_model(nm1=point.nm1)
    cache = prepare_spectrum(model, point.couplings, settings)
    standard = solve_spectrum(cache, point.couplings.mu; keep_vectors=true)
    packed = _pack_spectrum(standard, :standard)
    if include_adjoint
        adjoint_settings = _generator_with_k(settings, adjoint_k)
        adjoint = _solve_generator_adjoint(model, point.couplings, adjoint_settings)
        append!(packed, _pack_spectrum(adjoint, :adjoint))
    end
    provenance["ed_identity_signature"] = identity_signature
    snapshot = GeneratorEDSnapshot(
        snapshot_id=point.point_id, created_at=string(now()), point=point,
        settings=settings, include_adjoint=include_adjoint, adjoint_k=adjoint_k,
        sectors=packed, model_summary=_model_summary(model), provenance=provenance,
    )
    ensure_output(directory)
    atomic_jldsave(path; snapshot=snapshot)
    _write_ed_identity(directory, identity_signature)
    materialize_generator_snapshot(directory, snapshot; point=point)
    @info "saved generator ED snapshot" point_id=point.point_id path
    return (snapshot=snapshot, path=path, directory=directory, reused=false)
end

function locate_generator_snapshot(
    point::GeneratorPoint,
    settings::SolverSettings;
    data_root::AbstractString=joinpath(PACKAGE_ROOT, "output", "generator"),
    include_adjoint::Bool=true,
    adjoint_k::Int=settings.k,
)
    provenance = _generator_provenance()
    identity_signature = generator_snapshot_id(
        point, settings; include_adjoint=include_adjoint,
        adjoint_k=adjoint_k, provenance=provenance,
    )
    directory = generator_snapshot_directory(data_root, point, identity_signature)
    path = joinpath(directory, "ed_snapshot.jld2")
    isfile(path) || throw(ArgumentError(
        "No ED snapshot for point '$(point.point_id)'. Run generator --point=$(point.point_id) first.",
    ))
    snapshot = load_generator_snapshot(path)
    _validate_snapshot_identity(
        snapshot, point, settings, provenance, identity_signature, directory;
        include_adjoint=include_adjoint, adjoint_k=adjoint_k,
    ) || throw(ArgumentError(
        "The ED stored under point '$(point.point_id)' does not match the current Hamiltonian/solver/code settings. " *
        "Run generator with a new point_id, or use --force to intentionally replace it.",
    ))
    return (snapshot=snapshot, path=path, directory=directory, reused=true)
end

function _select_snapshot_state(
    snapshot::GeneratorEDSnapshot,
    label::Symbol,
    specification::AbstractDict,
)
    family = Symbol(lowercase(String(get(specification, "family", "standard"))))
    family in (:standard, :adjoint) || throw(ArgumentError(
        "State $label family must be standard or adjoint",
    ))
    states = snapshot_states(snapshot, family)
    isempty(states) && throw(ArgumentError("Snapshot has no $family states for $label"))
    l2 = Int(specification["l2"])
    c2 = Int(specification["c2"])
    rank = Int(get(specification, "rank", 1))
    catalog, _ = level_catalog(
        states; quantum_tol=snapshot.settings.quantum_tol,
        degeneracy_tol=snapshot.settings.degeneracy_tol,
    )
    levels = get(catalog, (l2, c2), PhysicalLevel[])
    length(levels) >= rank || throw(ArgumentError(
        "Missing state $label: family=$family (L2,C2)=($l2,$c2) physical rank=$rank; inspect physical_levels_$(family).csv",
    ))
    members = levels[rank].members
    member, member_index = if haskey(specification, "z") || haskey(specification, "r")
        z = Int(get(specification, "z", first(members).sector.z))
        r = Int(get(specification, "r", first(members).sector.r))
        matches = filter(state -> state.sector.z == z && state.sector.r == r, members)
        isempty(matches) && throw(ArgumentError(
            "State $label rank $rank has no member with Z=$z,R=$r; inspect physical_levels_$(family).csv",
        ))
        chosen = first(matches)
        chosen, findfirst(candidate -> candidate === chosen, members)
    else
        member_index = Int(get(specification, "member", 1))
        1 <= member_index <= length(members) || throw(ArgumentError(
            "State $label member=$member_index is outside 1:$(length(members))",
        ))
        members[member_index], member_index
    end
    return ConformalState(
        label=label, name=String(get(specification, "name", String(label))),
        l2=l2, c2=c2, rank=rank,
        state=member.vector::Vector{Float64}, basis=member.basis,
        energy=member.energy,
    ), family, member, member_index
end

function _ensure_analysis_state!(
    store::ConformalStateStore,
    selections::Dict{Symbol,NamedTuple},
    snapshot::GeneratorEDSnapshot,
    state_specs::AbstractDict,
    label::Symbol,
)
    haskey(store.states, label) && return store[label]
    raw = get(state_specs, String(label), nothing)
    raw === nothing && throw(ArgumentError("tower config has no [states.$label] definition"))
    state, family, member, member_index = _select_snapshot_state(snapshot, label, raw)
    store.states[label] = state
    selections[label] = (
        label=String(label), name=state.name, family=String(family),
        l2=state.l2, c2=state.c2, physical_rank=state.rank,
        member_index=member_index, energy=state.energy,
        z=member.sector.z, r=member.sector.r,
        sector_rank=member.rank,
    )
    return state
end

function _same_angular_overlap(
    input::ConformalState,
    targets::Vector{ConformalState},
    lambda_terms,
    model::ModelParameters;
    target_l::Int,
    ladder::Symbol=:minus,
)
    isempty(targets) && throw(ArgumentError("At least one target is required"))
    basis = first(targets).basis
    all(target -> target.basis === basis, targets) || throw(ArgumentError(
        "All same-angular targets must use the same adjoint basis member",
    ))
    ladder in (:minus, :plus) || throw(ArgumentError("ladder must be minus or plus"))
    lz2 = ladder == :minus ? -2 : 2
    shifted_basis = Basis(Confs(model.no, [model.no1, lz2, 2, 0], model.qnd))
    outward_terms = ladder == :minus ? model.lm : model.lp
    return_terms = ladder == :minus ? model.lp : model.lm
    shifted = real.(Operator(input.basis, shifted_basis, outward_terms) * input.state)
    shifted = real.(Operator(shifted_basis, shifted_basis, lambda_terms) * shifted)
    generated = real.(Operator(shifted_basis, basis, return_terms) * shifted)
    ell = round(Int, (sqrt(1 + 4input.l2) - 1) / 2)
    possible = ell == 0 ? [1] : [ell - 1, ell, ell + 1]
    generated = project_angular_momentum(
        generated, float_opmat(Operator(basis, model.l2)), target_l, possible,
    )
    norm2 = real(dot(generated, generated))
    overlaps = Dict(
        target.label => (norm2 > 0 ? abs2(dot(target.state, generated)) / norm2 : 0.0)
        for target in targets
    )
    return (overlaps=overlaps, total=sum(values(overlaps)), norm2=norm2, generated=generated)
end

function _generator_fit_config_signature(path::AbstractString, snapshot_id)
    isfile(path) || throw(ArgumentError("Generator fit configuration not found: $path"))
    # Lambda 一旦拟合后就是 point 的固定数据。输出/打印代码改变不应迫使用户重拟合；
    # 这里只检查 ED 身份和训练配置。明确 --refit 时再用新实现替换固定 Lambda。
    return stable_id("generator-fit-config-v1", snapshot_id, read(path))
end

"""读取已经固定在 ED 快照旁的 microscopic Lambda；不存在时明确要求先跑 generator。"""
function load_generator_fit(snapshot_directory::AbstractString)
    output = joinpath(snapshot_directory, "generator")
    path = joinpath(output, "generator_fit.jld2")
    isfile(path) || throw(ArgumentError(
        "No fixed generator fit beside this ED snapshot. Run generator --point=POINT_ID first.",
    ))
    data = JLD2.load(path)
    for key in ("fit", "generator_fit_id", "snapshot_id", "fit_signature")
        haskey(data, key) || error("Invalid fixed generator file: missing '$key'")
    end
    return (
        fit=data["fit"], generator_fit_id=String(data["generator_fit_id"]),
        snapshot_id=String(data["snapshot_id"]),
        fit_signature=String(data["fit_signature"]),
        fit_config_signature=get(data, "fit_config_signature", nothing),
        fit_config=get(data, "fit_config", nothing),
        selections=get(data, "selected_states", Dict{Symbol,NamedTuple}()),
        output=output, path=path, reused=true,
    )
end

"""
在一个已经保存的 ED 快照上拟合一次 `Lambda|S>≈|dS>`，并把 Lambda 固定保存。

普通重跑会复用同一拟合；若训练态配置发生变化，必须显式传 `force=true`
（CLI 的 `--refit`），避免 tower 分析过程中无意改变生成元定义。
"""
function run_generator_fit(
    snapshot::GeneratorEDSnapshot,
    point::GeneratorPoint,
    fit_config_path::AbstractString;
    snapshot_directory::AbstractString,
    force::Bool=false,
)
    config = TOML.parsefile(fit_config_path)
    state_specs = get(config, "states", Dict{String,Any}())
    fit_config = get(config, "fit", Dict{String,Any}())
    source_label = Symbol(get(fit_config, "source", "S"))
    target_label = Symbol(get(fit_config, "target", "dS"))
    svd_rtol = Float64(get(fit_config, "svd_rtol", 1.0e-10))
    svd_rtol > 0 || throw(ArgumentError("fit.svd_rtol must be positive"))
    ed_identity_signature = String(get(
        snapshot.provenance, "ed_identity_signature", snapshot.snapshot_id,
    ))
    fit_config_signature = _generator_fit_config_signature(
        fit_config_path, ed_identity_signature,
    )
    generator_fit_id = point.point_id
    output = joinpath(snapshot_directory, "generator")
    fit_path = joinpath(output, "generator_fit.jld2")
    if isfile(fit_path) && !force
        saved = load_generator_fit(snapshot_directory)
        saved.snapshot_id == snapshot.snapshot_id || error("Fixed generator snapshot ID mismatch")
        saved.generator_fit_id == generator_fit_id || error(
            "Fixed generator belongs to another point",
        )
        legacy_config_path = joinpath(output, "resolved_generator_fit.toml")
        same_config = if saved.fit_config_signature !== nothing
            String(saved.fit_config_signature) == fit_config_signature
        elseif saved.fit_config !== nothing
            saved.fit_config == config
        elseif isfile(legacy_config_path)
            TOML.parsefile(legacy_config_path) == config
        else
            false
        end
        same_config || throw(ArgumentError(
            "generator_fit.toml changed after Lambda was fixed. Pass --refit only if you intentionally want to replace it.",
        ))
        @info "reusing fixed generator" point_id=point.point_id generator_fit_id output
        return saved
    end

    ensure_output(output)
    model = build_model(nm1=point.nm1)
    store = ConformalStateStore()
    selections = Dict{Symbol,NamedTuple}()
    source = _ensure_analysis_state!(store, selections, snapshot, state_specs, source_label)
    target = _ensure_analysis_state!(store, selections, snapshot, state_specs, target_label)
    candidates = generator_candidates(model)
    fit = fit_generator(source, target, candidates; svd_rtol=svd_rtol)
    fit_signature = stable_id(
        "fixed-generator-v2", snapshot.snapshot_id, fit_config_signature,
        fit.coefficients, fit.fidelity, fit.numerical_rank,
    )
    atomic_csv(
        joinpath(output, "generator_coefficients.csv"),
        DataFrame(name=fit.names, coefficient=fit.coefficients),
    )
    atomic_csv(joinpath(output, "generator_summary.csv"), DataFrame([(
        point_id=point.point_id, snapshot_id=snapshot.snapshot_id,
        generator_fit_id=generator_fit_id, nm1=point.nm1, muc=point.couplings.mu,
        fit_source=String(source_label), fit_target=String(target_label),
        fidelity=fit.fidelity, numerical_rank=fit.numerical_rank,
        candidate_count=length(candidates), svd_rtol=svd_rtol,
    )]))
    atomic_csv(
        joinpath(output, "generator_selected_states.csv"),
        DataFrame([selections[label] for label in sort!(collect(keys(selections)); by=String)]),
    )
    atomic_jldsave(
        fit_path; fit=fit, selected_states=selections, fit_signature=fit_signature,
        fit_config_signature=fit_config_signature,
        generator_fit_id=generator_fit_id, snapshot_id=snapshot.snapshot_id,
        point=point, fit_config=config,
    )
    @info "saved fixed generator" point_id=point.point_id generator_fit_id output
    return (
        fit=fit, generator_fit_id=generator_fit_id,
        snapshot_id=snapshot.snapshot_id, fit_signature=fit_signature,
        selections=selections,
        output=output, path=fit_path, reused=false,
    )
end

function _tower_config_signature(
    path::AbstractString, point::GeneratorPoint, snapshot_id, generator_fit_id,
)
    isfile(path) || throw(ArgumentError("Tower configuration not found: $path"))
    factor = point.factor === nothing ? "none" : repr(point.factor)
    analysis_source_signature = _file_set_signature([
        joinpath(PACKAGE_ROOT, "src", "Conformal.jl"),
        joinpath(PACKAGE_ROOT, "src", "GeneratorWorkflow.jl"),
    ])
    return stable_id(
        "tower-analysis-v2", snapshot_id, generator_fit_id, factor,
        analysis_source_signature, read(path),
    )
end

function _tower_output_directory(snapshot_directory::AbstractString, signature::AbstractString)
    root = ensure_output(joinpath(snapshot_directory, "tower"))
    maximum_index = 0
    for name in readdir(root)
        startswith(name, "tower_") || continue
        suffix = name[7:end]
        isempty(suffix) && continue
        all(isdigit, suffix) || continue
        index = parse(Int, suffix)
        maximum_index = max(maximum_index, index)
        directory = joinpath(root, name)
        # 新结果把内部签名放进隐藏文件；仍识别旧的可见文件，保证已有 tower 可复用。
        identity_path = joinpath(directory, ".tower_identity.toml")
        legacy_identity_path = joinpath(directory, "tower_identity.toml")
        saved_identity_path = isfile(identity_path) ? identity_path : legacy_identity_path
        if isfile(saved_identity_path)
            saved = TOML.parsefile(saved_identity_path)
            get(saved, "analysis_signature", "") == signature &&
                return (directory=directory, analysis_id=name)
        end
    end
    analysis_id = "tower_" * lpad(maximum_index + 1, 2, '0')
    directory = ensure_output(joinpath(root, analysis_id))
    atomic_toml(joinpath(directory, ".tower_identity.toml"), Dict(
        "analysis_id" => analysis_id,
        "analysis_signature" => String(signature),
        "created_at" => string(now()),
    ))
    return (directory=directory, analysis_id=analysis_id)
end

_tower_value(value) = ismissing(value) ? "--" : @sprintf("%.4f", Float64(value))

function _tower_state_names(state_specs::AbstractDict)
    names = Dict{String,String}()
    for (label, specification) in state_specs
        names[String(label)] = String(get(specification, "name", String(label)))
    end
    return names
end

_tower_display_name(label, names::AbstractDict) = get(names, String(label), String(label))

function _tower_target_cell(row, names::AbstractDict)
    name = _tower_display_name(row.target, names)
    # 原脚本显示的是 (E_target-E_input)/factor。factor 不存在时才退回原始能量差。
    energy = !ismissing(row.scaled_delta) ? row.scaled_delta : row.delta_energy
    energy_text = ismissing(energy) ? "--" : @sprintf("%.2f", Float64(energy))
    return @sprintf("%s(%s)", name, energy_text)
end

function _print_tower_border(io::IO, target_columns::Int)
    print(io, "|--------|----|")
    for _ in 1:target_columns
        print(io, "---------------------|")
    end
    println(io, "--------|")
end

function _print_tower_header(io::IO, target_columns::Int)
    print(io, "| Input  | l' |")
    for index in 1:target_columns
        @printf(io, " %-19s |", "Target $(index) Ovlp")
    end
    println(io, " Total  |")
end

"""
把 tower overlap 打印成旧 `conformal_generator.jl` 的固定格子：
`Input | l' | Target(Δ/f) overlap | ... | Total`。

CSV 仍保留 relation、mode、能量等完整机器可读信息；这里只改变终端展示。
"""
function _print_tower_summary(
    table::DataFrame;
    point_id::AbstractString,
    analysis_id::AbstractString,
    state_names::AbstractDict=Dict{String,String}(),
    io::IO=stdout,
)
    println(io)
    println(io, "Tower overlaps: point=$(point_id), $(analysis_id)")
    if nrow(table) == 0
        println(io, "没有启用的 overlap relation。")
        return nothing
    end

    groups = collect(groupby(table, :relation; sort=false))
    target_columns = max(2, maximum(nrow, groups))
    _print_tower_header(io, target_columns)
    _print_tower_border(io, target_columns)

    skipped = Pair{String,String}[]
    for group in groups
        first_row = first(eachrow(group))
        input = _tower_display_name(first_row.input, state_names)
        target_l = ismissing(first_row.target_l) ? "--" : string(Int(first_row.target_l))
        @printf(io, "| %-6s | %-2s |", input, target_l)
        if String(first_row.status) != "ok"
            @printf(io, " %-12s %6s |", "SKIPPED", "--")
            for _ in 2:target_columns
                print(io, "          -          |")
            end
            @printf(io, " %6s |\n", "--")
            push!(skipped, String(first_row.relation) => String(first_row.message))
            continue
        end
        for (index, row) in enumerate(eachrow(group))
            @printf(io, " %-12s %6s |", _tower_target_cell(row, state_names), _tower_value(row.overlap))
        end
        for _ in (nrow(group) + 1):target_columns
            print(io, "          -          |")
        end
        @printf(io, " %6s |\n", _tower_value(first_row.total_overlap))
    end
    _print_tower_border(io, target_columns)
    for (relation, message) in skipped
        println(io, "Skipped $(relation): $(message)")
    end
    return nothing
end

"""
只读取 ED 快照和已经固定的 Lambda，进行任意多条 tower overlap。

态的 `(L²,C₂,physical rank,member)` 与每条关系由 tower TOML 决定。修改 TOML
会进入新的 `tower_01/tower_02` 连续目录，但不会重新对角化或重新拟合 Lambda。
"""
function run_tower_analysis(
    snapshot::GeneratorEDSnapshot,
    point::GeneratorPoint,
    tower_config_path::AbstractString;
    snapshot_directory::AbstractString,
    generator_fit,
    force::Bool=false,
)
    config = TOML.parsefile(tower_config_path)
    state_specs = get(config, "states", Dict{String,Any}())
    generator_fit.snapshot_id == snapshot.snapshot_id || error(
        "Fixed generator belongs to a different ED snapshot",
    )
    fit = generator_fit.fit
    generator_fit_id = generator_fit.generator_fit_id
    fit_signature = generator_fit.fit_signature
    analysis_signature = _tower_config_signature(
        tower_config_path, point, snapshot.snapshot_id, fit_signature,
    )
    selected_output = _tower_output_directory(snapshot_directory, analysis_signature)
    analysis_id = selected_output.analysis_id
    output = selected_output.directory
    completion = joinpath(output, "analysis_metadata.toml")
    overlap_path = joinpath(output, "tower_overlaps.csv")
    if isfile(completion) && isfile(overlap_path) && !force
        metadata = TOML.parsefile(completion)
        if get(metadata, "status", "") == "complete"
            overlap_table = CSV.read(overlap_path, DataFrame)
            _print_tower_summary(
                overlap_table; point_id=point.point_id, analysis_id=analysis_id,
                state_names=_tower_state_names(state_specs),
            )
            @info "reusing tower analysis" point_id=point.point_id analysis_id output
            return (output=output, analysis_id=analysis_id, reused=true)
        end
    end
    ensure_output(output)
    model = build_model(nm1=point.nm1)
    store = ConformalStateStore()
    selections = Dict{Symbol,NamedTuple}()

    overlap_rows = NamedTuple[]
    saved_vectors = Dict{String,Vector{Float64}}()
    save_vectors = Bool(get(get(config, "output", Dict{String,Any}()), "save_generated_vectors", false))
    relations = get(config, "overlaps", Any[])
    for (relation_index, relation) in enumerate(relations)
        Bool(get(relation, "enabled", true)) || continue
        name = String(get(relation, "name", "relation_$relation_index"))
        required = Bool(get(relation, "required", true))
        input_label = Symbol(relation["input"])
        target_labels = Symbol.(relation["targets"])
        mode = Symbol(lowercase(String(get(relation, "mode", "regular"))))
        target_l = haskey(relation, "target_l") ? Int(relation["target_l"]) : nothing
        try
            input = _ensure_analysis_state!(
                store, selections, snapshot, state_specs, input_label,
            )
            targets = [
                _ensure_analysis_state!(store, selections, snapshot, state_specs, label)
                for label in target_labels
            ]
            result = if mode == :regular
                generator_overlap(
                    input, targets, fit.terms;
                    target_l=target_l, l2_terms=target_l === nothing ? nothing : model.l2,
                )
            elseif mode == :same_angular
                target_l === nothing && throw(ArgumentError(
                    "same_angular relation '$name' requires target_l",
                ))
                _same_angular_overlap(
                    input, targets, fit.terms, model; target_l=target_l,
                    ladder=Symbol(lowercase(String(get(relation, "ladder", "minus")))),
                )
            else
                throw(ArgumentError("Unknown overlap mode '$mode'"))
            end
            save_vectors && (saved_vectors[name] = result.generated)
            for target_state in targets
                delta_energy = target_state.energy - input.energy
                push!(overlap_rows, (
                    status="ok", relation=name, mode=String(mode),
                    input=String(input_label), target=String(target_state.label),
                    target_l=target_l === nothing ? missing : target_l,
                    overlap=result.overlaps[target_state.label],
                    total_overlap=result.total, generated_norm2=result.norm2,
                    input_energy=input.energy, target_energy=target_state.energy,
                    delta_energy=delta_energy,
                    scaled_delta=point.factor === nothing ? missing : delta_energy / point.factor,
                    message="",
                ))
            end
        catch err
            required && rethrow()
            push!(overlap_rows, (
                status="skipped", relation=name, mode=String(mode),
                input=String(input_label), target=join(String.(target_labels), ","),
                target_l=target_l === nothing ? missing : target_l,
                overlap=missing, total_overlap=missing, generated_norm2=missing,
                input_energy=missing, target_energy=missing,
                delta_energy=missing, scaled_delta=missing,
                message=sanitize_error(err),
            ))
        end
    end
    overlap_table = DataFrame(overlap_rows)
    atomic_csv(overlap_path, overlap_table)
    atomic_csv(
        joinpath(output, "selected_states.csv"),
        DataFrame([selections[label] for label in sort!(collect(keys(selections)); by=String)]),
    )
    if save_vectors
        atomic_jldsave(
            joinpath(output, "tower_analysis.jld2");
            selected_states=selections, generated_vectors=saved_vectors,
            snapshot_id=snapshot.snapshot_id, generator_fit_id=generator_fit_id,
            fit_signature=fit_signature, point=point, tower_config=config,
        )
    end
    metadata = Dict{String,Any}(
        "status" => "complete", "created_at" => string(now()),
        "point_id" => point.point_id, "snapshot_id" => snapshot.snapshot_id,
        "analysis_id" => analysis_id, "tower_config" => abspath(tower_config_path),
        "analysis_signature" => analysis_signature,
        "generator_fit_id" => generator_fit_id,
        "fit_signature" => fit_signature,
        "generator_fidelity" => fit.fidelity,
        "generator_numerical_rank" => fit.numerical_rank,
        "factor" => point.factor === nothing ? "not_provided" : point.factor,
        "tower_config_snapshot" => config,
    )
    atomic_toml(completion, metadata)
    _print_tower_summary(
        overlap_table; point_id=point.point_id, analysis_id=analysis_id,
        state_names=_tower_state_names(state_specs),
    )
    @info "saved tower analysis" point_id=point.point_id analysis_id output
    return (
        output=output, analysis_id=analysis_id,
        generator_fit_id=generator_fit_id, store=store, reused=false,
    )
end
