module FastED

using CSV
using DataFrames
using Dates
using FuzzifiED
using JLD2
using LinearAlgebra
using MottJainED
using Printf
using SHA
using SparseArrays
using TOML

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const CACHE_SCHEMA_VERSION = 1
const RESULT_SCHEMA_VERSION = 1
const SECTOR_ORDER = [(z=1, r=1), (z=1, r=-1), (z=-1, r=1), (z=-1, r=-1)]

section(config, name) = get(config, String(name), Dict{String,Any}())
getvalue(values, name, default) = get(values, String(name), default)

project_path(path::AbstractString) = isabspath(path) ? normpath(path) :
    normpath(joinpath(PROJECT_ROOT, path))

function safe_label(label::AbstractString)
    cleaned = replace(strip(String(label)), r"[^A-Za-z0-9_.-]+" => "_")
    cleaned = strip(cleaned, ['_', '.', '-'])
    return isempty(cleaned) ? "run" : cleaned
end

function solver_with_k(settings::SolverSettings, k::Int)
    return SolverSettings(
        k=k,
        eig_tol=settings.eig_tol,
        energy_tol=settings.energy_tol,
        quantum_tol=settings.quantum_tol,
        degeneracy_tol=settings.degeneracy_tol,
        dense_cutoff=settings.dense_cutoff,
        ncv_extra=settings.ncv_extra,
        warm_start=settings.warm_start,
    )
end

function solver_values(settings::SolverSettings)
    return [
        settings.k,
        settings.eig_tol,
        settings.energy_tol,
        settings.quantum_tol,
        settings.degeneracy_tol,
        settings.dense_cutoff,
        settings.ncv_extra,
        settings.warm_start,
    ]
end

function solver_dict(settings::SolverSettings)
    return Dict{String,Any}(
        "k" => settings.k,
        "eig_tol" => settings.eig_tol,
        "energy_tol" => settings.energy_tol,
        "quantum_tol" => settings.quantum_tol,
        "degeneracy_tol" => settings.degeneracy_tol,
        "dense_cutoff" => settings.dense_cutoff,
        "ncv_extra" => settings.ncv_extra,
        "warm_start" => settings.warm_start,
    )
end

function coupling_dict(couplings::Couplings; include_mu::Bool=true)
    values = Dict{String,Any}(
        "Uf" => couplings.Uf,
        "Uf0" => couplings.Uf0,
        "U0" => couplings.U0,
        "Vf" => couplings.Vf,
        "Vf0" => couplings.Vf0,
        "V0" => couplings.V0,
        "t" => couplings.t,
    )
    include_mu && (values["mu"] = couplings.mu)
    return values
end

function files_digest(paths; root::Union{Nothing,AbstractString}=nothing)
    buffer = IOBuffer()
    for path in sort!(String.(collect(paths)))
        label = root === nothing ? basename(path) : relpath(path, root)
        write(buffer, label)
        write(buffer, UInt8(0))
        write(buffer, read(path))
        write(buffer, UInt8(0))
    end
    return bytes2hex(sha1(take!(buffer)))
end

function manifest_dependency_identity(name::AbstractString)
    manifest_path = joinpath(
        PROJECT_ROOT, "Manifest-v$(VERSION.major).$(VERSION.minor).toml",
    )
    isfile(manifest_path) || return Dict{String,Any}(
        "version" => "unknown", "git_tree_sha1" => "unknown",
    )
    manifest = TOML.parsefile(manifest_path)
    entries = get(get(manifest, "deps", Dict{String,Any}()), String(name), Any[])
    isempty(entries) && return Dict{String,Any}(
        "version" => "unknown", "git_tree_sha1" => "unknown",
    )
    entry = only(entries)
    return Dict{String,Any}(
        "version" => String(get(entry, "version", "unknown")),
        "git_tree_sha1" => String(get(entry, "git-tree-sha1", "unknown")),
    )
end

function julia_source_files(root::AbstractString)
    files = String[]
    isdir(root) || return files
    for (directory, _, names) in walkdir(root)
        for name in names
            endswith(name, ".jl") && push!(files, joinpath(directory, name))
        end
    end
    return sort!(files)
end

"""Fingerprint only the source that defines the cached matrices and eigenspectrum."""
function source_identity()
    project_files = [
        joinpath(PROJECT_ROOT, "src", "Types.jl"),
        joinpath(PROJECT_ROOT, "src", "Model.jl"),
        joinpath(PROJECT_ROOT, "src", "Spectrum.jl"),
    ]
    fuzz_source = pathof(FuzzifiED)
    fuzz_root = normpath(joinpath(dirname(fuzz_source), ".."))
    fuzz_files = julia_source_files(joinpath(fuzz_root, "src"))
    return Dict{String,Any}(
        "project_git_revision" => MottJainED.git_revision(PROJECT_ROOT),
        "project_source_hash" => files_digest(project_files; root=PROJECT_ROOT),
        "fuzzified_version" => string(Base.pkgversion(FuzzifiED)),
        "fuzzified_path" => fuzz_root,
        "fuzzified_git_revision" => MottJainED.git_revision(fuzz_root),
        "fuzzified_source_hash" => files_digest(fuzz_files; root=fuzz_root),
        "fuzzified_jll" => manifest_dependency_identity("FuzzifiED_jll"),
        "julia_version" => string(VERSION),
    )
end

function load_spec(path::AbstractString)
    config_path = abspath(path)
    isfile(config_path) || throw(ArgumentError("Configuration file not found: $config_path"))
    config = TOML.parsefile(config_path)
    fast = section(config, :fast_ed)
    isempty(fast) && throw(ArgumentError("Profile must define a [fast_ed] section"))

    nm1 = Int(getvalue(section(config, :model), :nm1, 7))
    couplings = Couplings(section(config, :hamiltonian))
    MottJainED.validate(couplings)
    base_solver = MottJainED._solver(config)
    solver = solver_with_k(base_solver, Int(getvalue(fast, :k, base_solver.k)))
    MottJainED.validate(solver)
    mus = Float64.(getvalue(fast, :mus, [couplings.mu]))
    isempty(mus) && throw(ArgumentError("fast_ed.mus cannot be empty"))
    all(isfinite, mus) || throw(ArgumentError("fast_ed.mus must contain only finite values"))
    length(unique(mus)) == length(mus) || throw(ArgumentError("fast_ed.mus must be unique"))

    definition = Symbol(lowercase(replace(
        String(getvalue(fast, :score_definition, "critical5")), '-' => '_',
    )))
    raw_terms = getvalue(fast, :score_terms, nothing)
    terms = raw_terms === nothing ? nothing : Symbol.(String.(raw_terms))
    metric = Symbol(lowercase(String(getvalue(fast, :score_metric, "q"))))
    identity = source_identity()
    fixed = MottJainED.with_coupling(couplings, :mu, 0.0)
    cache_id = MottJainED.stable_id(
        "fast-ed-cache-v$(CACHE_SCHEMA_VERSION)",
        nm1,
        MottJainED.coupling_vector(fixed),
        identity["project_source_hash"],
        identity["fuzzified_version"],
        identity["fuzzified_source_hash"],
        identity["fuzzified_jll"],
    )
    settings_id = MottJainED.stable_id(
        "fast-ed-result-v$(RESULT_SCHEMA_VERSION)", solver_values(solver),
    )
    cache_root = project_path(String(getvalue(fast, :cache_root, "output/fast_ed/cache")))
    result_root = project_path(String(getvalue(fast, :result_root, "output/fast_ed/runs")))
    run_name = safe_label(String(getvalue(fast, :run_name, "n$(nm1)")))
    cache_directory = joinpath(cache_root, "nm$(nm1)_$(cache_id)")
    result_directory = joinpath(
        result_root, run_name, "nm$(nm1)_$(cache_id)", "k$(solver.k)_$(settings_id)",
    )

    return (
        config=config,
        config_path=config_path,
        nm1=nm1,
        couplings=couplings,
        solver=solver,
        mus=mus,
        definition=definition,
        terms=terms,
        metric=metric,
        identity=identity,
        cache_id=cache_id,
        settings_id=settings_id,
        cache_directory=cache_directory,
        result_directory=result_directory,
        run_name=run_name,
    )
end

sector_label(z::Integer, r::Integer) = "z$(z > 0 ? "pos" : "neg")_r$(r > 0 ? "pos" : "neg")"
sector_cache_path(spec, z::Integer, r::Integer) =
    joinpath(spec.cache_directory, "$(sector_label(z, r)).jld2")
cache_manifest_path(spec) = joinpath(spec.cache_directory, "cache_manifest.toml")

function mu_id(mu::Real)
    return MottJainED.stable_id("mu", Float64(mu))
end

mu_directory(spec, mu::Real) = joinpath(spec.result_directory, "mu_$(mu_id(mu))")
sector_result_path(spec, mu::Real, z::Integer, r::Integer) =
    joinpath(mu_directory(spec, mu), "$(sector_label(z, r)).csv")

function validate_cache_data(data, spec, z::Integer, r::Integer)
    Int(data["schema_version"]) == CACHE_SCHEMA_VERSION ||
        throw(ArgumentError("Unsupported sector-cache schema"))
    String(data["cache_id"]) == spec.cache_id ||
        throw(ArgumentError("Sector cache does not match this Hamiltonian/source identity"))
    Int(data["nm1"]) == spec.nm1 || throw(ArgumentError("Sector cache has the wrong nm1"))
    Int(data["z"]) == z || throw(ArgumentError("Sector cache has the wrong Z label"))
    Int(data["r"]) == r || throw(ArgumentError("Sector cache has the wrong R label"))
    return data
end

function cache_header(path::AbstractString)
    names = (
        "schema_version", "cache_id", "nm1", "z", "r", "dimension",
        "h0_nnz", "number_f_nnz", "l2_nnz", "c2_nnz", "build_seconds",
    )
    return JLD2.jldopen(path, "r") do file
        Dict{String,Any}(name => file[name] for name in names)
    end
end

function cache_entry(path::AbstractString, spec, z::Integer, r::Integer)
    data = validate_cache_data(cache_header(path), spec, z, r)
    return Dict{String,Any}(
        "z" => z,
        "r" => r,
        "file" => basename(path),
        "dimension" => Int(data["dimension"]),
        "h0_nnz" => Int(data["h0_nnz"]),
        "number_f_nnz" => Int(data["number_f_nnz"]),
        "l2_nnz" => Int(data["l2_nnz"]),
        "c2_nnz" => Int(data["c2_nnz"]),
        "build_seconds" => Float64(data["build_seconds"]),
        "bytes" => filesize(path),
    )
end

function cache_manifest(spec, entries; complete::Bool)
    return Dict{String,Any}(
        "schema_version" => CACHE_SCHEMA_VERSION,
        "cache_id" => spec.cache_id,
        "complete" => complete,
        "nm1" => spec.nm1,
        "created_at" => string(now()),
        "config_path" => spec.config_path,
        "hamiltonian" => coupling_dict(spec.couplings; include_mu=false),
        "source" => spec.identity,
        "sectors" => entries,
    )
end

function build_sector_cache(model, h0terms, spec, z::Integer, r::Integer, path::AbstractString)
    started = time()
    basis = FuzzifiED.Basis(model.cfs[0], [z, r], model.qnf)
    basis.dim == 0 && return nothing
    h0 = MottJainED.lower_sparse(MottJainED.float_opmat(FuzzifiED.Operator(basis, h0terms)))
    number_f = MottJainED.lower_sparse(MottJainED.float_opmat(
        FuzzifiED.Operator(basis, model.number_f),
    ))
    l2 = MottJainED.lower_sparse(MottJainED.float_opmat(FuzzifiED.Operator(basis, model.l2)))
    c2 = MottJainED.lower_sparse(MottJainED.float_opmat(FuzzifiED.Operator(basis, model.c2)))
    build_seconds = time() - started
    MottJainED.atomic_jldsave(
        path;
        schema_version=CACHE_SCHEMA_VERSION,
        cache_id=spec.cache_id,
        nm1=spec.nm1,
        z=Int(z),
        r=Int(r),
        dimension=Int(basis.dim),
        h0_nnz=nnz(h0),
        number_f_nnz=nnz(number_f),
        l2_nnz=nnz(l2),
        c2_nnz=nnz(c2),
        build_seconds=build_seconds,
        h0=h0,
        number_f=number_f,
        l2=l2,
        c2=c2,
    )
    return cache_entry(path, spec, z, r)
end

"""
Build one persistent sparse-matrix file per non-empty `(Z,R)` sector.

The sectors are built and released sequentially so N=7 never keeps four copies of
H0/Nf/L2/C2 in one Julia process. Existing matching files are reused unless `force=true`.
"""
function prepare_caches(spec; force::Bool=false)
    mkpath(spec.cache_directory)
    if !force && isfile(cache_manifest_path(spec))
        try
            manifest = load_cache_manifest(spec)
            for entry in manifest["sectors"]
                path = joinpath(spec.cache_directory, String(entry["file"]))
                isfile(path) || throw(ArgumentError("Missing sector cache: $path"))
                cache_entry(path, spec, Int(entry["z"]), Int(entry["r"]))
            end
            @info "reusing complete matrix cache" path=spec.cache_directory
            return manifest
        catch err
            @warn "complete cache could not be reused; checking sector files individually" exception=(err, catch_backtrace())
        end
    end
    model = build_model(nm1=spec.nm1)
    fixed = MottJainED.with_coupling(spec.couplings, :mu, 0.0)
    h0terms = hamiltonian_terms(model, fixed; include_mu=false)
    entries = Dict{String,Any}[]
    for key in SECTOR_ORDER
        path = sector_cache_path(spec, key.z, key.r)
        entry = if isfile(path) && !force
            cache_entry(path, spec, key.z, key.r)
        else
            @info "building persistent sector cache" nm1=spec.nm1 z=key.z r=key.r path
            build_sector_cache(model, h0terms, spec, key.z, key.r, path)
        end
        entry === nothing || push!(entries, entry)
        MottJainED.atomic_toml(
            cache_manifest_path(spec), cache_manifest(spec, entries; complete=false),
        )
        GC.gc()
    end
    isempty(entries) && error("No non-empty symmetry sectors were found")
    MottJainED.atomic_toml(
        cache_manifest_path(spec), cache_manifest(spec, entries; complete=true),
    )
    return TOML.parsefile(cache_manifest_path(spec))
end

function load_cache_manifest(spec)
    path = cache_manifest_path(spec)
    isfile(path) || throw(ArgumentError("Cache manifest not found: $path"))
    manifest = TOML.parsefile(path)
    Bool(get(manifest, "complete", false)) || throw(ArgumentError(
        "Cache manifest is incomplete; rerun the prepare command",
    ))
    String(manifest["cache_id"]) == spec.cache_id || throw(ArgumentError(
        "Cache manifest does not match this Hamiltonian/source identity",
    ))
    return manifest
end

function manifest_sector(spec, sector_index::Integer)
    manifest = load_cache_manifest(spec)
    sectors = manifest["sectors"]
    1 <= sector_index <= length(sectors) || throw(ArgumentError(
        "sector-index must be between 1 and $(length(sectors))",
    ))
    return manifest, sectors[sector_index]
end

function result_is_current(path::AbstractString, spec, mu::Real, z::Integer, r::Integer)
    isfile(path) || return false
    try
        data = CSV.read(path, DataFrame)
        nrow(data) > 0 || return false
        dimension = Int(data.sector_dimension[1])
        expected = dimension <= spec.solver.dense_cutoff ? min(spec.solver.k, dimension) :
                   min(spec.solver.k, dimension - 2)
        sort(Int.(data.rank)) == collect(1:expected) || return false
        all(isfinite, data.energy) && all(isfinite, data.l2) && all(isfinite, data.c2) || return false
        return all(String.(data.cache_id) .== spec.cache_id) &&
               all(String.(data.settings_id) .== spec.settings_id) &&
               all(Int.(data.nm1) .== spec.nm1) && all(Int.(data.k) .== spec.solver.k) &&
               all(Int.(data.sector_dimension) .== dimension) &&
               all(isapprox.(Float64.(data.mu), Float64(mu); atol=0.0, rtol=0.0)) &&
               all(Int.(data.z) .== z) && all(Int.(data.r) .== r)
    catch
        return false
    end
end

"""Load one sector once for a sequence of μ values, without keeping eigenvectors."""
function load_sector_solver(spec, sector_index::Integer)
    _, entry = manifest_sector(spec, sector_index)
    z, r = Int(entry["z"]), Int(entry["r"])
    cache_path = joinpath(spec.cache_directory, String(entry["file"]))
    data = validate_cache_data(JLD2.load(cache_path), spec, z, r)
    sector = MottJainED.SectorCache(
        SectorKey(z, r),
        nothing,
        data["h0"],
        data["number_f"],
        MottJainED.hermitian_opmat(data["l2"]),
        MottJainED.hermitian_opmat(data["c2"]),
        Float64[],
    )
    return (sector=sector, dimension=Int(data["dimension"]), cache_path=cache_path,
            cache_id=spec.cache_id, sector_index=Int(sector_index), z=z, r=r)
end

function solve_sector(spec, mu::Real, sector_index::Integer; force::Bool=false, resident=nothing)
    isfinite(mu) || throw(ArgumentError("mu must be finite"))
    _, entry = manifest_sector(spec, sector_index)
    z, r = Int(entry["z"]), Int(entry["r"])
    output = sector_result_path(spec, mu, z, r)
    if !force && result_is_current(output, spec, mu, z, r)
        @info "reusing completed sector result" mu z r output
        return output
    end
    loaded = resident === nothing ? load_sector_solver(spec, sector_index) : resident
    loaded.cache_id == spec.cache_id && loaded.sector_index == sector_index &&
        loaded.z == z && loaded.r == r || throw(ArgumentError("Resident sector identity mismatch"))
    sector = loaded.sector
    cache_path = loaded.cache_path
    started = time()
    energies, vectors = MottJainED._eigensystem(sector, Float64(mu), spec.solver)
    order = sortperm(energies)
    energies = energies[order]
    vectors = vectors[:, order]
    MottJainED._resolve_quantum_numbers!(
        energies, vectors, sector.l2, sector.c2, spec.solver,
    )
    order = sortperm(energies)
    energies = energies[order]
    vectors = vectors[:, order]
    solve_seconds = time() - started
    job_id = MottJainED.stable_id(
        "fast-ed-sector-v$(RESULT_SCHEMA_VERSION)", spec.cache_id,
        spec.settings_id, Float64(mu), z, r,
    )
    rows = [(
        schema_version=RESULT_SCHEMA_VERSION,
        job_id=job_id,
        cache_id=spec.cache_id,
        settings_id=spec.settings_id,
        nm1=spec.nm1,
        mu=Float64(mu),
        z=z,
        r=r,
        rank=rank,
        energy=Float64(energies[rank]),
        l2=Float64(real(dot(vectors[:, rank], sector.l2 * vectors[:, rank]))),
        c2=Float64(real(dot(vectors[:, rank], sector.c2 * vectors[:, rank]))),
        sector_dimension=loaded.dimension,
        k=spec.solver.k,
        eig_tol=spec.solver.eig_tol,
        solve_seconds=solve_seconds,
        julia_threads=Threads.nthreads(),
    ) for rank in eachindex(energies)]
    MottJainED.atomic_csv(output, DataFrame(rows))
    metadata = Dict{String,Any}(
        "schema_version" => RESULT_SCHEMA_VERSION,
        "job_id" => job_id,
        "cache_id" => spec.cache_id,
        "settings_id" => spec.settings_id,
        "nm1" => spec.nm1,
        "mu" => Float64(mu),
        "z" => z,
        "r" => r,
        "states" => length(energies),
        "sector_dimension" => loaded.dimension,
        "solve_seconds" => solve_seconds,
        "julia_threads" => Threads.nthreads(),
        "solver" => solver_dict(spec.solver),
        "cache_file" => cache_path,
        "config_path" => spec.config_path,
        "completed_at" => string(now()),
    )
    MottJainED.atomic_toml(replace(output, r"\.csv$" => ".toml"), metadata)
    return output
end

function solve_sector(spec, mu_index::Integer, sector_index::Integer; force::Bool=false)
    1 <= mu_index <= length(spec.mus) || throw(ArgumentError(
        "mu-index must be between 1 and $(length(spec.mus))",
    ))
    return solve_sector(spec, spec.mus[mu_index], sector_index; force=force)
end

function rows_to_states(data::DataFrame)
    return [SpectrumState(
        Float64(row.energy), Float64(row.l2), Float64(row.c2),
        SectorKey(Int(row.z), Int(row.r)), Int(row.rank), nothing, nothing,
    ) for row in eachrow(data)]
end

function score_summary_row(spec, mu::Real, score::CFTScore, state_count::Integer)
    return (
        cache_id=spec.cache_id,
        settings_id=spec.settings_id,
        nm1=spec.nm1,
        mu=Float64(mu),
        k=spec.solver.k,
        score_valid=score.valid,
        definition=String(score.definition),
        score_terms=join(String.(score.terms), ","),
        metric=String(score.metric),
        objective=score.objective,
        q=score.q,
        cost=score.cost,
        factor=score.factor,
        delta_s=score.delta_s,
        delta_o=score.delta_o,
        state_count=Int(state_count),
        reason=score.reason,
    )
end

function relation_rows(spec, mu::Real, score::CFTScore)
    return [(
        cache_id=spec.cache_id,
        settings_id=spec.settings_id,
        nm1=spec.nm1,
        mu=Float64(mu),
        term=String(score.terms[index]),
        label=score.labels[index],
        raw_gap=score.raw_gaps[index],
        target_gap=score.target_gaps[index],
        scaled_gap=score.raw_gaps[index] / score.factor,
        residual=score.raw_gaps[index] / score.factor - score.target_gaps[index],
    ) for index in eachindex(score.labels)]
end

function collect_mu(spec, mu::Real)
    manifest = load_cache_manifest(spec)
    parts = DataFrame[]
    missing = String[]
    for entry in manifest["sectors"]
        z, r = Int(entry["z"]), Int(entry["r"])
        path = sector_result_path(spec, mu, z, r)
        if !result_is_current(path, spec, mu, z, r)
            push!(missing, "$(sector_label(z, r)): $path")
            continue
        end
        push!(parts, CSV.read(path, DataFrame))
    end
    isempty(missing) || throw(ArgumentError(
        "Missing or incompatible sector results:\n$(join(missing, "\n"))",
    ))
    merged = reduce(vcat, parts)
    sort!(merged, [:energy, :z, :r, :rank])
    states = rows_to_states(merged)
    score = cft_score(
        states;
        settings=spec.solver,
        definition=spec.definition,
        terms=spec.terms,
        metric=spec.metric,
    )
    directory = mu_directory(spec, mu)
    MottJainED.atomic_csv(joinpath(directory, "merged_spectrum.csv"), merged)
    summary = DataFrame([score_summary_row(spec, mu, score, length(states))])
    MottJainED.atomic_csv(joinpath(directory, "score_summary.csv"), summary)
    relations = score.valid ? DataFrame(relation_rows(spec, mu, score)) : DataFrame(
        cache_id=String[], settings_id=String[], nm1=Int[], mu=Float64[],
        term=String[], label=String[], raw_gap=Float64[], target_gap=Float64[],
        scaled_gap=Float64[], residual=Float64[],
    )
    MottJainED.atomic_csv(joinpath(directory, "score_relations.csv"), relations)
    return (states=states, score=score, summary=summary, merged=merged)
end

function collect_mu(spec, mu_index::Integer)
    1 <= mu_index <= length(spec.mus) || throw(ArgumentError(
        "mu-index must be between 1 and $(length(spec.mus))",
    ))
    return collect_mu(spec, spec.mus[mu_index])
end

function collect_all(spec; allow_incomplete::Bool=false)
    rows = NamedTuple[]
    failures = String[]
    for (index, mu) in enumerate(spec.mus)
        result = try
            collect_mu(spec, mu)
        catch err
            allow_incomplete || rethrow()
            push!(failures, "mu-index=$index mu=$mu: $(sprint(showerror, err))")
            continue
        end
        push!(rows, score_summary_row(spec, mu, result.score, length(result.states)))
    end
    if !isempty(rows)
        summary = DataFrame(rows)
        MottJainED.atomic_csv(joinpath(spec.result_directory, "scan_summary.csv"), summary)
        valid = findall(row -> Bool(row.score_valid) && isfinite(row.objective), eachrow(summary))
        if !isempty(valid)
            best_index = valid[argmin(summary.objective[valid])]
            best = summary[best_index:best_index, :]
            best_mu = Float64(best.mu[1])
            MottJainED.atomic_csv(joinpath(spec.result_directory, "best_summary.csv"), best)
            best_directory = mu_directory(spec, best_mu)
            for (source_name, target_name) in (
                ("score_relations.csv", "best_relations.csv"),
                ("merged_spectrum.csv", "best_spectrum.csv"),
            )
                source = joinpath(best_directory, source_name)
                isfile(source) && MottJainED.atomic_csv(
                    joinpath(spec.result_directory, target_name), CSV.read(source, DataFrame),
                )
            end
            MottJainED.atomic_toml(
                joinpath(spec.result_directory, "collection_manifest.toml"),
                Dict{String,Any}(
                    "complete" => isempty(failures) && length(rows) == length(spec.mus),
                    "collected_mu_count" => length(rows),
                    "requested_mu_count" => length(spec.mus),
                    "best_found" => true,
                    "best_mu" => best_mu,
                    "best_mu_id" => mu_id(best_mu),
                    "cache_id" => spec.cache_id,
                    "settings_id" => spec.settings_id,
                    "config_path" => spec.config_path,
                    "collected_at" => string(now()),
                ),
            )
        end
    end
    return (rows=rows, failures=failures)
end

"""Validate that every requested μ was collected before a cache may be released."""
function validate_collection(spec; require_interior::Bool=false)
    manifest_path = joinpath(spec.result_directory, "collection_manifest.toml")
    isfile(manifest_path) || throw(ArgumentError(
        "Collection manifest not found: $manifest_path",
    ))
    manifest = TOML.parsefile(manifest_path)
    Bool(get(manifest, "complete", false)) || throw(ArgumentError(
        "Collection is incomplete; keep the cache and finish the missing sector results",
    ))
    Bool(get(manifest, "best_found", false)) || throw(ArgumentError(
        "Collection has no valid best μ; keep the cache for diagnosis",
    ))
    String(get(manifest, "cache_id", "")) == spec.cache_id || throw(ArgumentError(
        "Collection cache identity does not match this configuration",
    ))
    String(get(manifest, "settings_id", "")) == spec.settings_id || throw(ArgumentError(
        "Collection solver identity does not match this configuration",
    ))
    Int(get(manifest, "requested_mu_count", -1)) == length(spec.mus) ||
        throw(ArgumentError("Collection requested-μ count does not match the configuration"))
    Int(get(manifest, "collected_mu_count", -1)) == length(spec.mus) ||
        throw(ArgumentError("Not every configured μ was collected"))

    required = (
        "scan_summary.csv", "best_summary.csv", "best_relations.csv", "best_spectrum.csv",
    )
    for name in required
        path = joinpath(spec.result_directory, name)
        isfile(path) || throw(ArgumentError("Required collected result is missing: $path"))
    end

    summary = CSV.read(joinpath(spec.result_directory, "scan_summary.csv"), DataFrame)
    nrow(summary) == length(spec.mus) || throw(ArgumentError(
        "scan_summary.csv has $(nrow(summary)) rows for $(length(spec.mus)) configured μ values",
    ))
    all(Bool.(summary.score_valid)) || throw(ArgumentError(
        "At least one μ has an invalid CFT score; keep the cache for diagnosis",
    ))
    all(isfinite, Float64.(summary.objective)) || throw(ArgumentError(
        "At least one μ has a non-finite objective; keep the cache for diagnosis",
    ))
    Set(Float64.(summary.mu)) == Set(spec.mus) || throw(ArgumentError(
        "scan_summary.csv μ values do not match the configuration",
    ))

    best = CSV.read(joinpath(spec.result_directory, "best_summary.csv"), DataFrame)
    nrow(best) == 1 || throw(ArgumentError("best_summary.csv must contain exactly one row"))
    best_mu = Float64(best.mu[1])
    best_mu == Float64(manifest["best_mu"]) || throw(ArgumentError(
        "best_summary.csv and collection_manifest.toml disagree on best μ",
    ))
    best_mu in spec.mus || throw(ArgumentError("Best μ is not part of the configured scan"))
    if require_interior && length(spec.mus) >= 3
        ordered = sort(spec.mus)
        best_mu != first(ordered) && best_mu != last(ordered) || throw(ArgumentError(
            "Best μ=$best_mu is on the scan boundary; keep the cache and extend the scan",
        ))
    end
    return (
        manifest=manifest,
        summary=summary,
        best=best,
        best_mu=best_mu,
        result_directory=spec.result_directory,
    )
end

"""
Delete one exact persistent matrix-cache directory only after collected results pass validation.

Profiles must opt in with `fast_ed.allow_cache_release=true`. This deliberately prevents the
retained central cache and older profiles from being removed by an accidental command.
"""
function release_cache(spec; require_interior::Bool=false)
    fast = section(spec.config, :fast_ed)
    Bool(getvalue(fast, :allow_cache_release, false)) || throw(ArgumentError(
        "This profile does not set fast_ed.allow_cache_release=true; cache release refused",
    ))
    validation = validate_collection(spec; require_interior=require_interior)
    load_cache_manifest(spec)

    cache_root = project_path(String(getvalue(fast, :cache_root, "output/fast_ed/cache")))
    cache_directory = normpath(spec.cache_directory)
    expected_name = "nm$(spec.nm1)_$(spec.cache_id)"
    basename(cache_directory) == expected_name || throw(ArgumentError(
        "Cache directory name is not the exact expected identity: $cache_directory",
    ))
    isdir(cache_root) || throw(ArgumentError("Cache root is missing: $cache_root"))
    isdir(cache_directory) || throw(ArgumentError("Cache directory is missing: $cache_directory"))
    islink(cache_directory) && throw(ArgumentError("Refusing to release a symlinked cache"))
    dirname(realpath(cache_directory)) == realpath(cache_root) || throw(ArgumentError(
        "Cache directory is not a direct child of the configured cache root",
    ))

    bytes = sum(filesize(joinpath(root, name)) for (root, _, names) in walkdir(cache_directory)
                for name in names)
    rm(cache_directory; recursive=true)
    return (
        cache_id=spec.cache_id,
        bytes=bytes,
        best_mu=validation.best_mu,
        result_directory=validation.result_directory,
    )
end

"""Compare an assembled cached result with the unchanged in-process solver."""
function compare_direct(spec, mu_index::Integer)
    spec.nm1 <= 6 || throw(ArgumentError(
        "Direct comparison intentionally stops at N=6; do not duplicate an N=7 run",
    ))
    1 <= mu_index <= length(spec.mus) || throw(ArgumentError(
        "mu-index must be between 1 and $(length(spec.mus))",
    ))
    mu = spec.mus[mu_index]
    cached = collect_mu(spec, mu)
    build_started = time()
    model = build_model(nm1=spec.nm1)
    direct_cache = prepare_spectrum(model, spec.couplings, spec.solver)
    direct_build_seconds = time() - build_started
    solve_started = time()
    direct_states = solve_spectrum(direct_cache, mu)
    direct_solve_seconds = time() - solve_started
    direct_score = cft_score(
        direct_states;
        settings=spec.solver,
        definition=spec.definition,
        terms=spec.terms,
        metric=spec.metric,
    )

    state_key(state) = (state.sector.z, state.sector.r, state.rank)
    cached_by_key = Dict(state_key(state) => state for state in cached.states)
    direct_by_key = Dict(state_key(state) => state for state in direct_states)
    keys_match = keys(cached_by_key) == keys(direct_by_key)
    common = intersect(Set(keys(cached_by_key)), Set(keys(direct_by_key)))
    energy_error = isempty(common) ? Inf : maximum(
        abs(cached_by_key[key].energy - direct_by_key[key].energy) for key in common
    )
    l2_error = isempty(common) ? Inf : maximum(
        abs(cached_by_key[key].l2 - direct_by_key[key].l2) for key in common
    )
    c2_error = isempty(common) ? Inf : maximum(
        abs(cached_by_key[key].c2 - direct_by_key[key].c2) for key in common
    )
    score_comparable = cached.score.valid && direct_score.valid
    q_error = score_comparable ? abs(cached.score.q - direct_score.q) : NaN
    delta_s_error = score_comparable ? abs(cached.score.delta_s - direct_score.delta_s) : NaN
    delta_o_error = score_comparable ? abs(cached.score.delta_o - direct_score.delta_o) : NaN
    passed = keys_match &&
             energy_error <= max(1.0e-10, 20spec.solver.eig_tol) &&
             l2_error <= max(1.0e-8, spec.solver.quantum_tol / 10) &&
             c2_error <= max(1.0e-8, spec.solver.quantum_tol / 10) &&
             cached.score.valid == direct_score.valid
    row = (
        nm1=spec.nm1,
        mu=mu,
        k=spec.solver.k,
        passed=passed,
        keys_match=keys_match,
        cached_states=length(cached.states),
        direct_states=length(direct_states),
        max_energy_abs_error=energy_error,
        max_l2_abs_error=l2_error,
        max_c2_abs_error=c2_error,
        cached_score_valid=cached.score.valid,
        direct_score_valid=direct_score.valid,
        q_abs_error=q_error,
        delta_s_abs_error=delta_s_error,
        delta_o_abs_error=delta_o_error,
        direct_build_seconds=direct_build_seconds,
        direct_solve_seconds=direct_solve_seconds,
    )
    MottJainED.atomic_csv(
        joinpath(mu_directory(spec, mu), "direct_comparison.csv"), DataFrame([row]),
    )
    return row
end

function plan(spec)
    return (
        config=spec.config_path,
        nm1=spec.nm1,
        cache_id=spec.cache_id,
        cache_directory=spec.cache_directory,
        result_directory=spec.result_directory,
        mus=spec.mus,
        sectors=SECTOR_ORDER,
        tasks=length(spec.mus) * length(SECTOR_ORDER),
        k=spec.solver.k,
        threads=Threads.nthreads(),
    )
end

export CACHE_SCHEMA_VERSION, RESULT_SCHEMA_VERSION, SECTOR_ORDER,
       load_spec, prepare_caches, load_cache_manifest, solve_sector,
       collect_mu, collect_all, validate_collection, release_cache,
       compare_direct, plan, sector_result_path

end
