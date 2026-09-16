module EndpointEntanglement

using CSV
using CairoMakie
using DataFrames
using Dates
using FuzzifiED
using JLD2
using LinearAlgebra
using MottJainED
using TOML

export EndpointSpec, load_spec, plan, run_endpoint_entanglement,
       fiqh_root_lz2, laughlin_root_lz2

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))

Base.@kwdef struct EndpointSpec
    config_path::String
    output::String
    nm1::Int
    couplings::Couplings
    solver::SolverSettings
    mu_left::Float64
    mu_right::Float64
    cut_x::Float64
    qa::Int
    f3a::Int
    f8a::Int
    lz2_values::Vector{Int}
    left_expected_counts::Vector{Int}
    right_expected_counts::Vector{Int}
    nm1_a::Int
    nm0_a::Int
    lambda_plot_min::Float64
    plot_xi_max::Float64
    phase_fraction_min::Float64
    raw_config::Dict{String,Any}
end

_section(config, name) = get(config, String(name), Dict{String,Any}())
_get(section, name, default) = get(section, String(name), default)
_project_path(path) = isabspath(path) ? normpath(path) : normpath(joinpath(PROJECT_ROOT, path))

function _solver(config)
    section = _section(config, :solver)
    return SolverSettings(
        k=Int(_get(section, :k, 2)),
        eig_tol=Float64(_get(section, :eig_tol, 1e-9)),
        energy_tol=Float64(_get(section, :energy_tol, 1e-7)),
        quantum_tol=Float64(_get(section, :quantum_tol, 2e-3)),
        degeneracy_tol=Float64(_get(section, :degeneracy_tol, 2e-6)),
        dense_cutoff=Int(_get(section, :dense_cutoff, 128)),
        ncv_extra=Int(_get(section, :ncv_extra, 12)),
        # The two endpoints are separate Hamiltonians.  Reusing one endpoint as
        # the Krylov seed of the other is unnecessary and complicates auditing.
        warm_start=false,
    )
end

"""The compact fIQH root has `n_per_flavor` lowest orbitals filled per flavor."""
function fiqh_root_lz2(nm1::Int, n_per_flavor::Int)
    0 <= n_per_flavor <= nm1 || throw(ArgumentError("invalid fIQH particle number"))
    orbitals = collect(-(nm1 - 1):2:(nm1 - 1))
    return 3sum(orbitals[1:n_per_flavor])
end

"""The fermionic 1/3 Laughlin root occupies orbital indices 1,4,7,... ."""
function laughlin_root_lz2(nm0::Int, particle_count::Int)
    particle_count >= 0 || throw(ArgumentError("invalid Laughlin particle number"))
    orbitals = collect(-(nm0 - 1):2:(nm0 - 1))
    positions = 1 .+ 3collect(0:particle_count-1)
    isempty(positions) && return 0
    maximum(positions) <= nm0 || throw(ArgumentError("Laughlin root does not fit in the orbitals"))
    return sum(orbitals[positions])
end

function _validate(spec::EndpointSpec)
    iseven(spec.nm1) || throw(ArgumentError("endpoint workflow currently requires even nm1"))
    spec.qa == 3spec.nm1 ÷ 2 || throw(ArgumentError(
        "QA must be half of the total charge 3*nm1",
    ))
    spec.qa % 3 == 0 || throw(ArgumentError("QA must be divisible by three"))
    spec.f3a == 0 && spec.f8a == 0 || throw(ArgumentError(
        "the audited endpoint workflow currently fixes F3A=F8A=0",
    ))
    0 < spec.cut_x < 1 || throw(ArgumentError("cut_x must lie between zero and one"))
    0 < spec.nm1_a < spec.nm1 || throw(ArgumentError("invalid charge-1 orbital cut"))
    nm0 = 3spec.nm1 - 2
    0 < spec.nm0_a < nm0 || throw(ArgumentError("invalid charge-3 orbital cut"))
    isfinite(spec.mu_left) && isfinite(spec.mu_right) || throw(ArgumentError(
        "endpoint chemical potentials must be finite",
    ))
    spec.mu_left < spec.mu_right || throw(ArgumentError("mu_left must be below mu_right"))
    length(spec.lz2_values) == length(spec.left_expected_counts) ==
        length(spec.right_expected_counts) || throw(ArgumentError(
        "lz2_values and expected-count lists must have equal lengths",
    ))
    all(>(0), spec.left_expected_counts) || throw(ArgumentError("left counts must be positive"))
    all(>(0), spec.right_expected_counts) || throw(ArgumentError("right counts must be positive"))
    issorted(spec.lz2_values) || throw(ArgumentError("lz2_values must be sorted"))
    all(diff(spec.lz2_values) .== 2) || throw(ArgumentError(
        "successive edge sectors must differ by two in 2LzA",
    ))
    n = spec.qa ÷ 3
    left_root = fiqh_root_lz2(spec.nm1, n)
    right_root = laughlin_root_lz2(nm0, n)
    first(spec.lz2_values) == left_root == right_root || throw(ArgumentError(
        "first Lz2A must be the common compact root; expected $left_root",
    ))
    0 < spec.lambda_plot_min < 1 || throw(ArgumentError("invalid lambda_plot_min"))
    spec.plot_xi_max > 0 || throw(ArgumentError("plot_xi_max must be positive"))
    0 <= spec.phase_fraction_min <= 1 || throw(ArgumentError("invalid phase_fraction_min"))
    spec.solver.k >= 2 || throw(ArgumentError("solver.k must be at least two"))
    return spec
end

function load_spec(path::AbstractString)
    config_path = abspath(path)
    isfile(config_path) || throw(ArgumentError("configuration file not found: $config_path"))
    config = TOML.parsefile(config_path)
    model = _section(config, :model)
    endpoint = _section(config, :endpoint_entanglement)
    output = String(_get(endpoint, :output, "output/endpoint_entanglement/n6_original_baseline"))
    spec = EndpointSpec(
        config_path=config_path,
        output=_project_path(output),
        nm1=Int(_get(model, :nm1, 6)),
        couplings=Couplings(_section(config, :hamiltonian)),
        solver=_solver(config),
        mu_left=Float64(_get(endpoint, :mu_left, 0.02)),
        mu_right=Float64(_get(endpoint, :mu_right, 0.22)),
        cut_x=Float64(_get(endpoint, :cut_x, 0.5)),
        qa=Int(_get(endpoint, :qa, 9)),
        f3a=Int(_get(endpoint, :f3a, 0)),
        f8a=Int(_get(endpoint, :f8a, 0)),
        lz2_values=Int.(_get(endpoint, :lz2_values, [-27, -25, -23, -21])),
        left_expected_counts=Int.(_get(endpoint, :left_expected_counts, [1, 3, 9, 22])),
        right_expected_counts=Int.(_get(endpoint, :right_expected_counts, [1, 1, 2, 3])),
        nm1_a=Int(_get(endpoint, :nm1_a, 3)),
        nm0_a=Int(_get(endpoint, :nm0_a, 8)),
        lambda_plot_min=Float64(_get(endpoint, :lambda_plot_min, 1e-14)),
        plot_xi_max=Float64(_get(endpoint, :plot_xi_max, 25.0)),
        phase_fraction_min=Float64(_get(endpoint, :phase_fraction_min, 0.95)),
        raw_config=config,
    )
    return _validate(spec)
end

function _run_id(spec::EndpointSpec)
    return MottJainED.stable_id(
        "endpoint-entanglement-v1", spec.nm1, MottJainED.coupling_vector(spec.couplings),
        spec.mu_left, spec.mu_right, spec.cut_x, spec.qa, spec.f3a, spec.f8a,
        spec.lz2_values, spec.left_expected_counts, spec.right_expected_counts,
        spec.nm1_a, spec.nm0_a, spec.lambda_plot_min,
    )
end

function plan(spec::EndpointSpec)
    nblocks = length(spec.lz2_values)
    return """
N=$(spec.nm1) endpoint entanglement plan
  Hamiltonian: Uf=$(spec.couplings.Uf), Uf0=$(spec.couplings.Uf0), U0=$(spec.couplings.U0), Vf=$(spec.couplings.Vf), Vf0=$(spec.couplings.Vf0), V0=$(spec.couplings.V0), t=$(spec.couplings.t)
  ground states: mu=$(spec.mu_left), $(spec.mu_right)
  fixed subsystem sector: QA=$(spec.qa), F3A=$(spec.f3a), F8A=$(spec.f8a)
  2LzA sectors: $(join(spec.lz2_values, ", "))
  left:  $nblocks sequential hemisphere RSES blocks, expected $(join(spec.left_expected_counts, ","))
  right: $nblocks sequential orbital ES blocks, expected $(join(spec.right_expected_counts, ","))
  output: $(spec.output)
"""
end

function _ground_state_audit(cache, mu::Float64, endpoint::String)
    candidates = NamedTuple[]
    best = nothing
    best_sector = nothing
    for sector in cache.sectors
        energies, vectors = MottJainED._eigensystem(sector, mu, cache.settings)
        order = sortperm(energies)
        energies = energies[order]
        vectors = vectors[:, order]
        for rank in eachindex(energies)
            push!(candidates, (
                endpoint=endpoint, mu=mu, z=sector.key.z, r=sector.key.r,
                rank=rank, energy=energies[rank],
            ))
        end
        index = argmin(energies)
        if isnothing(best) || energies[index] < best.energy
            best = (
                energy=energies[index], vector=copy(vectors[:, index]),
                basis=sector.basis, key=sector.key,
            )
            best_sector = sector
        end
    end
    isnothing(best) && error("no ground-state candidate at mu=$mu")
    table = DataFrame(candidates)
    table.energy_above_ground = table.energy .- best.energy
    nf = real(dot(best.vector, MottJainED.hermitian_opmat(best_sector.number_f) * best.vector))
    norm2 = real(dot(best.vector, best.vector))
    return merge(best, (Nf=nf, norm2=norm2,)), table
end

function _prepare_ground_cache(model, couplings, settings)
    # The standard spectrum cache also materializes L2 and C2.  Endpoint work
    # needs only H0 and Nf, so omitting those two large matrices substantially
    # lowers N=6 setup time and peak memory without changing the Hamiltonian.
    h0terms = hamiltonian_terms(model, couplings; include_mu=false)
    sectors = MottJainED.SectorCache[]
    for z in (1, -1), r in (1, -1)
        key = SectorKey(z, r)
        basis = FuzzifiED.Basis(model.cfs[0], [z, r], model.qnf)
        basis.dim == 0 && continue
        h0 = MottJainED.lower_sparse(MottJainED.float_opmat(
            FuzzifiED.Operator(basis, h0terms),
        ))
        number_f = MottJainED.lower_sparse(MottJainED.float_opmat(
            FuzzifiED.Operator(basis, model.number_f),
        ))
        push!(sectors, MottJainED.SectorCache(
            key, basis, h0, number_f, nothing, nothing, Float64[],
        ))
    end
    isempty(sectors) && error("no non-empty ground-state sectors")
    return ModelCache(model, couplings, settings, sectors)
end

function _save_ground_checkpoint(path, run_id, left, right)
    MottJainED.atomic_jldsave(
        path;
        run_id=run_id,
        left_vector=left.vector,
        left_energy=left.energy,
        left_z=left.key.z,
        left_r=left.key.r,
        left_Nf=left.Nf,
        left_norm2=left.norm2,
        right_vector=right.vector,
        right_energy=right.energy,
        right_z=right.key.z,
        right_r=right.key.r,
        right_Nf=right.Nf,
        right_norm2=right.norm2,
    )
    return path
end

function _basis_for_key(model, z::Int, r::Int)
    basis = FuzzifiED.Basis(model.cfs[0], [z, r], model.qnf)
    basis.dim > 0 || error("saved ground-state sector Z=$z,R=$r is empty")
    return basis
end

function _load_ground_checkpoint(spec, model, run_id)
    state_path = joinpath(spec.output, "ground_states.jld2")
    energy_path = joinpath(spec.output, "ground_sector_energies.csv")
    isfile(state_path) && isfile(energy_path) || return nothing
    saved = JLD2.load(state_path)
    get(saved, "run_id", "") == run_id || error("stale ground-state checkpoint")
    function restore(prefix)
        z = Int(saved["$(prefix)_z"])
        r = Int(saved["$(prefix)_r"])
        vector = Vector{Float64}(saved["$(prefix)_vector"])
        basis = _basis_for_key(model, z, r)
        length(vector) == basis.dim || error("saved $prefix vector has the wrong dimension")
        norm2 = real(dot(vector, vector))
        abs(norm2 - Float64(saved["$(prefix)_norm2"])) <= 1e-10 || error(
            "saved $prefix vector failed its norm check",
        )
        return (
            energy=Float64(saved["$(prefix)_energy"]), vector=vector,
            basis=basis, key=SectorKey(z, r), Nf=Float64(saved["$(prefix)_Nf"]),
            norm2=norm2,
        )
    end
    @info "reusing checkpointed endpoint ground states" path=state_path
    return (
        left=restore("left"), right=restore("right"),
        candidates=CSV.read(energy_path, DataFrame),
    )
end

function _orbital_cut(model, spec::EndpointSpec)
    orbitals_a = Int[]
    for m in 0:spec.nm1_a-1, flavor in 1:model.nf1
        push!(orbitals_a, m * model.nf1 + flavor)
    end
    append!(orbitals_a, model.no1 .+ collect(1:spec.nm0_a))
    set_a = Set(orbitals_a)
    orbitals_b = [orbital for orbital in 1:model.no if orbital ∉ set_a]
    amplitudes = ComplexF64[orbital in set_a ? 1 : 0 for orbital in 1:model.no]
    return orbitals_a, orbitals_b, amplitudes
end

function _hemisphere_amplitudes(model, x)
    return ComplexF64.(MottJainED._hemisphere_amplitudes(model, x))
end

function _block_path(spec::EndpointSpec, side::String, delta_l::Int)
    return joinpath(spec.output, "blocks", "$(side)_deltaL$(delta_l).csv")
end

function _compute_block(
    ground, model, qnd_a, qnd_b, amp_a, amp_b,
    sec_a::Vector{Int64}, sec_b::Vector{Int64};
    side::String, delta_l::Int, expected_count::Int, run_id::String,
)
    nor = ground.basis.cfs.nor
    cfsa = FuzzifiED.Confs(model.no, sec_a, qnd_a; nor=nor, num_th=1, disp_std=false)
    cfsb = FuzzifiED.Confs(model.no, sec_b, qnd_b; nor=nor, num_th=1, disp_std=false)
    cfsa.ncf > 0 || error("empty A configuration sector: $sec_a")
    cfsb.ncf > 0 || error("empty B configuration sector: $sec_b")
    qnf = FuzzifiED.QNOffd[]
    bsa = FuzzifiED.Basis(cfsa, ComplexF64[], qnf; num_th=1, disp_std=false)
    bsb = FuzzifiED.Basis(cfsb, ComplexF64[], qnf; num_th=1, disp_std=false)
    bsa.dim > 0 || error("empty A basis sector: $sec_a")
    bsb.dim > 0 || error("empty B basis sector: $sec_b")
    matrix_mib = 16.0 * bsa.dim * bsb.dim / 1024.0^2
    @info "entanglement block" side delta_l Lz2A=sec_a[2] dim_a=bsa.dim dim_b=bsb.dim matrix_mib
    decomposition = FuzzifiED.StateDecompMat(
        ground.vector, ground.basis, bsa, bsb, amp_a, amp_b,
    )
    singular_values = svdvals(decomposition)
    lambdas = sort!(Float64.(abs2.(singular_values)); rev=true)
    rows = NamedTuple[]
    for (rank, lambda) in enumerate(lambdas)
        lambda > 0 && isfinite(lambda) || continue
        push!(rows, (
            run_id=run_id, side=side, QA=sec_a[1], Lz2A=sec_a[2],
            F3A=sec_a[3], F8A=sec_a[4], delta_L=delta_l,
            rank_in_sector=rank, lambda=lambda, xi=-log(lambda),
            expected_count=expected_count, expected_branch_candidate=rank <= expected_count,
            dim_a=Int(bsa.dim), dim_b=Int(bsb.dim), matrix_mib=matrix_mib,
        ))
    end
    isempty(rows) && error("no nonzero Schmidt values for $side deltaL=$delta_l")
    return DataFrame(rows)
end

function _load_or_compute_block(
    spec, ground, model, qnd_a, qnd_b, amp_a, amp_b,
    sec_a, sec_b; side, delta_l, expected_count, run_id, force,
)
    path = _block_path(spec, side, delta_l)
    if isfile(path) && !force
        saved = CSV.read(path, DataFrame)
        required = ["run_id", "side", "Lz2A", "delta_L"]
        all(in(names(saved)), required) || error("invalid checkpoint: $path")
        all(saved.run_id .== run_id) || error("stale checkpoint identity: $path")
        all(saved.side .== side) || error("wrong checkpoint side: $path")
        all(saved.Lz2A .== sec_a[2]) || error("wrong checkpoint sector: $path")
        @info "reusing completed entanglement block" side delta_l path
        return saved
    end
    data = _compute_block(
        ground, model, qnd_a, qnd_b, amp_a, amp_b,
        sec_a, sec_b; side, delta_l, expected_count, run_id,
    )
    MottJainED.atomic_csv(path, data)
    data = nothing
    GC.gc()
    return CSV.read(path, DataFrame)
end

function _side_spectrum(spec, ground, model, side::String, expected_counts; force=false)
    total_charge = 3spec.nm1
    if side == "left_rses"
        qnd_a = model.qnd
        qnd_b = model.qnd
        amp_a = _hemisphere_amplitudes(model, spec.cut_x)
    elseif side == "right_oes"
        orbitals_a, orbitals_b, amp_a = _orbital_cut(model, spec)
        qnd_a = [model.qnd; FuzzifiED.GetPinOrbQNDiag(model.no, orbitals_b)]
        qnd_b = [model.qnd; FuzzifiED.GetPinOrbQNDiag(model.no, orbitals_a)]
    else
        throw(ArgumentError("unknown side: $side"))
    end
    amp_b = sqrt.(max.(0.0, 1 .- abs2.(amp_a)))
    run_id = _run_id(spec)
    blocks = DataFrame[]
    for (delta_l, (lz2, expected)) in enumerate(zip(spec.lz2_values, expected_counts))
        delta_l -= 1
        sec_a = Int64[spec.qa, lz2, spec.f3a, spec.f8a]
        sec_b = Int64[total_charge - spec.qa, -lz2, -spec.f3a, -spec.f8a]
        if side == "right_oes"
            push!(sec_a, 0)
            push!(sec_b, 0)
        end
        push!(blocks, _load_or_compute_block(
            spec, ground, model, qnd_a, qnd_b, amp_a, amp_b, sec_a, sec_b;
            side, delta_l, expected_count=expected, run_id, force,
        ))
    end
    data = vcat(blocks...; cols=:union)
    data.xi_shifted = data.xi .- minimum(data.xi)
    sort!(data, [:delta_L, :rank_in_sector])
    return data
end

function _diagnostics(data::DataFrame)
    rows = NamedTuple[]
    for block in groupby(data, :delta_L)
        sort!(block, :rank_in_sector)
        expected = Int(first(block.expected_count))
        available = nrow(block)
        top = available >= expected ? block.xi_shifted[expected] : NaN
        next = available > expected ? block.xi_shifted[expected + 1] : NaN
        push!(rows, (
            side=String(first(block.side)), delta_L=Int(first(block.delta_L)),
            Lz2A=Int(first(block.Lz2A)), expected_count=expected,
            levels_available=available, dim_a=Int(first(block.dim_a)),
            dim_b=Int(first(block.dim_b)), matrix_mib=Float64(first(block.matrix_mib)),
            target_block_weight=sum(block.lambda), xi_expected_top=top, xi_next=next,
            gap_after_expected=isfinite(next) && isfinite(top) ? next - top : NaN,
        ))
    end
    return sort!(DataFrame(rows), :delta_L)
end

function _plot_panel!(axis, data, spec, color, title)
    visible = filter(row -> row.lambda >= spec.lambda_plot_min &&
        row.xi_shifted <= spec.plot_xi_max, data)
    background = filter(row -> !row.expected_branch_candidate, visible)
    candidate = filter(row -> row.expected_branch_candidate, visible)
    nrow(background) > 0 && scatter!(
        axis, background.delta_L, background.xi_shifted;
        color=(:gray, 0.42), markersize=5, label="other levels",
    )
    nrow(candidate) > 0 && scatter!(
        axis, candidate.delta_L, candidate.xi_shifted;
        color=color, markersize=10, label="lowest expected-count levels",
    )
    axis.title = title
    axis.xlabel = "ΔL"
    axis.ylabel = "ξ − ξₘᵢₙ"
    axis.xticks = (collect(0:length(spec.lz2_values)-1), string.(0:length(spec.lz2_values)-1))
    ylims!(axis, 0, spec.plot_xi_max)
    axislegend(axis; position=:rt, framevisible=false)
    return axis
end

function _save_plots(left, right, spec)
    left_counts = join(spec.left_expected_counts, ",")
    right_counts = join(spec.right_expected_counts, ",")
    left_fig = Figure(size=(700, 520))
    _plot_panel!(
        Axis(left_fig[1, 1]), left, spec, :dodgerblue,
        "left endpoint: hemisphere RSES (target $left_counts)",
    )
    save(joinpath(spec.output, "left_fiqh_rses.png"), left_fig)

    right_fig = Figure(size=(700, 520))
    _plot_panel!(
        Axis(right_fig[1, 1]), right, spec, :darkorange,
        "right endpoint: orbital ES (target $right_counts)",
    )
    save(joinpath(spec.output, "right_laughlin_oes.png"), right_fig)

    combined = Figure(size=(1320, 520))
    _plot_panel!(
        Axis(combined[1, 1]), left, spec, :dodgerblue,
        "left μ=$(spec.mu_left): RSES ($left_counts)",
    )
    _plot_panel!(
        Axis(combined[1, 2]), right, spec, :darkorange,
        "right μ=$(spec.mu_right): OES ($right_counts)",
    )
    save(joinpath(spec.output, "endpoint_entanglement_comparison.png"), combined)
    return nothing
end

function _ground_summary(spec, endpoint, cut, mu, ground)
    n0 = (3spec.nm1 - ground.Nf) / 3
    phase_fraction = endpoint == "left_fIQH" ? ground.Nf / (3spec.nm1) : n0 / spec.nm1
    return (
        endpoint=endpoint, cut=cut, mu=mu, ground_energy=ground.energy,
        ground_z=ground.key.z, ground_r=ground.key.r, state_norm2=ground.norm2,
        Nf=ground.Nf, N0=n0, target_phase_fraction=phase_fraction,
        phase_fraction_min=spec.phase_fraction_min,
        phase_fraction_pass=phase_fraction >= spec.phase_fraction_min,
    )
end

"""
Compute only the four audited edge sectors at each endpoint.

The left endpoint uses a hemisphere real-space cut; the right endpoint uses an
orbital cut.  Blocks are checkpointed separately, so a restarted job does not
repeat completed dense SVDs.  The highlighted levels are the lowest `K` levels
for the *hypothesized* counting, while `block_diagnostics.csv` records the gap
after level K; highlighting is not treated as proof of the counting.
"""
function run_endpoint_entanglement(spec::EndpointSpec; force::Bool=false)
    _validate(spec)
    FuzzifiED.NumThreads = Threads.nthreads()
    BLAS.set_num_threads(1)
    MottJainED.ensure_output(spec.output)
    MottJainED.ensure_output(joinpath(spec.output, "blocks"))
    run_id = _run_id(spec)
    identity_path = joinpath(spec.output, "case_identity.toml")
    if isfile(identity_path)
        saved = TOML.parsefile(identity_path)
        get(saved, "run_id", "") == run_id || error(
            "output directory belongs to a different endpoint configuration: $(spec.output)",
        )
    else
        MottJainED.atomic_toml(identity_path, Dict(
            "run_id" => run_id, "created_at" => string(now()),
            "config_path" => spec.config_path,
        ))
    end
    MottJainED.atomic_toml(joinpath(spec.output, "resolved_config.toml"), spec.raw_config)
    MottJainED.write_run_metadata(
        spec.output; command="endpoint-entanglement", config_path=spec.config_path,
    )

    model = build_model(nm1=spec.nm1)
    checkpoint = force ? nothing : _load_ground_checkpoint(spec, model, run_id)
    cache = nothing
    if isnothing(checkpoint)
        @info "building lean N=$(spec.nm1) Hamiltonian cache (H0 and Nf only)"
        cache = _prepare_ground_cache(model, spec.couplings, spec.solver)
        @info "solving left endpoint ground state" mu=spec.mu_left
        left_ground, left_candidates = _ground_state_audit(cache, spec.mu_left, "left_fIQH")
        @info "solving right endpoint ground state" mu=spec.mu_right
        right_ground, right_candidates = _ground_state_audit(cache, spec.mu_right, "right_Laughlin")
        candidate_table = vcat(left_candidates, right_candidates)
        MottJainED.atomic_csv(joinpath(spec.output, "ground_sector_energies.csv"), candidate_table)
        _save_ground_checkpoint(
            joinpath(spec.output, "ground_states.jld2"), run_id, left_ground, right_ground,
        )
    else
        left_ground = checkpoint.left
        right_ground = checkpoint.right
        candidate_table = checkpoint.candidates
    end
    ground_summary = DataFrame([
        _ground_summary(spec, "left_fIQH", "hemisphere_RSES", spec.mu_left, left_ground),
        _ground_summary(spec, "right_Laughlin", "orbital_ES", spec.mu_right, right_ground),
    ])
    MottJainED.atomic_csv(joinpath(spec.output, "ground_summary.csv"), ground_summary)
    for row in eachrow(ground_summary)
        row.phase_fraction_pass || @warn(
            "endpoint is not deep enough by the configured density check",
            endpoint=row.endpoint, fraction=row.target_phase_fraction,
            required=row.phase_fraction_min,
        )
    end

    # The sparse Hamiltonian cache is the largest resident object.  Entanglement
    # only needs each selected vector and its basis, so release all operator matrices.
    cache = nothing
    GC.gc()

    left = _side_spectrum(
        spec, left_ground, model, "left_rses", spec.left_expected_counts; force,
    )
    MottJainED.atomic_csv(joinpath(spec.output, "left_fiqh_rses_spectrum.csv"), left)
    left_diag = _diagnostics(left)
    MottJainED.atomic_csv(joinpath(spec.output, "left_fiqh_rses_diagnostics.csv"), left_diag)
    left_ground = nothing
    GC.gc()

    right = _side_spectrum(
        spec, right_ground, model, "right_oes", spec.right_expected_counts; force,
    )
    MottJainED.atomic_csv(joinpath(spec.output, "right_laughlin_oes_spectrum.csv"), right)
    right_diag = _diagnostics(right)
    MottJainED.atomic_csv(joinpath(spec.output, "right_laughlin_oes_diagnostics.csv"), right_diag)
    right_ground = nothing
    GC.gc()

    diagnostics = vcat(left_diag, right_diag)
    MottJainED.atomic_csv(joinpath(spec.output, "block_diagnostics.csv"), diagnostics)
    _save_plots(left, right, spec)
    MottJainED.atomic_toml(joinpath(spec.output, "completed.toml"), Dict(
        "run_id" => run_id,
        "completed_at" => string(now()),
        "phase_checks_pass" => all(ground_summary.phase_fraction_pass),
        "left_target_sector_weight" => sum(left.lambda),
        "right_target_sector_weight" => sum(right.lambda),
    ))
    @info "endpoint entanglement workflow complete" output=spec.output
    return (
        ground=ground_summary, left=left, right=right,
        diagnostics=diagnostics, output=spec.output,
    )
end

end # module
