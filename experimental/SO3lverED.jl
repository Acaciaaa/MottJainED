module SO3lverED

using FuzzifiED
using FuzzifiED.SO3lver
using LinearAlgebra
using MottJainED

const COMPONENT_ORDER = (:Uf, :Uf0, :U0, :Vf, :Vf0, :V0, :t, :mu)
const CFT_BLOCK_KEYS = (
    (:singlet, 0), (:singlet, 1), (:singlet, 2),
    (:adjoint, 0), (:adjoint, 1), (:adjoint, 2),
)
const CFT_SCORE_LABELS = ("dS-S", "J", "curlJ", "dJ(rank1)", "T(rank1)")
const CFT_SCORE_TARGETS = Float64[1, 2, 3, 3, 3]
const REPRESENTATIONS = Dict(
    :singlet => (c2=0.0, f3=0, f8=0),
    # The adjoint highest weight occurs once per SU(3) octet.  This avoids the
    # two zero-weight copies present in the old Cartan-zero calculation.
    :adjoint => (c2=3.0, f3=1, f8=3),
)

"""Size- and representation-dependent data used by the SO(3)-adapted solver."""
struct SO3Model
    nm1::Int
    nm0::Int
    no1::Int
    no0::Int
    representation::Symbol
    c2::Float64
    weight::NTuple{2,Int}
    qnd_light::Vector{QNDiag}
    qnd_heavy::Vector{QNDiag}
    sec_light::Matrix{Int64}
    sec_heavy::Matrix{Int64}
    tms_lzlp_light::Tuple{Terms,Terms}
    tms_lzlp_heavy::Tuple{Terms,Terms}
    tms_c2::Terms
    components::NamedTuple
    decompositions::NamedTuple
    all_decompositions::CoupleDecomps
    channel_component::Vector{Symbol}
    channel_coefficient::Vector{ComplexF64}
end

"""Reusable segment spaces and reduced matrix elements for one SU(3) irrep."""
struct SO3Workspace
    model::SO3Model
    light_space::SegSpace{Float64}
    heavy_space::SegSpace{Float64}
    segment_operators::Matrix{SegOperator}
    build_seconds::Float64
end

"""A fixed-(L, SU(3)) matrix-free Hamiltonian whose couplings can be retuned."""
mutable struct SO3Hamiltonian
    workspace::SO3Workspace
    ell::Int
    space::CompSpace{Float64}
    operator::CompOperator{Float64}
    couplings::Couplings
end

representation_data(representation::Symbol) = get(REPRESENTATIONS, representation) do
    throw(ArgumentError("representation must be :singlet or :adjoint"))
end

function _flavor_populations(ne::Int, f3::Int, f8::Int)
    # Use explicit multiplication: Julia parses `3f3` as a Float32 literal.
    numerators = (
        2 * ne + 3 * f3 + f8,
        2 * ne - 3 * f3 + f8,
        2 * ne - 2 * f8,
    )
    all(value -> value % 6 == 0, numerators) || return nothing
    return ntuple(i -> numerators[i] ÷ 6, 3)
end

function _light_sectors(nm1::Int, f3::Int, f8::Int)
    sectors = Vector{Int64}[]
    for ne in 0:3:3nm1
        populations = _flavor_populations(ne, f3, f8)
        isnothing(populations) && continue
        all(n -> 0 <= n <= nm1, populations) || continue
        # A segment uses m=0 for integer L and m=1/2 for half-integer L.
        lz2 = ((nm1 + 1) * ne) % 2
        push!(sectors, Int64[ne, lz2, f3, f8])
    end
    isempty(sectors) && error("No light-fermion sectors for weight ($f3,$f8)")
    return reduce(hcat, sectors)
end

function _heavy_sectors(nm0::Int, nm1::Int)
    sectors = [
        Int64[3ne, ((nm0 + 1) * ne) % 2, 0, 0]
        for ne in 0:nm1
    ]
    return reduce(hcat, sectors)
end

function _combine_decompositions(decompositions::NamedTuple)
    all_decompositions = CoupleDecomp[]
    channel_component = Symbol[]
    channel_coefficient = ComplexF64[]
    for component in COMPONENT_ORDER
        for decomposition in getproperty(decompositions, component)
            # Laplacian observables have an identically zero l=0 component.
            # ContactCouple currently still emits that formal channel; passing
            # its empty Terms to FuzzifiED.Operator raises an empty-reduction
            # error.  Removing it is exact because the whole tensor-product
            # channel is zero.
            all(eachindex(decomposition.amd)) do segment
                modes = decomposition.amd[segment]
                modes === :Identity && return true
                rank = decomposition.ch[1, segment]
                return any(
                    angular_rank == rank && !isempty(terms)
                    for ((angular_rank, _), terms) in modes.comps
                )
            end || continue
            push!(all_decompositions, decomposition)
            push!(channel_component, component)
            push!(channel_coefficient, decomposition.coeff)
        end
    end
    return all_decompositions, channel_component, channel_coefficient
end

"""
Construct the exact same eight Hamiltonian components as `MottJainED.build_model`,
but decomposed into charge-1 and charge-3 segments for SO(3)lver.
"""
function build_so3_model(; nm1::Int, representation::Symbol=:singlet)
    nm1 >= 2 || throw(ArgumentError("nm1 must be at least 2"))
    irrep = representation_data(representation)
    nm0 = 3nm1 - 2
    no1 = 3nm1
    no0 = nm0
    radius2 = Float64(nm1)
    FuzzifiED.ElementType = Float64
    FuzzifiED.ObsNormRadSq = radius2
    FuzzifiED.ObsMomIncr = true

    qnd_light = QNDiag[
        GetNeQNDiag(no1),
        GetLz2QNDiag(nm1, 3),
        GetFlavQNDiag(nm1, 3, [1, -1, 0]),
        GetFlavQNDiag(nm1, 3, [1, 1, -2]),
    ]
    qnd_heavy = QNDiag[
        3GetNeQNDiag(no0),
        GetLz2QNDiag(nm0, 1),
        zero(QNDiag, no0),
        zero(QNDiag, no0),
    ]
    sec_light = _light_sectors(nm1, irrep.f3, irrep.f8)
    sec_heavy = _heavy_sectors(nm0, nm1)
    tms_lzlp_light = GetLzLpTerms(nm1, 3)
    tms_lzlp_heavy = GetLzLpTerms(nm0, 1)
    tms_c2 = GetC2Terms(nm1, 3, :SU)

    fields = [GetElectronObs(nm1, 3, flavor; norm_r2=radius2) for flavor in 1:3]
    heavy_field = GetElectronObs(nm0, 1, 1; norm_r2=radius2)
    density_light = GetDensityObs(nm1, 3; norm_r2=radius2)
    density_heavy = heavy_field' * heavy_field
    triplet = fields[1] * fields[2] * fields[3]
    lap_light = Laplacian(density_light; norm_r2=radius2)
    lap_heavy = Laplacian(density_heavy; norm_r2=radius2)

    # Keeping the same-segment Terms is useful for exact regression tests
    # against the conventional Fock-basis implementation.
    components = (
        Uf=SimplifyTerms(GetIntegral(density_light * density_light; norm_r2=radius2)),
        Uf0=nothing,
        U0=SimplifyTerms(GetIntegral(density_heavy * density_heavy; norm_r2=radius2)),
        Vf=SimplifyTerms(GetIntegral(density_light * lap_light; norm_r2=radius2)),
        Vf0=nothing,
        V0=SimplifyTerms(GetIntegral(density_heavy * lap_heavy; norm_r2=radius2)),
        t=nothing,
        mu=GetPolTerms(nm1, 3),
    )

    zero_shift = zeros(Int64, 4)
    no_shift = zeros(Int64, 4, 2)
    light_to_heavy = zeros(Int64, 4, 2)
    light_to_heavy[1, :] .= (-3, 3)
    heavy_to_light = -light_to_heavy
    decompositions = (
        Uf=SingleSegCouple(2, 1, components.Uf, zero_shift),
        Uf0=ContactCouple([density_light, density_heavy], no_shift),
        U0=SingleSegCouple(2, 2, components.U0, zero_shift),
        Vf=SingleSegCouple(2, 1, components.Vf, zero_shift),
        Vf0=ContactCouple([lap_light, density_heavy], no_shift),
        V0=SingleSegCouple(2, 2, components.V0, zero_shift),
        # Both factors are fermion-odd.  The graded tensor product therefore
        # contributes a relative minus sign between a conversion channel and
        # its Hermitian conjugate.
        t=-ContactCouple([triplet, heavy_field'], light_to_heavy) +
            ContactCouple([triplet', heavy_field], heavy_to_light),
        mu=SingleSegCouple(2, 1, components.mu, zero_shift),
    )
    all_decompositions, channel_component, channel_coefficient =
        _combine_decompositions(decompositions)

    return SO3Model(
        nm1, nm0, no1, no0, representation, irrep.c2,
        (irrep.f3, irrep.f8), qnd_light, qnd_heavy,
        sec_light, sec_heavy, tms_lzlp_light, tms_lzlp_heavy,
        tms_c2, components, decompositions, all_decompositions,
        channel_component, channel_coefficient,
    )
end

"""Build the expensive segment decomposition once and reuse it across L and couplings."""
function build_workspace(
    model::SO3Model;
    l2c2_ratio::Float64=0.1 / model.nm1^2,
    nst_max_light::Vector{Int}=zeros(Int, size(model.sec_light, 2)),
    nst_max_heavy::Vector{Int}=zeros(Int, size(model.sec_heavy, 2)),
    heavy_space::Union{Nothing,SegSpace{Float64}}=nothing,
    disp_std::Bool=true,
)
    started = time()
    light_space = BuildSegSpace(
        model.no1, model.sec_light, model.qnd_light,
        model.tms_lzlp_light, model.tms_c2, [model.c2];
        l2c2_ratio, nst_max=Int64.(nst_max_light), disp_std,
    )
    if isnothing(heavy_space)
        heavy_space = BuildSegSpace(
            model.no0, model.sec_heavy, model.qnd_heavy,
            model.tms_lzlp_heavy;
            nst_max=Int64.(nst_max_heavy), disp_std,
        )
    else
        heavy_space.sec == model.sec_heavy || throw(ArgumentError(
            "reused heavy segment has incompatible quantum-number sectors",
        ))
    end
    segment_operators = BuildSegOperators(
        [light_space, heavy_space], model.all_decompositions; disp_std,
    )
    return SO3Workspace(
        model, light_space, heavy_space, segment_operators, time() - started,
    )
end

_coupling(couplings::Couplings, component::Symbol) = getfield(couplings, component)

function _retune_operator!(operator::CompOperator{Float64}, model::SO3Model, couplings::Couplings)
    MottJainED.validate(couplings)
    for channel in eachindex(operator.coeff)
        value = model.channel_coefficient[channel] *
                _coupling(couplings, model.channel_component[channel])
        abs(imag(value)) <= 1e-10 || error(
            "SO(3)lver channel $channel is unexpectedly complex: $value",
        )
        operator.coeff[channel] = real(value)
    end
    return operator
end

"""Build a matrix-free Hamiltonian in a fixed physical angular-momentum sector."""
function build_hamiltonian(
    workspace::SO3Workspace,
    ell::Int,
    couplings::Couplings;
    disp_std::Bool=true,
)
    ell >= 0 || throw(ArgumentError("ell must be non-negative"))
    model = workspace.model
    total_sector = Int64[3model.nm1, 0, model.weight[1], model.weight[2]]
    space = BuildCompSpace(
        [workspace.light_space, workspace.heavy_space], total_sector, 2ell;
        disp_std,
    )
    # Use the full method because the convenience wrapper in FuzzifiED 2.0.1
    # does not forward `disp_std`.
    operator = BuildCompOperator(
        space, space, model.all_decompositions,
        workspace.segment_operators, 0; disp_std,
    )
    _retune_operator!(operator, model, couplings)
    return SO3Hamiltonian(workspace, ell, space, operator, couplings)
end

function retune!(hamiltonian::SO3Hamiltonian, couplings::Couplings)
    _retune_operator!(hamiltonian.operator, hamiltonian.workspace.model, couplings)
    hamiltonian.couplings = couplings
    return hamiltonian
end

"""Return the lowest energies (and optionally vectors) in one exact (L,C2) block."""
function solve(
    hamiltonian::SO3Hamiltonian;
    k::Int=10,
    tol::Float64=1e-8,
    ncv::Int=max(2k, k + 10),
    vectors::Bool=false,
    initvec::Union{Nothing,AbstractVector{<:Real}}=nothing,
    dense_cutoff::Int=128,
    disp_std::Bool=true,
)
    dimension = hamiltonian.space.dim
    dimension > 0 || return vectors ? (Float64[], zeros(Float64, 0, 0)) : Float64[]
    count = min(k, dimension)
    if dimension <= dense_cutoff || count == dimension
        decomposition = eigen(Symmetric(Matrix(hamiltonian.operator; disp_std=false)))
        order = sortperm(decomposition.values)[1:count]
        energies = Float64.(decomposition.values[order])
        states = Matrix{Float64}(decomposition.vectors[:, order])
    else
        count = min(count, dimension - 1)
        krylov_dimension = min(dimension, max(count + 1, ncv))
        kwargs = (; tol, ncv=krylov_dimension, issymmetric=true, disp_std)
        if isnothing(initvec)
            energies, states = GetEigensystem(
                hamiltonian.operator, count; kwargs...,
            )
        else
            length(initvec) == dimension || throw(DimensionMismatch(
                "initial vector has length $(length(initvec)); expected $dimension",
            ))
            energies, states = GetEigensystem(
                hamiltonian.operator, count;
                initvec=Float64.(initvec), kwargs...,
            )
        end
        # KrylovKit may return an extra converged Ritz value; keep the public
        # contract exact and deterministic.
        order = sortperm(energies)[1:count]
        energies = Float64.(energies[order])
        states = Matrix{Float64}(states[:, order])
    end
    return vectors ? (energies, states) : energies
end

sector_dimension(hamiltonian::SO3Hamiltonian) = hamiltonian.space.dim

function _block_energies(blocks, representation::Symbol, ell::Int, count::Int)
    key = (representation, ell)
    haskey(blocks, key) || throw(ArgumentError(
        "missing SO(3)lver spectrum block representation=$representation L=$ell",
    ))
    energies = sort!(Float64.(collect(blocks[key])))
    length(energies) >= count || throw(ArgumentError(
        "SO(3)lver spectrum block representation=$representation L=$ell " *
        "contains $(length(energies)) levels; need at least $count",
    ))
    all(isfinite, energies) || throw(ArgumentError(
        "non-finite energy in representation=$representation L=$ell",
    ))
    return energies
end

"""
Evaluate the project's fixed five-relation CFT locator from exact physical blocks.

The old Fock-basis workflow selected raw ranks before merging the two Cartan-zero
copies of every SU(3) octet.  An adjoint highest-weight block contains each octet
once, so old raw rank 3 for curl-J becomes physical rank 2 here.  The other five
relations and the fitted scale factor are unchanged.
"""
function score_cft_blocks(blocks)
    singlet_l0 = _block_energies(blocks, :singlet, 0, 2)
    singlet_l1 = _block_energies(blocks, :singlet, 1, 1)
    singlet_l2 = _block_energies(blocks, :singlet, 2, 1)
    adjoint_l0 = _block_energies(blocks, :adjoint, 0, 1)
    adjoint_l1 = _block_energies(blocks, :adjoint, 1, 2)
    adjoint_l2 = _block_energies(blocks, :adjoint, 2, 1)

    block_minima = [
        (representation=rep, ell=ell,
         energy=first(_block_energies(blocks, rep, ell, 1)))
        for (rep, ell) in CFT_BLOCK_KEYS
    ]
    ground_block = block_minima[argmin(getproperty.(block_minima, :energy))]
    ground = ground_block.energy

    raw_gaps = Float64[
        singlet_l1[1] - singlet_l0[2],
        adjoint_l1[1] - ground,
        adjoint_l1[2] - ground,
        adjoint_l2[1] - ground,
        singlet_l2[1] - ground,
    ]
    factor = dot(raw_gaps, CFT_SCORE_TARGETS) / dot(CFT_SCORE_TARGETS, CFT_SCORE_TARGETS)
    isfinite(factor) && factor > 0 || error("fitted CFT scale factor is not positive")
    scaled_gaps = raw_gaps ./ factor
    residuals = scaled_gaps .- CFT_SCORE_TARGETS
    q = sqrt(sum(abs2, residuals) / length(residuals))
    return (
        q=q,
        factor=factor,
        delta_s=(singlet_l0[2] - ground) / factor,
        delta_o=(adjoint_l0[1] - ground) / factor,
        raw_gaps=raw_gaps,
        target_gaps=copy(CFT_SCORE_TARGETS),
        scaled_gaps=scaled_gaps,
        residuals=residuals,
        labels=collect(CFT_SCORE_LABELS),
        ground_energy=ground,
        ground_representation=ground_block.representation,
        ground_ell=ground_block.ell,
        ground_is_singlet_l0=(ground_block.representation == :singlet && ground_block.ell == 0),
    )
end

export SO3Model, SO3Workspace, SO3Hamiltonian,
       build_so3_model, build_workspace, build_hamiltonian,
       retune!, solve, sector_dimension, representation_data,
       score_cft_blocks, CFT_BLOCK_KEYS

end
