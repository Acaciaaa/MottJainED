module SO3lverED

using FuzzifiED
using FuzzifiED.JackToolkit
using FuzzifiED.SO3lver
using LinearAlgebra
using MottJainED
using Optim

const COMPONENT_ORDER = (:Uf, :Uf0, :U0, :Vf, :Vf0, :V0, :t, :mu)
const CFT_BLOCK_KEYS = (
    (:singlet, 0), (:singlet, 1), (:singlet, 2),
    (:adjoint, 0), (:adjoint, 1), (:adjoint, 2),
)
const CFT_SCORE_TERMS = (:ds_s, :j, :curlj, :dj_rank1, :t_rank1)
const CFT_STABLE_SIX_TERMS = (
    :ds_s, :dds_ds, :j, :curlj, :dj_rank1, :t_rank1,
)
const CFT_AUDITED_SEVEN_TERMS = (
    :ds_s, :dds_ds, :boxs_s, :j, :curlj, :dj_rank1, :t_rank1,
)
const CFT_RELATION_SPECS = Dict{Symbol,NamedTuple}(
    :ds_s => (
        label="dS-S", upper=(:singlet, 1, 1), lower=(:singlet, 0, 2), target=1.0,
    ),
    :dds_ds => (
        label="ddS-dS", upper=(:singlet, 2, 2), lower=(:singlet, 1, 1), target=1.0,
    ),
    :boxs_s => (
        label="boxS-S", upper=(:singlet, 0, 3), lower=(:singlet, 0, 2), target=2.0,
    ),
    :boxo_o => (
        label="boxO-O", upper=(:adjoint, 0, 2), lower=(:adjoint, 0, 1), target=2.0,
    ),
    :j => (
        label="J", upper=(:adjoint, 1, 1), lower=nothing, target=2.0,
    ),
    :curlj => (
        label="curlJ", upper=(:adjoint, 1, 2), lower=nothing, target=3.0,
    ),
    :boxj_j => (
        label="boxJ-J", upper=(:adjoint, 1, 3), lower=(:adjoint, 1, 1), target=2.0,
    ),
    :dj_rank1 => (
        label="dJ(rank1)", upper=(:adjoint, 2, 1), lower=nothing, target=3.0,
    ),
    :t_rank1 => (
        label="T(rank1)", upper=(:singlet, 2, 1), lower=nothing, target=3.0,
    ),
)
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
    heavy_space_mode::Symbol
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

const HEAVY_SPACE_MODES = (:full, :laughlin13)

function _check_heavy_space_mode(mode::Symbol)
    mode in HEAVY_SPACE_MODES || throw(ArgumentError(
        "heavy_space_mode must be :full or :laughlin13",
    ))
    return mode
end

"""The positive fermionic `V₁` parent Hamiltonian for the heavy segment."""
heavy_v1_terms(model::SO3Model) = GetDenIntTerms(model.nm0, 1, [0.0, 1.0])

function _laughlin13_jack(
    no::Int64,
    sector::Vector{Int64},
    basis::Basis,
    l2v1_matrix::OpMat,
    requested_count::Int64,
)
    particle_count = sector[1] ÷ 3
    states = GetJackStates(
        basis, no, particle_count, 1, 3, sector[2],
    )
    size(states, 2) == requested_count || error(
        "Laughlin 1/3 Jack count mismatch in sector $sector: " *
        "constructed $(size(states, 2)), requested $requested_count",
    )
    return OrganizeJackStates(states, l2v1_matrix)
end

function laughlin13_root_counts(model::SO3Model)
    return Int64[
        length(GetJackRoots(
            model.nm0, model.sec_heavy[1, sector] ÷ 3,
            1, 3, model.sec_heavy[2, sector],
        ))
        for sector in axes(model.sec_heavy, 2)
    ]
end

"""
Build the exact fermionic Laughlin-1/3 quasihole space of the heavy segment.

The `V₁` parent Hamiltonian is supplied as the auxiliary Casimir.  This is
important even though Jack states are used in the iterative branch: for small
sectors SO3lver deliberately switches to dense diagonalization, and the
zero-Casimir filter then still removes every non-quasihole state.
"""
function build_laughlin13_heavy_space(
    model::SO3Model;
    l2v1_ratio::Float64=0.1 / model.nm1^2,
    disp_std::Bool=true,
)
    counts = laughlin13_root_counts(model)
    all(>(0), counts) || error(
        "a required heavy sector has no Laughlin-1/3 quasihole state",
    )
    return BuildSegSpace(
        model.no0, model.sec_heavy, model.qnd_heavy,
        model.tms_lzlp_heavy, heavy_v1_terms(model), [0.0];
        l2c2_ratio=l2v1_ratio,
        nst_max=counts,
        diag_method=_laughlin13_jack,
        disp_std,
    )
end

"""Build the expensive segment decomposition once and reuse it across L and couplings."""
function build_workspace(
    model::SO3Model;
    l2c2_ratio::Float64=0.1 / model.nm1^2,
    nst_max_light::Vector{Int}=zeros(Int, size(model.sec_light, 2)),
    nst_max_heavy::Vector{Int}=zeros(Int, size(model.sec_heavy, 2)),
    heavy_space::Union{Nothing,SegSpace{Float64}}=nothing,
    heavy_space_mode::Symbol=:full,
    disp_std::Bool=true,
)
    _check_heavy_space_mode(heavy_space_mode)
    heavy_space_mode == :laughlin13 && any(x -> !iszero(x), nst_max_heavy) &&
        throw(ArgumentError(
            "nst_max_heavy is determined exactly by the Laughlin root count",
        ))
    started = time()
    light_space = BuildSegSpace(
        model.no1, model.sec_light, model.qnd_light,
        model.tms_lzlp_light, model.tms_c2, [model.c2];
        l2c2_ratio, nst_max=Int64.(nst_max_light), disp_std,
    )
    if isnothing(heavy_space)
        if heavy_space_mode == :full
            heavy_space = BuildSegSpace(
                model.no0, model.sec_heavy, model.qnd_heavy,
                model.tms_lzlp_heavy;
                nst_max=Int64.(nst_max_heavy), disp_std,
            )
        else
            heavy_space = build_laughlin13_heavy_space(model; disp_std)
        end
    else
        heavy_space.sec == model.sec_heavy || throw(ArgumentError(
            "reused heavy segment has incompatible quantum-number sectors",
        ))
        if heavy_space_mode == :laughlin13
            size.(heavy_space.sts, 2) == laughlin13_root_counts(model) ||
                throw(ArgumentError(
                    "reused heavy segment does not have the exact " *
                    "Laughlin-1/3 root counts",
                ))
        end
    end
    segment_operators = BuildSegOperators(
        [light_space, heavy_space], model.all_decompositions; disp_std,
    )
    return SO3Workspace(
        model, light_space, heavy_space, heavy_space_mode,
        segment_operators, time() - started,
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

"""
Build the heavy-segment `V₁` parent Hamiltonian in an existing composite
block.  It is identically zero (up to roundoff) in a `:laughlin13` workspace
and is also useful for explicit `P H P` and large-`V₁` validation tests.
"""
function build_heavy_v1_operator(
    hamiltonian::SO3Hamiltonian;
    disp_std::Bool=true,
)
    model = hamiltonian.workspace.model
    decomposition = SingleSegCouple(
        2, 2, heavy_v1_terms(model), zeros(Int64, 4),
    )
    spaces = [
        hamiltonian.workspace.light_space,
        hamiltonian.workspace.heavy_space,
    ]
    segment_operators = BuildSegOperators(
        spaces, decomposition; disp_std,
    )
    return BuildCompOperator(
        hamiltonian.space, hamiltonian.space,
        decomposition, segment_operators, 0; disp_std,
    )
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

const GENERATOR_CANDIDATE_NAMES = (
    :Uf, :Vf, :U0, :V0, :Uf0, :Vf0, :t, :mu,
)
const DEFAULT_CONFORMAL_PRIMARY_SPECS = (
    (label=:S, representation=:singlet, ell=0, rank=2),
    (label=:O, representation=:adjoint, ell=0, rank=1),
    (label=:J, representation=:adjoint, ell=1, rank=1),
    (label=:T, representation=:singlet, ell=2, rank=1),
)

"""The eight SU(3)-singlet microscopic rank-one tensors used to fit `Λ=P+K`."""
struct SO3GeneratorCandidates
    names::NTuple{8,Symbol}
    decompositions::Vector{CoupleDecomps}
end

"""The same generator candidates represented between two exact SO(3) blocks."""
struct SO3GeneratorOperators
    names::NTuple{8,Symbol}
    initial_space::CompSpace{Float64}
    final_space::CompSpace{Float64}
    operators::Vector{CompOperator{Float64}}
end

function _real_generator_decomposition(cpd::CoupleDecomps, name::Symbol)
    reference_index = findfirst(channel -> abs(channel.coeff) > 1.0e-13, cpd)
    isnothing(reference_index) && error("generator candidate $name is identically zero")
    reference = cpd[reference_index].coeff
    # Odd-rank Hermitian tensors can carry one common imaginary phase in the
    # spherical-tensor convention.  SO3lver's Float64 spaces require removing
    # that convention-dependent phase from the whole candidate.
    phase = abs(imag(reference)) <= 1.0e-12 ? 1.0 + 0.0im :
            conj(reference / abs(reference))
    output = phase * cpd
    for channel in output
        value = channel.coeff
        abs(imag(value)) <= 1.0e-11 * max(1.0, abs(real(value))) || error(
            "generator candidate $name has incompatible channel phases",
        )
        channel.coeff = ComplexF64(real(value))
    end
    return output
end

"""
Construct the microscopic conformal-generator ansatz as native SO(3)lver
rank-one tensors.

The old Fock-space ansatz had 18 columns because it kept three flavour copies
of `Uf`, `Vf`, `Uf0`, `Vf0`, and `mu`.  Their fitted coefficients are equal by
SU(3), so an exact singlet calculation contains only the eight independent
component families returned here.  Same-segment pair projectors retain their
original angular-momentum formula.  Mixed light-heavy projectors are recoupled
from the pairing channel to segment density tensors with the SO(3)lver 9j
convention.
"""
function build_generator_candidates(model::SO3Model)
    model.nm1 >= 4 || throw(ArgumentError(
        "native generator candidates require nm1 >= 4",
    ))
    s = (model.nm1 - 1) / 2
    charge1 = [GetElectronMod(model.nm1, 3, flavor) for flavor in 1:3]
    charge3 = GetElectronMod(model.nm0, 1, 1)

    ff_u = [
        FilterL2(charge1[i] * charge1[j], 2s)
        for i in 1:3 for j in i+1:3
    ]
    ff_v = [
        FilterL2(charge1[i] * charge1[j], 2s - 2)
        for i in 1:3 for j in i+1:3
    ]
    pair00 = charge3 * charge3
    pair00_u = FilterL2(pair00, 6s - 1)
    pair00_v = FilterL2(pair00, 6s - 3)

    zero_shift = zeros(Int64, 4)
    no_shift = zeros(Int64, 4, 2)
    light_to_heavy = zeros(Int64, 4, 2)
    light_to_heavy[1, :] .= (-3, 3)
    heavy_to_light = -light_to_heavy

    density_light = GetDensityMod(
        model.nm1, 3, Matrix{Float64}(I, 3, 3),
    )
    density_heavy = GetDensityMod(model.nm0, 1, ones(1, 1))

    function mixed_pair_decomposition(pair_l::Real)
        pair_l2 = Int64(2pair_l)
        pair_channel = [Int64[pair_l2 pair_l2; pair_l2 2]]
        channels, coefficients = RecoupleAngMom(
            s, 3s, 3s, s, pair_channel, ComplexF64[1],
        )
        return CoupleDecomps(
            [density_light, density_heavy], channels, coefficients, no_shift,
        )
    end

    radius2 = Float64(model.nm1)
    fields = [
        GetElectronObs(model.nm1, 3, flavor; norm_r2=radius2)
        for flavor in 1:3
    ]
    heavy_field = GetElectronObs(model.nm0, 1, 1; norm_r2=radius2)
    triplet = fields[1] * fields[2] * fields[3]

    decompositions = CoupleDecomps[
        SingleSegCouple(
            2, 1,
            reduce(+, FilterL2(mode' * mode, 1) for mode in ff_u),
            2, zero_shift,
        ),
        SingleSegCouple(
            2, 1,
            reduce(+, FilterL2(mode' * mode, 1) for mode in ff_v),
            2, zero_shift,
        ),
        SingleSegCouple(
            2, 2, FilterL2(pair00_u' * pair00_u, 1), 2, zero_shift,
        ),
        SingleSegCouple(
            2, 2, FilterL2(pair00_v' * pair00_v, 1), 2, zero_shift,
        ),
        mixed_pair_decomposition(4s),
        mixed_pair_decomposition(4s - 2),
        -ContactCouple([triplet, heavy_field'], light_to_heavy, 2) +
            ContactCouple([triplet', heavy_field], heavy_to_light, 2),
        SingleSegCouple(
            2, 1, FilterL2(density_light, 1), 2, zero_shift,
        ),
    ]
    real_decompositions = [
        _real_generator_decomposition(decompositions[i], GENERATOR_CANDIDATE_NAMES[i])
        for i in eachindex(decompositions)
    ]
    return SO3GeneratorCandidates(
        GENERATOR_CANDIDATE_NAMES, real_decompositions,
    )
end

"""
Build all eight generator candidates between exact initial and final SO(3)
spaces.  Segment reduced matrix elements are built once for the concatenated
channel list and then shared by the individual composite operators.
"""
function build_generator_operators(
    initial_space::CompSpace{Float64},
    final_space::CompSpace{Float64},
    candidates::SO3GeneratorCandidates;
    disp_std::Bool=true,
    num_th::Int=FuzzifiED.NumThreads,
)
    all_decompositions = CoupleDecomp[]
    ranges = UnitRange{Int}[]
    for decomposition in candidates.decompositions
        first_channel = length(all_decompositions) + 1
        append!(all_decompositions, decomposition)
        push!(ranges, first_channel:length(all_decompositions))
    end
    segment_operators = BuildSegOperators(
        initial_space.sgsp, final_space.sgsp, all_decompositions;
        num_th, disp_std,
    )
    operators = CompOperator{Float64}[]
    for (decomposition, range) in zip(candidates.decompositions, ranges)
        push!(operators, BuildCompOperator(
            initial_space, final_space, decomposition,
            segment_operators[:, range], 2; disp_std,
        ))
    end
    return SO3GeneratorOperators(
        candidates.names, initial_space, final_space, operators,
    )
end

"""Fit the eight native SO(3)lver tensors to `Λ|source> ≈ |target>`."""
function fit_so3_generator(
    source::AbstractVector{<:Real},
    target::AbstractVector{<:Real},
    operator_set::SO3GeneratorOperators;
    svd_rtol::Real=1.0e-10,
)
    length(source) == operator_set.initial_space.dim || throw(DimensionMismatch(
        "source has length $(length(source)); expected $(operator_set.initial_space.dim)",
    ))
    length(target) == operator_set.final_space.dim || throw(DimensionMismatch(
        "target has length $(length(target)); expected $(operator_set.final_space.dim)",
    ))
    source_vector = Vector{Float64}(source)
    target_vector = Vector{Float64}(target)
    design = hcat((operator * source_vector for operator in operator_set.operators)...)
    decomposition = svd(design; full=false)
    cutoff = Float64(svd_rtol) * maximum(decomposition.S; init=0.0)
    inverse_values = [value > cutoff ? inv(value) : 0.0 for value in decomposition.S]
    coefficients = decomposition.V * Diagonal(inverse_values) *
                   decomposition.U' * target_vector
    generated = design * coefficients
    generated_norm2 = real(dot(generated, generated))
    target_norm2 = real(dot(target_vector, target_vector))
    fidelity = generated_norm2 > 0 && target_norm2 > 0 ?
        abs2(dot(target_vector, generated)) / (target_norm2 * generated_norm2) : 0.0
    return (
        names=operator_set.names,
        coefficients=coefficients,
        generated=generated,
        fidelity=fidelity,
        numerical_rank=count(>(cutoff), decomposition.S),
        singular_values=decomposition.S,
    )
end

function apply_so3_generator(
    input::AbstractVector{<:Real},
    operator_set::SO3GeneratorOperators,
    coefficients::AbstractVector{<:Real},
)
    length(input) == operator_set.initial_space.dim || throw(DimensionMismatch(
        "input has length $(length(input)); expected $(operator_set.initial_space.dim)",
    ))
    length(coefficients) == length(operator_set.operators) || throw(DimensionMismatch(
        "received $(length(coefficients)) coefficients for " *
        "$(length(operator_set.operators)) generator candidates",
    ))
    input_vector = Vector{Float64}(input)
    return sum(
        coefficients[i] * (operator_set.operators[i] * input_vector)
        for i in eachindex(coefficients)
    )
end

function _check_generator_hamiltonians(
    operator_set::SO3GeneratorOperators,
    initial_hamiltonian::SO3Hamiltonian,
    final_hamiltonian::SO3Hamiltonian,
)
    operator_set.initial_space === initial_hamiltonian.space || throw(ArgumentError(
        "generator initial space does not match the initial Hamiltonian",
    ))
    operator_set.final_space === final_hamiltonian.space || throw(ArgumentError(
        "generator final space does not match the final Hamiltonian",
    ))
    initial_hamiltonian.workspace.model.representation ==
        final_hamiltonian.workspace.model.representation || throw(ArgumentError(
        "a singlet generator cannot change the internal representation",
    ))
    initial_hamiltonian.couplings == final_hamiltonian.couplings || throw(ArgumentError(
        "initial and final Hamiltonians must use the same couplings",
    ))
    return nothing
end

"""
Apply `Lambda=P+K`, `P`, or `K` without constructing a full-Fock operator.

The commutator with `D=(H-E0)/factor` is evaluated by exact matrix-free
Hamiltonian actions in the initial and final SO(3) blocks.  The vacuum energy
cancels from `[D,Lambda]`, so this works for arbitrary input vectors and does
not require a truncated sum over intermediate energy eigenstates.
"""
function apply_so3_conformal_generator(
    input::AbstractVector{<:Real},
    operator_set::SO3GeneratorOperators,
    coefficients::AbstractVector{<:Real},
    initial_hamiltonian::SO3Hamiltonian,
    final_hamiltonian::SO3Hamiltonian;
    factor::Real,
    generator::Symbol=:lambda,
)
    generator in (:lambda, :p, :k) || throw(ArgumentError(
        "generator must be :lambda, :p, or :k",
    ))
    factor_value = Float64(factor)
    isfinite(factor_value) && factor_value > 0 || throw(ArgumentError(
        "the dilatation scale factor must be positive and finite",
    ))
    _check_generator_hamiltonians(
        operator_set, initial_hamiltonian, final_hamiltonian,
    )
    input_vector = Vector{Float64}(input)
    lambda = apply_so3_generator(input_vector, operator_set, coefficients)
    generator == :lambda && return lambda
    h_input = initial_hamiltonian.operator * input_vector
    commutator = (
        final_hamiltonian.operator * lambda -
        apply_so3_generator(h_input, operator_set, coefficients)
    ) ./ factor_value
    return generator == :p ? (lambda .+ commutator) ./ 2 :
           (lambda .- commutator) ./ 2
end

"""Return `[D,P]-P` or `[D,K]+K` on an arbitrary SO(3)-adapted state."""
function so3_dilatation_residual(
    input::AbstractVector{<:Real},
    operator_set::SO3GeneratorOperators,
    coefficients::AbstractVector{<:Real},
    initial_hamiltonian::SO3Hamiltonian,
    final_hamiltonian::SO3Hamiltonian;
    factor::Real,
    generator::Symbol=:p,
)
    generator in (:p, :k) || throw(ArgumentError(
        "dilatation residual is defined here only for :p or :k",
    ))
    input_vector = Vector{Float64}(input)
    action = apply_so3_conformal_generator(
        input_vector, operator_set, coefficients,
        initial_hamiltonian, final_hamiltonian;
        factor, generator,
    )
    h_input = initial_hamiltonian.operator * input_vector
    action_on_h_input = apply_so3_conformal_generator(
        h_input, operator_set, coefficients,
        initial_hamiltonian, final_hamiltonian;
        factor, generator,
    )
    commutator = (
        final_hamiltonian.operator * action - action_on_h_input
    ) ./ Float64(factor)
    return generator == :p ? commutator .- action : commutator .+ action
end

function _operator_action_matrix(
    states::AbstractMatrix{<:Real},
    operator_set::SO3GeneratorOperators,
    coefficients::AbstractVector{<:Real},
)
    return hcat((
        apply_so3_generator(view(states, :, column), operator_set, coefficients)
        for column in axes(states, 2)
    )...)
end

function _candidate_action_matrix(
    state::AbstractVector{<:Real},
    operator_set::SO3GeneratorOperators,
)
    state_vector = Vector{Float64}(state)
    return hcat((operator * state_vector for operator in operator_set.operators)...)
end

function _apply_hamiltonian_columns(
    hamiltonian::SO3Hamiltonian,
    matrix::AbstractMatrix{<:Real},
)
    return hcat((
        hamiltonian.operator * Vector{Float64}(view(matrix, :, column))
        for column in axes(matrix, 2)
    )...)
end

function _minimum_generalized_eigenpair(
    numerator::AbstractMatrix{<:Real},
    denominator::AbstractMatrix{<:Real};
    rtol::Real=1.0e-10,
)
    denominator_decomposition = eigen(Symmetric(Matrix{Float64}(denominator)))
    maximum_value = maximum(denominator_decomposition.values; init=0.0)
    cutoff = Float64(rtol) * maximum_value
    kept = findall(>(cutoff), denominator_decomposition.values)
    isempty(kept) && error("generator normalization Gram matrix has zero numerical rank")
    whitening = denominator_decomposition.vectors[:, kept] *
                Diagonal(inv.(sqrt.(denominator_decomposition.values[kept])))
    reduced = Symmetric(whitening' * Matrix{Float64}(numerator) * whitening)
    decomposition = eigen(reduced)
    index = argmin(decomposition.values)
    coefficients = whitening * decomposition.vectors[:, index]
    coefficients ./= sqrt(real(dot(coefficients, denominator * coefficients)))
    pivot = argmax(abs.(coefficients))
    coefficients[pivot] < 0 && (coefficients .*= -1)
    return (
        value=max(0.0, Float64(decomposition.values[index])),
        coefficients=coefficients,
        spectrum=Float64.(decomposition.values),
        normalization_rank=length(kept),
        normalization_eigenvalues=Float64.(denominator_decomposition.values),
    )
end

function _fit_primary_annihilating_generator(
    channels,
    vacuum_design::AbstractMatrix{<:Real};
    factor_bounds::Tuple{<:Real,<:Real},
    fit_primary_labels,
    vacuum_weight::Real=1.0,
    dilatation_weight::Real=1.0,
    descendant_weight::Real=1.0,
    gram_rtol::Real=1.0e-10,
)
    lower, upper = Float64.(factor_bounds)
    0 < lower < upper || throw(ArgumentError(
        "factor_bounds must be positive and increasing",
    ))
    labels = Set(Symbol.(collect(fit_primary_labels)))
    available = Set(getproperty.(channels, :label))
    isempty(labels) && throw(ArgumentError("fit_primary_labels cannot be empty"))
    issubset(labels, available) || throw(ArgumentError(
        "unknown fit primary labels: $(collect(setdiff(labels, available)))",
    ))
    selected = filter(channel -> channel.label in labels, channels)
    isfinite(dilatation_weight) && dilatation_weight >= 0 || throw(ArgumentError(
        "dilatation_weight must be finite and non-negative",
    ))
    isfinite(descendant_weight) && descendant_weight >= 0 || throw(ArgumentError(
        "descendant_weight must be finite and non-negative",
    ))
    candidate_count = size(vacuum_design, 2)
    denominator = zeros(Float64, candidate_count, candidate_count)
    for channel in selected
        denominator .+= channel.weight .* (channel.lambda_design' * channel.lambda_design)
    end
    vacuum_gram = Float64(vacuum_weight) .* (vacuum_design' * vacuum_design)

    function fit_at_factor(factor)
        numerator = copy(vacuum_gram)
        for channel in selected
            k_design = (
                channel.lambda_design .- channel.commutator_numerator ./ factor
            ) ./ 2
            numerator .+= channel.weight .* (k_design' * k_design)
            if dilatation_weight > 0
                dilatation_design = (
                    channel.second_commutator_numerator ./ factor^2 .-
                    channel.lambda_design
                ) ./ 2
                numerator .+= Float64(dilatation_weight) .* channel.weight .*
                    (dilatation_design' * dilatation_design)
            end
            if descendant_weight > 0
                p_leakage_design = (
                    channel.lambda_leakage_design .+
                    channel.commutator_leakage_numerator ./ factor
                ) ./ 2
                numerator .+= Float64(descendant_weight) .* channel.weight .*
                    (p_leakage_design' * p_leakage_design)
            end
        end
        return _minimum_generalized_eigenpair(
            numerator, denominator; rtol=gram_rtol,
        )
    end

    log_scan = collect(range(log(lower), log(upper); length=41))
    scan_factors = exp.(log_scan)
    scan_values = [fit_at_factor(factor).value for factor in scan_factors]
    scan_best = argmin(scan_values)
    factors = Float64[lower, scan_factors[scan_best], upper]
    optimizer_converged = true
    if 1 < scan_best < length(log_scan)
        optimization = Optim.optimize(
            log_factor -> fit_at_factor(exp(log_factor)).value,
            log_scan[scan_best - 1], log_scan[scan_best + 1], Optim.Brent();
            rel_tol=1.0e-10, abs_tol=1.0e-12,
        )
        push!(factors, exp(Optim.minimizer(optimization)))
        optimizer_converged = Optim.converged(optimization)
    end
    fits = map(fit_at_factor, factors)
    best_index = argmin(getproperty.(fits, :value))
    factor = factors[best_index]
    fit = fits[best_index]
    boundary_distance = min(log(factor / lower), log(upper / factor))
    return merge(fit, (
        factor=factor,
        factor_bounds=(lower, upper),
        factor_at_boundary=boundary_distance < 1.0e-5,
        factor_scan=scan_factors,
        factor_scan_values=scan_values,
        fit_primary_labels=sort!(collect(labels)),
        vacuum_weight=Float64(vacuum_weight),
        dilatation_weight=Float64(dilatation_weight),
        descendant_weight=Float64(descendant_weight),
        optimizer_converged=optimizer_converged,
    ))
end

_vector_target_ells(ell::Int) = ell == 0 ? [1] : collect((ell - 1):(ell + 1))

function _default_conformal_block_counts(primary_specs)
    counts = Dict{Tuple{Symbol,Int},Int}((:singlet, 0) => 4)
    for spec in primary_specs
        key = (spec.representation, spec.ell)
        counts[key] = max(get(counts, key, 0), spec.rank + 2)
    end
    return counts
end

"""
Audit a common microscopic `Lambda`, `P`, and `K` against several primaries.

The generator coefficients and the cylinder energy scale are determined by a
basis-invariant generalized eigenproblem that minimizes the full-block vacuum
and `K|primary>` norms.  All target-state norms use exact SO(3)lver operator
actions; no full-Fock construction and no low-energy intermediate-state
truncation is used.  In addition to fixed candidate primaries, the function
diagonalizes `K^dagger K` inside each requested low-energy source subspace.
"""
function analyze_so3_conformal_algebra(
    singlet_workspace::SO3Workspace,
    adjoint_workspace::SO3Workspace,
    couplings::Couplings;
    primary_specs=DEFAULT_CONFORMAL_PRIMARY_SPECS,
    block_counts=_default_conformal_block_counts(primary_specs),
    fit_primary_labels=getproperty.(primary_specs, :label),
    factor_bounds::Tuple{<:Real,<:Real}=(0.01, 0.10),
    vacuum_weight::Real=1.0,
    dilatation_weight::Real=1.0,
    descendant_weight::Real=1.0,
    descendant_state_count::Int=6,
    eig_tol::Real=1.0e-8,
    ncv::Int=18,
    gram_rtol::Real=1.0e-10,
    disp_std::Bool=true,
)
    singlet_workspace.model.representation == :singlet || throw(ArgumentError(
        "the first workspace must be the singlet representation",
    ))
    adjoint_workspace.model.representation == :adjoint || throw(ArgumentError(
        "the second workspace must be the adjoint representation",
    ))
    singlet_workspace.model.nm1 == adjoint_workspace.model.nm1 || throw(ArgumentError(
        "singlet and adjoint workspaces must use the same system size",
    ))
    singlet_workspace.heavy_space_mode == adjoint_workspace.heavy_space_mode ||
        throw(ArgumentError("singlet and adjoint workspaces must use the same projection"))
    isempty(primary_specs) && throw(ArgumentError("primary_specs cannot be empty"))
    descendant_state_count >= 1 || throw(ArgumentError(
        "descendant_state_count must be positive",
    ))
    MottJainED.validate(couplings)

    workspaces = Dict(
        :singlet => singlet_workspace,
        :adjoint => adjoint_workspace,
    )
    required_blocks = Set{Tuple{Symbol,Int}}(((:singlet, 0),))
    for spec in primary_specs
        representation_data(spec.representation)
        spec.ell >= 0 || throw(ArgumentError("primary angular momentum must be non-negative"))
        spec.rank >= 1 || throw(ArgumentError("primary rank must be positive"))
        push!(required_blocks, (spec.representation, spec.ell))
        for target_ell in _vector_target_ells(spec.ell)
            push!(required_blocks, (spec.representation, target_ell))
        end
    end
    push!(required_blocks, (:singlet, 1))

    hamiltonians = Dict{Tuple{Symbol,Int},SO3Hamiltonian}()
    for (representation, ell) in sort!(collect(required_blocks))
        hamiltonians[(representation, ell)] = build_hamiltonian(
            workspaces[representation], ell, couplings; disp_std,
        )
    end

    requested_counts = Dict{Tuple{Symbol,Int},Int}(
        (Symbol(key[1]), Int(key[2])) => Int(value)
        for (key, value) in pairs(block_counts)
    )
    k2_source_keys = Set(keys(requested_counts))
    requested_counts[(:singlet, 0)] = max(get(requested_counts, (:singlet, 0), 0), 2)
    for spec in primary_specs
        key = (spec.representation, spec.ell)
        requested_counts[key] = max(get(requested_counts, key, 0), spec.rank)
    end
    for key in required_blocks
        requested_counts[key] = max(
            get(requested_counts, key, 0), descendant_state_count,
        )
    end
    energies = Dict{Tuple{Symbol,Int},Vector{Float64}}()
    states = Dict{Tuple{Symbol,Int},Matrix{Float64}}()
    for (key, count) in requested_counts
        haskey(hamiltonians, key) || continue
        count >= 1 || throw(ArgumentError("block_counts must be positive"))
        block_energies, block_states = solve(
            hamiltonians[key]; k=count, tol=Float64(eig_tol), ncv,
            vectors=true, disp_std,
        )
        _canonicalize_state_signs!(block_states)
        energies[key] = block_energies
        states[key] = block_states
    end

    candidates = Dict(
        representation => build_generator_candidates(workspace.model)
        for (representation, workspace) in workspaces
    )
    operator_cache = Dict{Tuple{Symbol,Int,Int},SO3GeneratorOperators}()
    function operators(representation::Symbol, initial_ell::Int, final_ell::Int)
        key = (representation, initial_ell, final_ell)
        return get!(operator_cache, key) do
            build_generator_operators(
                hamiltonians[(representation, initial_ell)].space,
                hamiltonians[(representation, final_ell)].space,
                candidates[representation]; disp_std,
            )
        end
    end

    vacuum = view(states[(:singlet, 0)], :, 1)
    vacuum_design = _candidate_action_matrix(
        vacuum, operators(:singlet, 0, 1),
    )
    channels = NamedTuple[]
    for spec in primary_specs
        source_key = (spec.representation, spec.ell)
        source_energy = energies[source_key][spec.rank]
        source = view(states[source_key], :, spec.rank)
        weight = inv(2spec.ell + 1)
        for target_ell in _vector_target_ells(spec.ell)
            final_hamiltonian = hamiltonians[(spec.representation, target_ell)]
            operator_set = operators(spec.representation, spec.ell, target_ell)
            lambda_design = _candidate_action_matrix(source, operator_set)
            commutator_numerator = _apply_hamiltonian_columns(
                final_hamiltonian, lambda_design,
            ) .- source_energy .* lambda_design
            second_commutator_numerator = _apply_hamiltonian_columns(
                final_hamiltonian, commutator_numerator,
            ) .- source_energy .* commutator_numerator
            target_states = states[(spec.representation, target_ell)]
            lambda_leakage_design = lambda_design .-
                target_states * (target_states' * lambda_design)
            commutator_leakage_numerator = commutator_numerator .-
                target_states * (target_states' * commutator_numerator)
            push!(channels, (
                label=spec.label,
                representation=spec.representation,
                source_ell=spec.ell,
                target_ell=target_ell,
                source_rank=spec.rank,
                source_energy=source_energy,
                weight=weight,
                lambda_design=lambda_design,
                commutator_numerator=commutator_numerator,
                second_commutator_numerator=second_commutator_numerator,
                lambda_leakage_design=lambda_leakage_design,
                commutator_leakage_numerator=commutator_leakage_numerator,
                target_states=target_states,
                operator_set=operator_set,
                initial_hamiltonian=hamiltonians[source_key],
                final_hamiltonian=final_hamiltonian,
            ))
        end
    end

    fit = _fit_primary_annihilating_generator(
        channels, vacuum_design;
        factor_bounds, fit_primary_labels, vacuum_weight,
        dilatation_weight, descendant_weight, gram_rtol,
    )
    factor = fit.factor
    direction = fit.coefficients

    ground_energy = energies[(:singlet, 0)][1]
    scalar_commutator_coefficients = Float64[]
    scalar_commutator_targets = Float64[]
    for spec in primary_specs
        spec.ell == 0 && spec.label in fit.fit_primary_labels || continue
        channel = only(filter(
            row -> row.label == spec.label && row.target_ell == 1,
            channels,
        ))
        lambda = channel.lambda_design * direction
        raw_commutator = channel.commutator_numerator * direction
        # For L=0,m=0 and the q=0 Cartesian/spherical component, the square
        # of the Wigner-Eckart 3j factor is 1/3.  Multiplying
        # <[Kz,Pz]>=2D by `factor` removes the cylinder scale from both sides.
        push!(scalar_commutator_coefficients,
              real(dot(lambda, raw_commutator)) / 3)
        push!(scalar_commutator_targets,
              2 * (channel.source_energy - ground_energy))
    end
    normalization_denominator = sum(abs2, scalar_commutator_coefficients)
    normalization_scale2 = normalization_denominator > eps(Float64) ?
        dot(scalar_commutator_coefficients, scalar_commutator_targets) /
        normalization_denominator : NaN
    commutator_normalization_valid = isfinite(normalization_scale2) &&
                                     normalization_scale2 > 0
    normalization_scale = commutator_normalization_valid ?
        sqrt(normalization_scale2) : 1.0
    coefficients = normalization_scale .* direction
    scalar_commutator_residuals = commutator_normalization_valid ?
        normalization_scale2 .* scalar_commutator_coefficients .-
            scalar_commutator_targets :
        fill(NaN, length(scalar_commutator_targets))
    fit = merge(fit, (
        coefficients=coefficients,
        direction_coefficients=direction,
        commutator_normalization_valid=commutator_normalization_valid,
        commutator_normalization_scale=normalization_scale,
        scalar_commutator_residuals=scalar_commutator_residuals,
    ))

    channel_rows = NamedTuple[]
    primary_rows = NamedTuple[]
    for spec in primary_specs
        selected = filter(channel -> channel.label == spec.label, channels)
        totals = Dict(
            :lambda => 0.0, :p => 0.0, :k => 0.0,
            :p_residual => 0.0, :k_residual => 0.0,
            :lambda_delta => 0.0,
        )
        for channel in selected
            lambda = channel.lambda_design * coefficients
            d_lambda = channel.commutator_numerator * coefficients ./ factor
            p_action = (lambda .+ d_lambda) ./ 2
            k_action = (lambda .- d_lambda) ./ 2
            delta_p = (
                channel.final_hamiltonian.operator * p_action .-
                channel.source_energy .* p_action
            ) ./ factor
            delta_k = (
                channel.final_hamiltonian.operator * k_action .-
                channel.source_energy .* k_action
            ) ./ factor
            p_residual = delta_p .- p_action
            k_residual = delta_k .+ k_action
            weight = channel.weight
            lambda_norm2 = weight * real(dot(lambda, lambda))
            p_norm2 = weight * real(dot(p_action, p_action))
            k_norm2 = weight * real(dot(k_action, k_action))
            p_residual_norm2 = weight * real(dot(p_residual, p_residual))
            k_residual_norm2 = weight * real(dot(k_residual, k_residual))
            projected_p = channel.target_states * (channel.target_states' * p_action)
            p_leakage_norm2 = weight * real(dot(
                p_action - projected_p, p_action - projected_p,
            ))
            lambda_delta = weight * real(dot(lambda, d_lambda))
            totals[:lambda] += lambda_norm2
            totals[:p] += p_norm2
            totals[:k] += k_norm2
            totals[:p_residual] += p_residual_norm2
            totals[:k_residual] += k_residual_norm2
            totals[:p_leakage] = get(totals, :p_leakage, 0.0) + p_leakage_norm2
            totals[:lambda_delta] += lambda_delta
            push!(channel_rows, (
                label=spec.label,
                representation=spec.representation,
                source_ell=spec.ell,
                source_rank=spec.rank,
                target_ell=channel.target_ell,
                source_energy=channel.source_energy,
                lambda_norm2=lambda_norm2,
                p_norm2=p_norm2,
                k_norm2=k_norm2,
                k_fraction=k_norm2 / max(lambda_norm2, eps(Float64)),
                p_dilatation_residual_norm2=p_residual_norm2,
                k_dilatation_residual_norm2=k_residual_norm2,
                p_low_energy_leakage_norm2=p_leakage_norm2,
            ))
        end
        scalar_commutator_lhs = spec.ell == 0 ?
            (totals[:p] - totals[:k]) / 3 : NaN
        scalar_commutator_target = spec.ell == 0 ?
            2 * (energies[(spec.representation, spec.ell)][spec.rank] - ground_energy) /
                factor : NaN
        push!(primary_rows, (
            label=spec.label,
            representation=spec.representation,
            ell=spec.ell,
            rank=spec.rank,
            energy=energies[(spec.representation, spec.ell)][spec.rank],
            lambda_norm2=totals[:lambda],
            p_norm2=totals[:p],
            k_norm2=totals[:k],
            k_fraction=totals[:k] / max(totals[:lambda], eps(Float64)),
            p_dilatation_fraction=totals[:p_residual] /
                max(totals[:p], eps(Float64)),
            k_dilatation_fraction=totals[:k_residual] /
                max(totals[:k], eps(Float64)),
            p_low_energy_leakage_fraction=totals[:p_leakage] /
                max(totals[:p], eps(Float64)),
            lambda_mean_scaled_gap=totals[:lambda_delta] /
                max(totals[:lambda], eps(Float64)),
            kp_commutator_lhs=scalar_commutator_lhs,
            kp_commutator_target=scalar_commutator_target,
            kp_commutator_fractional_residual=spec.ell == 0 ?
                (scalar_commutator_lhs - scalar_commutator_target) /
                    max(abs(scalar_commutator_target), eps(Float64)) : NaN,
            used_for_fit=spec.label in fit.fit_primary_labels,
        ))
    end

    vacuum_norm2 = real(dot(vacuum_design * coefficients, vacuum_design * coefficients))
    k2_results = Dict{Tuple{Symbol,Int},NamedTuple}()
    for (source_key, source_states) in states
        source_key in k2_source_keys || continue
        representation, source_ell = source_key
        lambda_gram = zeros(Float64, size(source_states, 2), size(source_states, 2))
        k_gram = zeros(Float64, size(source_states, 2), size(source_states, 2))
        for target_ell in _vector_target_ells(source_ell)
            haskey(hamiltonians, (representation, target_ell)) || continue
            operator_set = operators(representation, source_ell, target_ell)
            lambda_map = _operator_action_matrix(
                source_states, operator_set, coefficients,
            )
            k_map = hcat((
                apply_so3_conformal_generator(
                    view(source_states, :, column), operator_set, coefficients,
                    hamiltonians[source_key], hamiltonians[(representation, target_ell)];
                    factor, generator=:k,
                )
                for column in axes(source_states, 2)
            )...)
            all(isfinite, lambda_map) || error(
                "non-finite Lambda action in representation=$representation " *
                "L=$source_ell->$target_ell",
            )
            all(isfinite, k_map) || error(
                "non-finite K action in representation=$representation " *
                "L=$source_ell->$target_ell",
            )
            multiplet_weight = inv(2source_ell + 1)
            lambda_gram .+= multiplet_weight .* (lambda_map' * lambda_map)
            k_gram .+= multiplet_weight .* (k_map' * k_map)
        end
        decomposition = eigen(Symmetric(k_gram))
        mode_lambda_norm2 = [
            real(dot(view(decomposition.vectors, :, mode),
                     lambda_gram * view(decomposition.vectors, :, mode)))
            for mode in axes(decomposition.vectors, 2)
        ]
        mode_energy = [
            sum(abs2.(view(decomposition.vectors, :, mode)) .* energies[source_key])
            for mode in axes(decomposition.vectors, 2)
        ]
        k2_results[source_key] = (
            eigenvalues=Float64.(decomposition.values),
            eigenvectors=Matrix{Float64}(decomposition.vectors),
            lambda_norm2=mode_lambda_norm2,
            k_fraction=Float64.(decomposition.values) ./
                max.(mode_lambda_norm2, eps(Float64)),
            energy_expectation=mode_energy,
            source_energies=energies[source_key],
        )
    end

    return (
        fit=fit,
        primary_rows=primary_rows,
        channel_rows=channel_rows,
        k2=k2_results,
        vacuum_norm2=vacuum_norm2,
        energies=energies,
        dimensions=Dict(key => hamiltonian.space.dim for (key, hamiltonian) in hamiltonians),
        heavy_space_mode=singlet_workspace.heavy_space_mode,
        generator_candidate_names=GENERATOR_CANDIDATE_NAMES,
    )
end

"""Apply a fitted generator and resolve its normalized weight over target states."""
function so3_generator_overlaps(
    input::AbstractVector{<:Real},
    targets::AbstractMatrix{<:Real},
    operator_set::SO3GeneratorOperators,
    coefficients::AbstractVector{<:Real},
)
    size(targets, 1) == operator_set.final_space.dim || throw(DimensionMismatch(
        "target vectors have length $(size(targets, 1)); " *
        "expected $(operator_set.final_space.dim)",
    ))
    generated = apply_so3_generator(input, operator_set, coefficients)
    norm2 = real(dot(generated, generated))
    values = [
        norm2 > 0 ? abs2(dot(view(targets, :, i), generated)) / norm2 : 0.0
        for i in axes(targets, 2)
    ]
    return (values=values, total=sum(values), norm2=norm2, generated=generated)
end

function _canonicalize_state_signs!(states::Matrix{Float64})
    for column in axes(states, 2)
        vector = view(states, :, column)
        pivot = argmax(abs.(vector))
        vector[pivot] < 0 && (vector .*= -1)
    end
    return states
end

"""
Solve the scalar `L=0,1,2` blocks and run the complete native generator audit.
The returned rank indices are already physical SO(3) multiplets; no `L²`
projection or raw-Cartan duplicate merging is involved.
"""
function analyze_scalar_generator(
    workspace::SO3Workspace,
    couplings::Couplings;
    l0_count::Int=6,
    l2_count::Int=6,
    eig_tol::Float64=1.0e-8,
    ncv::Int=max(18, 2 * max(l0_count, l2_count)),
    svd_rtol::Float64=1.0e-10,
    disp_std::Bool=true,
)
    l0_count >= 3 || throw(ArgumentError("l0_count must be at least 3"))
    l2_count >= 2 || throw(ArgumentError("l2_count must be at least 2"))
    workspace.model.representation == :singlet || throw(ArgumentError(
        "scalar generator analysis requires a singlet SO3Workspace",
    ))
    h0 = build_hamiltonian(workspace, 0, couplings; disp_std)
    h1 = build_hamiltonian(workspace, 1, couplings; disp_std)
    h2 = build_hamiltonian(workspace, 2, couplings; disp_std)
    energies0, states0 = solve(
        h0; k=l0_count, tol=eig_tol, ncv, vectors=true, disp_std,
    )
    energies1, states1 = solve(
        h1; k=1, tol=eig_tol, ncv, vectors=true, disp_std,
    )
    energies2, states2 = solve(
        h2; k=l2_count, tol=eig_tol, ncv, vectors=true, disp_std,
    )
    _canonicalize_state_signs!(states0)
    _canonicalize_state_signs!(states1)
    _canonicalize_state_signs!(states2)

    candidates = build_generator_candidates(workspace.model)
    operators01 = build_generator_operators(
        h0.space, h1.space, candidates; disp_std,
    )
    fit = fit_so3_generator(
        view(states0, :, 2), view(states1, :, 1), operators01; svd_rtol,
    )
    operators10 = build_generator_operators(
        h1.space, h0.space, candidates; disp_std,
    )
    overlap0 = so3_generator_overlaps(
        view(states1, :, 1), states0, operators10, fit.coefficients,
    )
    operators12 = build_generator_operators(
        h1.space, h2.space, candidates; disp_std,
    )
    overlap2 = so3_generator_overlaps(
        view(states1, :, 1), states2, operators12, fit.coefficients,
    )
    return (
        fit=fit,
        l0_overlap=overlap0,
        l2_overlap=overlap2,
        energies=(l0=energies0, l1=energies1, l2=energies2),
        dimensions=(l0=h0.space.dim, l1=h1.space.dim, l2=h2.space.dim),
    )
end

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

function _relation_energy(blocks, state::Tuple{Symbol,Int,Int})
    representation, ell, rank = state
    return _block_energies(blocks, representation, ell, rank)[rank]
end

function _normalize_score_terms(terms)
    selected = Symbol.(collect(terms))
    isempty(selected) && throw(ArgumentError("CFT score needs at least one relation"))
    length(unique(selected)) == length(selected) || throw(ArgumentError(
        "CFT score terms must not contain duplicates",
    ))
    unknown = filter(term -> !haskey(CFT_RELATION_SPECS, term), selected)
    isempty(unknown) || throw(ArgumentError(
        "unknown SO(3)lver CFT score terms: $(join(String.(unknown), ", "))",
    ))
    return selected
end

"""
Evaluate selected CFT energy relations from exact physical blocks.

The default remains the project's fixed five-relation locator.  Provisional
higher descendants can be requested explicitly only by callers that separately
audit their state identity.  The old Fock-basis workflow selected raw ranks
before merging the two Cartan-zero copies of every SU(3) octet.  An adjoint
highest-weight block contains each octet once, so old raw rank 3 for curl-J
becomes physical rank 2 here.
"""
function score_cft_blocks(blocks; terms=CFT_SCORE_TERMS)
    selected_terms = _normalize_score_terms(terms)
    singlet_l0 = _block_energies(blocks, :singlet, 0, 2)
    adjoint_l0 = _block_energies(blocks, :adjoint, 0, 1)

    block_minima = [
        (representation=rep, ell=ell,
         energy=first(_block_energies(blocks, rep, ell, 1)))
        for (rep, ell) in CFT_BLOCK_KEYS
    ]
    ground_block = block_minima[argmin(getproperty.(block_minima, :energy))]
    ground = ground_block.energy

    raw_gaps = Float64[]
    targets = Float64[]
    labels = String[]
    for term in selected_terms
        spec = CFT_RELATION_SPECS[term]
        upper = _relation_energy(blocks, spec.upper)
        lower = isnothing(spec.lower) ? ground : _relation_energy(blocks, spec.lower)
        push!(raw_gaps, upper - lower)
        push!(targets, spec.target)
        push!(labels, spec.label)
    end
    factor = dot(raw_gaps, targets) / dot(targets, targets)
    isfinite(factor) && factor > 0 || error("fitted CFT scale factor is not positive")
    scaled_gaps = raw_gaps ./ factor
    residuals = scaled_gaps .- targets
    q = sqrt(sum(abs2, residuals) / length(residuals))
    return (
        q=q,
        factor=factor,
        delta_s=(singlet_l0[2] - ground) / factor,
        delta_o=(adjoint_l0[1] - ground) / factor,
        raw_gaps=raw_gaps,
        target_gaps=targets,
        scaled_gaps=scaled_gaps,
        residuals=residuals,
        labels=labels,
        terms=selected_terms,
        ground_energy=ground,
        ground_representation=ground_block.representation,
        ground_ell=ground_block.ell,
        ground_is_singlet_l0=(ground_block.representation == :singlet && ground_block.ell == 0),
    )
end

export SO3Model, SO3Workspace, SO3Hamiltonian,
       build_so3_model, build_workspace, build_hamiltonian,
       build_laughlin13_heavy_space, build_heavy_v1_operator,
       heavy_v1_terms, laughlin13_root_counts, HEAVY_SPACE_MODES,
       retune!, solve, sector_dimension, representation_data,
       SO3GeneratorCandidates, SO3GeneratorOperators,
       build_generator_candidates, build_generator_operators,
       fit_so3_generator, apply_so3_generator, so3_generator_overlaps,
       apply_so3_conformal_generator, so3_dilatation_residual,
       analyze_scalar_generator, analyze_so3_conformal_algebra,
       GENERATOR_CANDIDATE_NAMES, DEFAULT_CONFORMAL_PRIMARY_SPECS,
       score_cft_blocks, CFT_BLOCK_KEYS, CFT_SCORE_TERMS,
       CFT_STABLE_SIX_TERMS, CFT_AUDITED_SEVEN_TERMS, CFT_RELATION_SPECS

end
