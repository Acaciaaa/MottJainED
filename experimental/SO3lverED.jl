module SO3lverED

using FuzzifiED
using FuzzifiED.JackToolkit
using FuzzifiED.SO3lver
using LinearAlgebra
using MottJainED

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
       analyze_scalar_generator, GENERATOR_CANDIDATE_NAMES,
       score_cft_blocks, CFT_BLOCK_KEYS, CFT_SCORE_TERMS,
       CFT_STABLE_SIX_TERMS, CFT_AUDITED_SEVEN_TERMS, CFT_RELATION_SPECS

end
