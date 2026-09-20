using Test
using MottJainED
using LinearAlgebra

include(joinpath(@__DIR__, "..", "experimental", "SO3lverED.jl"))
using .SO3lverED
include(joinpath(@__DIR__, "..", "experimental", "SO3ParameterSearch.jl"))
using .SO3ParameterSearch

function conventional_catalog(couplings)
    model = build_model(nm1=2)
    settings = SolverSettings(k=100, dense_cutoff=10_000, warm_start=false)
    cache = prepare_spectrum(model, couplings, settings)
    states = solve_spectrum(cache, couplings.mu)
    catalog, rejected = level_catalog(states)
    @test isempty(rejected)
    return catalog
end

function reference_energies(catalog, ell, c2)
    levels = get(catalog, (ell * (ell + 1), c2), nothing)
    return isnothing(levels) ? Float64[] : Float64[level.energy for level in levels]
end

function isolated_component(component)
    return Couplings(; (
        name => (name == component ? 1.0 : 0.0)
        for name in fieldnames(Couplings)
    )...)
end

@testset "SO(3)lver exact small-size spectra" begin
    # Non-round coefficients make this a simultaneous regression test of all
    # density, Laplacian, conversion, and chemical-potential components.
    couplings = Couplings(
        Uf=0.37, Uf0=-0.21, U0=0.83,
        Vf=0.19, Vf0=-0.17, V0=0.11,
        t=0.29, mu=-0.07,
    )
    catalog = conventional_catalog(couplings)

    singlet_model = build_so3_model(nm1=2, representation=:singlet)
    @test singlet_model.weight == (0, 0)
    @test singlet_model.c2 == 0.0
    # The exact-zero l=0 Laplacian contact channel must be absent.
    @test length(singlet_model.all_decompositions) == 10
    singlet_workspace = build_workspace(singlet_model; disp_std=false)
    singlet_hamiltonians = [
        build_hamiltonian(singlet_workspace, ell, couplings; disp_std=false)
        for ell in 0:2
    ]
    # Inspect the actual full CompOperator matrix.  This must not be replaced
    # by Symmetric(Matrix(...)): doing so would hide a wrong relative sign
    # between the two fermion-odd conversion channels.
    for hamiltonian in singlet_hamiltonians
        for component in fieldnames(Couplings)
            retune!(hamiltonian, isolated_component(component))
            matrix = Matrix(hamiltonian.operator; disp_std=false)
            @test matrix ≈ matrix' atol=2e-12 rtol=2e-12
        end
        retune!(hamiltonian, couplings)
    end
    for (ell, hamiltonian) in enumerate(singlet_hamiltonians)
        physical_ell = ell - 1
        reference = reference_energies(catalog, physical_ell, 0)
        @test sector_dimension(hamiltonian) == length(reference)
        @test solve(hamiltonian; k=100, disp_std=false) ≈ reference atol=2e-11
    end
    # Force the Krylov/matrix-free branch in a three-dimensional block.  The
    # first two values must match conventional ED without dense symmetrization.
    matrix_free = solve(
        singlet_hamiltonians[1]; k=2, dense_cutoff=0, disp_std=false,
    )
    @test matrix_free ≈ reference_energies(catalog, 0, 0)[1:2] atol=2e-11
    warm_matrix_free = solve(
        singlet_hamiltonians[1]; k=2, dense_cutoff=0,
        initvec=ones(sector_dimension(singlet_hamiltonians[1])), disp_std=false,
    )
    @test warm_matrix_free ≈ matrix_free atol=2e-11
    @test_throws DimensionMismatch solve(
        singlet_hamiltonians[1]; k=2, dense_cutoff=0,
        initvec=ones(2), disp_std=false,
    )

    adjoint_model = build_so3_model(nm1=2, representation=:adjoint)
    @test adjoint_model.weight == (1, 3)
    @test adjoint_model.c2 == 3.0
    adjoint_workspace = build_workspace(
        adjoint_model; heavy_space=singlet_workspace.heavy_space, disp_std=false,
    )
    @test adjoint_workspace.heavy_space === singlet_workspace.heavy_space
    for ell in 0:2
        hamiltonian = build_hamiltonian(
            adjoint_workspace, ell, couplings; disp_std=false,
        )
        reference = reference_energies(catalog, ell, 3)
        @test sector_dimension(hamiltonian) == length(reference)
        @test solve(hamiltonian; k=100, disp_std=false) ≈ reference atol=2e-11
    end

    # Retuning changes only the channel coefficients; the expensive segment
    # states and reduced matrix elements remain valid.
    retuned = Couplings(
        Uf=0.51, Uf0=1.27, U0=3.9,
        Vf=-0.03, Vf0=0.44, V0=0.62,
        t=-0.23, mu=0.08,
    )
    retuned_catalog = conventional_catalog(retuned)
    hamiltonian = singlet_hamiltonians[1]
    original_space = hamiltonian.space
    retune!(hamiltonian, retuned)
    @test hamiltonian.space === original_space
    @test solve(hamiltonian; k=100, disp_std=false) ≈
          reference_energies(retuned_catalog, 0, 0) atol=2e-11
end

@testset "SO(3)lver physical-block CFT score" begin
    blocks = Dict{Tuple{Symbol,Int},Vector{Float64}}(
        (:singlet, 0) => [0.0, 1.2, 3.2],
        (:singlet, 1) => [2.2],
        (:singlet, 2) => [3.0, 3.2],
        (:adjoint, 0) => [1.6, 3.6],
        (:adjoint, 1) => [2.0, 3.0, 4.0],
        (:adjoint, 2) => [3.0],
    )
    score = score_cft_blocks(blocks)
    @test score.q < 1e-14
    @test score.factor ≈ 1.0
    @test score.delta_s ≈ 1.2
    @test score.delta_o ≈ 1.6
    @test score.ground_is_singlet_l0
    @test score.raw_gaps ≈ [1, 2, 3, 3, 3]
    @test score.labels == ["dS-S", "J", "curlJ", "dJ(rank1)", "T(rank1)"]
    @test score.terms == collect(CFT_SCORE_TERMS)
    seven = score_cft_blocks(blocks; terms=CFT_AUDITED_SEVEN_TERMS)
    @test seven.q < 1e-14
    @test seven.raw_gaps ≈ [1, 1, 2, 2, 3, 3, 3]
    @test seven.terms == collect(CFT_AUDITED_SEVEN_TERMS)
    holdout = score_cft_blocks(blocks; terms=[:boxo_o, :boxj_j])
    @test holdout.q < 1e-14
    @test holdout.raw_gaps ≈ [2, 2]
    @test_throws ArgumentError score_cft_blocks(blocks; terms=Symbol[])
    @test_throws ArgumentError score_cft_blocks(blocks; terms=[:j, :j])
    @test_throws ArgumentError score_cft_blocks(blocks; terms=[:not_a_relation])
    @test_throws ArgumentError score_cft_blocks(delete!(copy(blocks), (:adjoint, 0)))
end

@testset "SO(3)lver parameter-search identity guards" begin
    center = [1.834, 0.41, 0.525, 0.1218143040052]
    scales = [0.05, 0.025, 0.01, 0.0015]
    search_point = [-1.0, 0.5, 2.0, -0.25]
    physical_point = parameter_search_to_physical(search_point, center, scales)
    @test physical_point ≈ [1.784, 0.4225, 0.545, 0.1214393040052]
    @test physical_to_parameter_search(physical_point, center, scales) ≈ search_point
    @test_throws DimensionMismatch parameter_search_to_physical([1.0], center, scales)
    @test_throws ArgumentError physical_to_parameter_search(center, center, -scales)

    specifications = (
        (label=:boxS, representation=:singlet, ell=0, rank=3),
        (label=:ddS, representation=:singlet, ell=2, rank=2),
    )
    reference = Dict{Tuple{Symbol,Int},Matrix{Float64}}(
        (:singlet, 0) => Matrix{Float64}(I, 4, 4),
        (:singlet, 2) => Matrix{Float64}(I, 4, 4),
    )
    unchanged = Dict(key => copy(value) for (key, value) in reference)
    stable = track_reference_states(
        reference, unchanged; specifications, minimum_overlap=0.9,
    )
    @test stable.passed
    @test all(row -> row.same_rank && row.expected_overlap ≈ 1, stable.rows)

    swapped = Dict(key => copy(value) for (key, value) in reference)
    swapped[(:singlet, 0)][:, [3, 4]] = swapped[(:singlet, 0)][:, [4, 3]]
    unstable = track_reference_states(
        reference, swapped; specifications, minimum_overlap=0.9,
    )
    @test !unstable.passed
    @test only(filter(row -> row.label == :boxS, unstable.rows)).best_rank == 4

    terms = [:a, :b, :c, :d, :e]
    plus = Dict{Symbol,Any}()
    minus = Dict{Symbol,Any}()
    parameters = [:p1, :p2, :p3, :p4]
    for (index, parameter) in enumerate(parameters)
        direction = zeros(5)
        direction[index] = 1
        plus[parameter] = (terms=terms, residuals=direction)
        minus[parameter] = (terms=terms, residuals=-direction)
    end
    sensitivity = normalized_residual_jacobian(plus, minus, parameters)
    @test sensitivity.numerical_rank == 4
    @test sensitivity.condition_number ≈ 1
end

@testset "SO(3)lver score matches conventional raw-rank score" begin
    couplings = Couplings(
        Uf=0.46, Uf0=1.834, U0=4.14,
        Vf=0.0, Vf0=0.41, V0=0.525,
        t=0.5, mu=0.1217,
    )
    model = build_model(nm1=3)
    settings = SolverSettings(k=100, dense_cutoff=10_000, warm_start=false)
    conventional = solve_spectrum(prepare_spectrum(model, couplings, settings), couplings.mu)
    legacy = cft_score(
        conventional;
        terms=[:ds_s, :j, :curlj, :dj_rank1, :t_rank1], metric=:q,
    )
    @test legacy.valid

    blocks = Dict{Tuple{Symbol,Int},Vector{Float64}}()
    for representation in (:singlet, :adjoint)
        so3model = build_so3_model(nm1=3, representation=representation)
        workspace = build_workspace(so3model; disp_std=false)
        for ell in 0:2
            hamiltonian = build_hamiltonian(workspace, ell, couplings; disp_std=false)
            blocks[(representation, ell)] = solve(
                hamiltonian; k=3, disp_std=false,
            )
        end
    end
    adapted = score_cft_blocks(blocks)
    @test adapted.ground_is_singlet_l0
    @test adapted.raw_gaps ≈ legacy.raw_gaps atol=2e-11
    @test adapted.factor ≈ legacy.factor atol=2e-11
    @test adapted.q ≈ legacy.q atol=2e-11
    @test adapted.delta_s ≈ legacy.delta_s atol=2e-11
    @test adapted.delta_o ≈ legacy.delta_o atol=2e-11
end
