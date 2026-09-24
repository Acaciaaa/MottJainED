using Test
using MottJainED
using LinearAlgebra
using Random

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

@testset "Laughlin-1/3 heavy-space projection" begin
    couplings = Couplings(
        Uf=0.37, Uf0=-0.21, U0=0.83,
        Vf=0.19, Vf0=-0.17, V0=0.11,
        t=0.29, mu=-0.07,
    )
    model = build_so3_model(nm1=3, representation=:singlet)
    @test laughlin13_root_counts(model) == [1, 1, 2, 1]
    @test_throws ArgumentError build_workspace(
        model; heavy_space_mode=:unknown, disp_std=false,
    )
    @test_throws ArgumentError build_workspace(
        model; heavy_space_mode=:laughlin13,
        nst_max_heavy=ones(Int, size(model.sec_heavy, 2)), disp_std=false,
    )

    full_workspace = build_workspace(model; disp_std=false)
    projected_workspace = build_workspace(
        model; heavy_space_mode=:laughlin13, disp_std=false,
    )
    @test projected_workspace.heavy_space_mode == :laughlin13
    @test size.(projected_workspace.heavy_space.sts, 2) ==
          laughlin13_root_counts(model)

    for ell in 0:2
        full_hamiltonian = build_hamiltonian(
            full_workspace, ell, couplings; disp_std=false,
        )
        projected_hamiltonian = build_hamiltonian(
            projected_workspace, ell, couplings; disp_std=false,
        )
        full_v1 = Matrix(
            build_heavy_v1_operator(full_hamiltonian; disp_std=false);
            disp_std=false,
        )
        projected_v1 = Matrix(
            build_heavy_v1_operator(projected_hamiltonian; disp_std=false);
            disp_std=false,
        )
        @test opnorm(projected_v1) < 2e-12

        # This is an independent explicit projection: find ker(V1) inside the
        # unprojected composite block, form P H P there, and compare its whole
        # spectrum with the Hamiltonian assembled directly from Jack states.
        v1_decomposition = eigen(Symmetric(full_v1))
        zero_indices = findall(abs.(v1_decomposition.values) .< 2e-10)
        @test length(zero_indices) == sector_dimension(projected_hamiltonian)
        projector_basis = v1_decomposition.vectors[:, zero_indices]
        full_matrix = Matrix(full_hamiltonian.operator; disp_std=false)
        explicit_php = Symmetric(projector_basis' * full_matrix * projector_basis)
        projected_matrix = Symmetric(Matrix(
            projected_hamiltonian.operator; disp_std=false,
        ))
        @test eigvals(explicit_php) ≈ eigvals(projected_matrix) atol=2e-10

        # The low-energy spectrum of H + lambda*V1 must approach P H P when
        # the parent-Hamiltonian penalty is increased.
        projected_gaps = eigvals(projected_matrix)
        projected_gaps .-= first(projected_gaps)
        convergence_errors = Float64[]
        for penalty in (20.0, 200.0, 2000.0)
            finite_v1 = eigvals(Symmetric(full_matrix + penalty * full_v1))[
                1:length(projected_gaps)
            ]
            finite_v1 .-= first(finite_v1)
            push!(convergence_errors, maximum(abs.(finite_v1 - projected_gaps)))
        end
        if convergence_errors[1] > 1e-10
            @test convergence_errors[3] < convergence_errors[2] < convergence_errors[1]
        else
            # Some small blocks are invariant under H already, so every
            # finite penalty agrees with P H P to floating-point precision.
            @test maximum(convergence_errors) < 1e-10
        end
        @test convergence_errors[3] < 2e-3
    end

    # Inside the exact V1 kernel, the old V0 component has no remaining
    # pseudopotential action.  It must therefore be only a common energy shift
    # plus a chemical-potential shift, so V0 is no longer a physical axis.
    redundancy_hamiltonian = build_hamiltonian(
        projected_workspace, 0, isolated_component(:V0); disp_std=false,
    )
    v0_matrix = Matrix(redundancy_hamiltonian.operator; disp_std=false)
    retune!(redundancy_hamiltonian, isolated_component(:mu))
    mu_matrix = Matrix(redundancy_hamiltonian.operator; disp_std=false)
    identity_matrix = Matrix{Float64}(I, size(v0_matrix)...)
    redundancy_coefficients = hcat(
        vec(identity_matrix), vec(mu_matrix),
    ) \ vec(v0_matrix)
    redundancy_residual = norm(
        v0_matrix - redundancy_coefficients[1] * identity_matrix -
        redundancy_coefficients[2] * mu_matrix,
    ) / norm(v0_matrix)
    @test redundancy_residual < 2e-12
end

@testset "Native SO(3)lver conformal generator" begin
    # These values were independently produced by the historical FuzzifiED
    # 1.2.1 full-Fock implementation.  The native calculation uses only the
    # eight SU(3)-singlet tensors but must span exactly the same generated
    # state and descendant projections.
    couplings = Couplings(
        Uf=0.46, U0=4.14, Uf0=1.6605860385658133, Vf=0.0,
        Vf0=0.3449772487682664, V0=0.6208991002533251,
        t=0.5, mu=0.10948169732213274,
    )
    model = build_so3_model(nm1=4, representation=:singlet)
    workspace = build_workspace(model; disp_std=false)
    result = analyze_scalar_generator(
        workspace, couplings;
        l0_count=4, l2_count=3, eig_tol=1.0e-10, disp_std=false,
    )

    @test result.fit.names == GENERATOR_CANDIDATE_NAMES
    @test result.fit.numerical_rank == 8
    @test result.fit.fidelity ≈ 0.999443221677711 atol=2e-12
    @test result.l0_overlap.values[2:3] ≈
          [0.625725483601669, 0.338092599947906] atol=2e-11
    @test sum(result.l0_overlap.values[2:3]) ≈
          0.963818083549575 atol=2e-11
    @test result.l2_overlap.values[1:2] ≈
          [0.201414560552697, 0.702469374066752] atol=2e-11
    @test sum(result.l2_overlap.values[1:2]) ≈
          0.903883934619449 atol=2e-11
end

@testset "Projected native SO(3) conformal algebra" begin
    couplings = Couplings(
        Uf=0.48, U0=4.32, Uf0=1.68, Vf=0.0,
        Vf0=0.3549772487682664, V0=0.0,
        t=0.5, mu=0.22881105175318384,
    )
    singlet_model = build_so3_model(nm1=4, representation=:singlet)
    singlet_workspace = build_workspace(
        singlet_model; heavy_space_mode=:laughlin13, disp_std=false,
    )
    adjoint_model = build_so3_model(nm1=4, representation=:adjoint)
    adjoint_workspace = build_workspace(
        adjoint_model;
        heavy_space=singlet_workspace.heavy_space,
        heavy_space_mode=:laughlin13,
        disp_std=false,
    )
    prepared_problem = build_so3_conformal_problem(
        singlet_workspace, adjoint_workspace, couplings; disp_std=false,
    )
    result = analyze_so3_conformal_algebra(
        singlet_workspace, adjoint_workspace, couplings;
        block_counts=Dict(
            (:singlet, 0) => 3,
            (:singlet, 2) => 2,
            (:adjoint, 0) => 2,
            (:adjoint, 1) => 2,
        ),
        factor_bounds=(0.005, 0.2),
        eig_tol=1.0e-9,
        ncv=12,
        prepared_problem,
        disp_std=false,
    )

    @test result.heavy_space_mode == :laughlin13
    @test result.fit.factor_bounds[1] <= result.fit.factor <= result.fit.factor_bounds[2]
    @test isfinite(result.fit.value)
    @test result.fit.value >= 0
    @test result.fit.normalization_rank >= 1
    @test Set(getproperty.(result.primary_rows, :label)) == Set((:S, :O, :J, :T))
    @test all(row -> isfinite(row.k_fraction) && row.k_fraction >= 0,
              result.primary_rows)
    @test Set(keys(result.k2)) == Set((
        (:singlet, 0), (:singlet, 2), (:adjoint, 0), (:adjoint, 1),
    ))
    @test all(sector -> issorted(sector.eigenvalues), values(result.k2))
    @test length(result.mixed_commutator_rows) == 4
    @test length(result.mixed_commutator_pair_rows) == 36
    @test all(row -> isfinite(row.fractional_residual) &&
                     row.fractional_residual >= 0 &&
                     isfinite(row.p_commutator_fraction) &&
                     row.p_commutator_fraction >= 0 &&
                     isfinite(row.k_commutator_fraction) &&
                     row.k_commutator_fraction >= 0,
              result.mixed_commutator_rows)
    @test haskey(result.dimensions, (:singlet, 4))
    @test haskey(result.dimensions, (:adjoint, 3))
    @test !haskey(result.energies, (:singlet, 4))
    @test !haskey(result.energies, (:adjoint, 3))
    @test !isempty(prepared_problem.operator_cache)
    objective = score_so3_conformal_algebra(
        result; labels=[:S, :O, :J], worst_weight=0.2,
    )
    @test isfinite(objective.objective) && objective.objective >= 0
    @test objective.labels == [:J, :O, :S]
    @test Set(getproperty.(objective.rows, :term)) ==
          Set(CONFORMAL_OBJECTIVE_TERMS)
    holdout = score_so3_conformal_algebra(
        result; labels=[:T], worst_weight=0.0,
    )
    @test holdout.labels == [:T]
    @test_throws ArgumentError score_so3_conformal_algebra(
        result; labels=[:missing],
    )

    # The full magnetic-substate reconstruction must agree with the simpler
    # scalar reduced-matrix-element identity used for generator normalization.
    for primary in filter(row -> row.ell == 0, result.primary_rows)
        diagonal_pairs = filter(
            row -> row.label == primary.label &&
                   row.first_axis == row.second_axis,
            result.mixed_commutator_pair_rows,
        )
        @test length(diagonal_pairs) == 3
        @test all(
            row -> isapprox(
                row.lhs_diagonal_expectation, primary.kp_commutator_lhs;
                atol=2e-10, rtol=2e-10,
            ),
            diagonal_pairs,
        )
        @test all(
            row -> isapprox(
                row.rhs_diagonal_expectation, primary.kp_commutator_target;
                atol=2e-10, rtol=2e-10,
            ),
            diagonal_pairs,
        )
    end

    # P and K are formed from exact cross-block Hamiltonian actions, not from
    # a low-energy spectral sum.  Verify their defining identities directly.
    h0 = build_hamiltonian(singlet_workspace, 0, couplings; disp_std=false)
    h1 = build_hamiltonian(singlet_workspace, 1, couplings; disp_std=false)
    _, states0 = solve(h0; k=2, vectors=true, disp_std=false)
    operator_set = build_generator_operators(
        h0.space, h1.space, build_generator_candidates(singlet_model);
        disp_std=false,
    )
    source = view(states0, :, 2)
    lambda = apply_so3_conformal_generator(
        source, operator_set, result.fit.coefficients, h0, h1;
        factor=result.fit.factor, generator=:lambda,
    )
    p_action = apply_so3_conformal_generator(
        source, operator_set, result.fit.coefficients, h0, h1;
        factor=result.fit.factor, generator=:p,
    )
    k_action = apply_so3_conformal_generator(
        source, operator_set, result.fit.coefficients, h0, h1;
        factor=result.fit.factor, generator=:k,
    )
    commutator = (
        h1.operator * lambda -
        apply_so3_generator(h0.operator * Vector(source), operator_set,
                            result.fit.coefficients)
    ) ./ result.fit.factor
    @test p_action + k_action ≈ lambda atol=2e-11 rtol=2e-11
    @test p_action - k_action ≈ commutator atol=2e-11 rtol=2e-11
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
    stable_six = score_cft_blocks(blocks; terms=CFT_STABLE_SIX_TERMS)
    @test stable_six.q < 1e-14
    @test stable_six.raw_gaps ≈ [1, 1, 2, 3, 3, 3]
    @test stable_six.terms == collect(CFT_STABLE_SIX_TERMS)
    holdout = score_cft_blocks(blocks; terms=[:boxo_o, :boxj_j])
    @test holdout.q < 1e-14
    @test holdout.raw_gaps ≈ [2, 2]
    @test_throws ArgumentError score_cft_blocks(blocks; terms=Symbol[])
    @test_throws ArgumentError score_cft_blocks(blocks; terms=[:j, :j])
    @test_throws ArgumentError score_cft_blocks(blocks; terms=[:not_a_relation])
    @test_throws ArgumentError score_cft_blocks(delete!(copy(blocks), (:adjoint, 0)))
end

@testset "SO(3)lver parameter-search identity guards" begin
    stable_labels = getproperty.(STABLE_SIX_TRACKED_STATE_SPECS, :label)
    @test :ddS in stable_labels
    @test :boxS ∉ stable_labels

    lhs = latin_hypercube_points(
        8, [1.0, -2.0], [3.0, 2.0], MersenneTwister(17),
    )
    @test size(lhs) == (8, 2)
    @test all(1.0 .<= lhs[:, 1] .<= 3.0)
    @test all(-2.0 .<= lhs[:, 2] .<= 2.0)
    for column in 1:2
        lower, upper = column == 1 ? (1.0, 3.0) : (-2.0, 2.0)
        bins = floor.(Int, 8 .* (lhs[:, column] .- lower) ./ (upper - lower))
        @test sort(bins) == collect(0:7)
    end
    @test lhs == latin_hypercube_points(
        8, [1.0, -2.0], [3.0, 2.0], MersenneTwister(17),
    )

    brackets = mu_refinement_brackets(
        collect(0.0:0.1:0.6),
        [4.0, 1.0, 3.0, 2.0, 0.5, 2.0, 4.0],
        trues(7); maximum_count=2,
    )
    @test getproperty.(brackets, :index) == [5, 2]
    @test brackets[1].lower ≈ 0.3
    @test brackets[1].upper ≈ 0.5
    invalid_middle = mu_refinement_brackets(
        collect(0.0:0.1:0.4), [3.0, 1.0, 2.0, 0.5, 3.0],
        Bool[true, true, false, true, true],
    )
    @test isempty(invalid_middle)

    center = [1.834, 0.41, 0.525, 0.1218143040052]
    scales = [0.05, 0.025, 0.01, 0.0015]
    search_point = [-1.0, 0.5, 2.0, -0.25]
    physical_point = parameter_search_to_physical(search_point, center, scales)
    @test physical_point ≈ [1.784, 0.4225, 0.545, 0.1214393040052]
    @test physical_to_parameter_search(physical_point, center, scales) ≈ search_point
    @test_throws DimensionMismatch parameter_search_to_physical([1.0], center, scales)
    @test_throws ArgumentError physical_to_parameter_search(center, center, -scales)

    linked_grid = linked_uf_grid_values(
        [0.44, 0.46], [1.64, 1.66], [0.345], [0.623, 0.64];
        u0_over_uf=9.0,
    )
    @test length(linked_grid) == 8
    @test all(point -> point.U0 ≈ 9 * point.Uf, linked_grid)
    @test first(linked_grid).Uf == 0.44
    @test first(linked_grid).U0 ≈ 3.96
    @test first(linked_grid).Uf0 == 1.64
    @test first(linked_grid).Vf0 == 0.345
    @test first(linked_grid).V0 == 0.623
    @test last(linked_grid).Uf == 0.46
    @test last(linked_grid).U0 ≈ 4.14
    @test last(linked_grid).Uf0 == 1.66
    @test last(linked_grid).Vf0 == 0.345
    @test last(linked_grid).V0 == 0.64
    @test_throws ArgumentError linked_uf_grid_values(
        [0.46, 0.46], [1.64], [0.345], [0.623],
    )
    @test_throws ArgumentError linked_uf_grid_values(
        [0.46], [1.64], [0.345], [0.623]; u0_over_uf=0,
    )

    generator_gate = assess_scalar_generator_overlaps(
        [0.02, 0.42, 0.40, 0.03, 0.02],
        [0.43, 0.38, 0.04, 0.02],
    )
    @test generator_gate.passed
    @test generator_gate.l0_unresolved_overlap_upper_bound ≈ 0.11
    @test generator_gate.l2_unresolved_overlap_upper_bound ≈ 0.13
    @test generator_gate.boxs_leading_margin ≈ 0.29
    @test generator_gate.dds_leading_margin ≈ 0.25
    truncated_failure = assess_scalar_generator_overlaps(
        [0.02, 0.30, 0.30, 0.03], [0.43, 0.38, 0.04];
        minimum_expected_subspace_overlap=0.0,
        minimum_expected_subspace_fraction=0.0,
    )
    @test !truncated_failure.passed
    @test !truncated_failure.boxs_is_leading_non_s
    @test_throws ArgumentError assess_scalar_generator_overlaps(
        [0.1, 0.4, 0.4], [0.4, 0.4, 0.1],
    )

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
