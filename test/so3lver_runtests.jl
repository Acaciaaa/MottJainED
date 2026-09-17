using Test
using MottJainED

include(joinpath(@__DIR__, "..", "experimental", "SO3lverED.jl"))
using .SO3lverED

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
    for (ell, hamiltonian) in enumerate(singlet_hamiltonians)
        physical_ell = ell - 1
        reference = reference_energies(catalog, physical_ell, 0)
        @test sector_dimension(hamiltonian) == length(reference)
        @test solve(hamiltonian; k=100, disp_std=false) ≈ reference atol=2e-11
    end

    adjoint_model = build_so3_model(nm1=2, representation=:adjoint)
    @test adjoint_model.weight == (1, 3)
    @test adjoint_model.c2 == 3.0
    adjoint_workspace = build_workspace(adjoint_model; disp_std=false)
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
