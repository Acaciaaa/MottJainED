using Test
using MottJainED

@testset "Couplings" begin
    base = Couplings()
    changed = MottJainED.with_coupling(base, :Uf0, 2.5)
    @test changed.Uf0 == 2.5
    @test changed.Uf == base.Uf
    @test_throws ArgumentError MottJainED.validate(Couplings(mu=NaN))
end

function state(energy, l2, c2; z=1, r=1)
    return SpectrumState(energy, l2, c2, SectorKey(z, r), 1, nothing, nothing)
end

@testset "CFT tower classification" begin
    states = SpectrumState[
        state(0.0, 0, 0), state(1.2, 0, 0), state(3.2, 0, 0),
        state(2.2, 2, 0), state(3.0, 6, 0), state(3.2, 6, 0),
        state(2.0, 2, 3; z=1), state(2.0, 2, 3; z=-1),
        state(3.0, 2, 3; z=1), state(3.0, 2, 3; z=-1),
        state(3.0, 6, 3), state(1.6, 0, 3),
    ]
    catalog, rejected = level_catalog(states)
    @test isempty(rejected)
    @test length(catalog[(2, 3)]) == 2
    @test catalog[(2, 3)][1].multiplicity == 2
    score = cft_score(states)
    @test score.valid
    @test score.q < 1e-12
    @test score.factor ≈ 1.0
    @test score.delta_s ≈ 1.2
    @test score.delta_o ≈ 1.6
end

@testset "Configuration" begin
    config = load_config(joinpath(dirname(@__DIR__), "config", "default.toml"))
    @test config["model"]["nm1"] == 5
end

@testset "Cached Hamiltonian is exact" begin
    couplings = Couplings(mu=0.03)
    settings = SolverSettings(k=8)
    model = build_model(nm1=2)
    cache = prepare_spectrum(model, couplings, settings)
    sector = first(cache.sectors)
    direct = MottJainED.float_opmat(
        MottJainED.FuzzifiED.Operator(
            sector.basis, hamiltonian_terms(model, couplings),
        ),
    )
    cached = MottJainED.hermitian_opmat(
        sector.h0 + couplings.mu * sector.number_f,
    )
    @test Matrix(direct) ≈ Matrix(cached) atol=1e-13
    states = solve_spectrum(cache, couplings.mu)
    @test first(states).l2 ≈ 0.0 atol=1e-10
    @test first(states).c2 ≈ 0.0 atol=1e-10
end
