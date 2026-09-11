using Test
using TOML
using MottJainED

include(joinpath(@__DIR__, "..", "experimental", "FastED.jl"))
using .FastED

@testset "Independent persistent-sector ED" begin
    mktempdir() do directory
        config_path = joinpath(directory, "fast_ed_test.toml")
        config = Dict{String,Any}(
            "model" => Dict("nm1" => 2),
            "hamiltonian" => Dict(
                "Uf" => 0.5,
                "Uf0" => 1.8,
                "U0" => 4.5,
                "Vf" => 0.0,
                "Vf0" => 0.4,
                "V0" => 1.0,
                "t" => 0.5,
                "mu" => 0.03,
            ),
            "solver" => Dict(
                "k" => 8,
                "dense_cutoff" => 128,
                "warm_start" => false,
            ),
            "fast_ed" => Dict(
                "run_name" => "test",
                "cache_root" => joinpath(directory, "cache"),
                "result_root" => joinpath(directory, "results"),
                "k" => 8,
                "mus" => [0.03],
                "score_terms" => ["j"],
                "score_metric" => "q",
            ),
        )
        open(config_path, "w") do io
            TOML.print(io, config; sorted=true)
        end

        spec = FastED.load_spec(config_path)
        manifest = FastED.prepare_caches(spec)
        @test manifest["complete"]
        @test length(manifest["sectors"]) == 4
        mtimes = Dict(
            String(entry["file"]) => mtime(joinpath(spec.cache_directory, String(entry["file"])))
            for entry in manifest["sectors"]
        )
        reused = FastED.prepare_caches(spec)
        @test reused["cache_id"] == manifest["cache_id"]
        @test all(
            mtime(joinpath(spec.cache_directory, file)) == timestamp
            for (file, timestamp) in mtimes
        )

        for sector_index in eachindex(manifest["sectors"])
            output = FastED.solve_sector(spec, 1, sector_index)
            @test isfile(output)
        end
        result = FastED.collect_mu(spec, 1)
        @test isfile(joinpath(FastED.sector_result_path(spec, 0.03, 1, 1)))
        @test isfile(joinpath(spec.result_directory, "mu_$(FastED.mu_id(0.03))", "merged_spectrum.csv"))

        model = build_model(nm1=2)
        direct_cache = prepare_spectrum(model, spec.couplings, spec.solver)
        direct = solve_spectrum(direct_cache, 0.03)
        @test length(result.states) == length(direct)
        @test [state.energy for state in result.states] ≈
              [state.energy for state in direct] atol=1e-11
        @test [state.l2 for state in result.states] ≈
              [state.l2 for state in direct] atol=1e-10
        @test [state.c2 for state in result.states] ≈
              [state.c2 for state in direct] atol=1e-10
        @test [(state.sector.z, state.sector.r, state.rank) for state in result.states] ==
              [(state.sector.z, state.sector.r, state.rank) for state in direct]

        collected = FastED.collect_all(spec)
        @test length(collected.rows) == 1
        @test isempty(collected.failures)
        @test isfile(joinpath(spec.result_directory, "scan_summary.csv"))
        if result.score.valid
            @test isfile(joinpath(spec.result_directory, "best_summary.csv"))
            @test isfile(joinpath(spec.result_directory, "best_relations.csv"))
            @test isfile(joinpath(spec.result_directory, "best_spectrum.csv"))
            collection = TOML.parsefile(joinpath(
                spec.result_directory, "collection_manifest.toml",
            ))
            @test collection["complete"]
            @test collection["best_mu"] == 0.03
        end
    end
end
