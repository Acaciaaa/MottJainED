# Exercise the fresh prepare -> five-point sector scan -> collect path at N=3.
# Slurm launches the sector calls independently; this local check uses one thread.
using Test, TOML, CSV, DataFrames
include(joinpath(@__DIR__, "..", "experimental", "FastED.jl"))

@testset "Fresh five-point scout integration" begin
    mktempdir() do directory
        config = TOML.parsefile(joinpath(@__DIR__, "..", "config", "fast_ed",
                                        "n7_uf0_165_k20_scout_restart.toml"))
        config["model"]["nm1"] = 3
        config["fast_ed"]["cache_root"] = joinpath(directory, "cache")
        config["fast_ed"]["result_root"] = joinpath(directory, "results")
        config["fast_ed"]["mus"] = 0.10853896038463 .+ [-0.01, -0.005, 0.0, 0.005, 0.01]
        path = joinpath(directory, "profile.toml")
        open(io -> TOML.print(io, config), path, "w")
        spec = FastED.load_spec(path)
        FastED.prepare_caches(spec; force=true)
        for task_id in 0:19
            FastED.solve_sector(spec, task_id ÷ 4 + 1, task_id % 4 + 1; force=true)
        end
        FastED.collect_all(spec)
        audit = FastED.validate_collection(spec; require_interior=true)
        @test audit.best_mu ≈ 0.10853896038463 atol=1e-12
        scan = CSV.read(joinpath(spec.result_directory, "scan_summary.csv"), DataFrame)
        @test nrow(scan) == 5
        @test all(isfinite, scan.q)
        manifest = TOML.parsefile(joinpath(spec.result_directory, "collection_manifest.toml"))
        @test manifest["complete"] && manifest["collected_mu_count"] == 5
        @test manifest["stage"] == "scout"
        @test manifest["requires_review"] && !manifest["muc_confirmed"]
        @test isdir(spec.cache_directory)
        @test_throws ArgumentError FastED.release_cache(spec)
    end
end
