# Historical automatic-search prototype regression; not the current Slurm path.
# Does not construct any N=7 matrix or touch production cache/result directories.
using Test, TOML, CSV, DataFrames
include(joinpath(@__DIR__, "..", "experimental", "FastED.jl"))
include(joinpath(@__DIR__, "..", "experimental", "FastMuSearch.jl"))

@testset "Sequential cached mu-search integration" begin
    mktempdir() do directory
        config = TOML.parsefile(joinpath(@__DIR__, "..", "config", "fast_ed", "n7_uf0_165_search.toml"))
        config["model"]["nm1"] = 3
        config["fast_ed"]["cache_root"] = joinpath(directory, "cache")
        config["fast_ed"]["result_root"] = joinpath(directory, "results")
        config["fast_ed"]["run_name"] = "integration"
        merge!(config["mu_search"], Dict(
            "mu_min"=>0.08, "mu_max"=>0.15, "wide_count"=>5,
            "guard_min"=>0.09, "guard_max"=>0.145, "guard_count"=>5,
            "seed_mu"=>0.12, "local_half_width"=>0.01, "local_count"=>5,
            "mu_abs_tol"=>0.00002, "max_evaluations"=>100,
        ))
        path = joinpath(directory, "profile.toml")
        open(io -> TOML.print(io, config), path, "w")
        spec = FastED.load_spec(path)
        result = FastMuSearch.run(spec)
        @test result["complete"] && result["audit_passed"]
        @test result["execution_mode"] == "serial_sectors" && result["solver_processes"] == 1
        @test !result["global_minimum_proven"] && !result["cache_released"]
        @test abs(result["best_mu"]-0.10853896038463) < 0.0001
        @test isdir(spec.cache_directory)
        trace = CSV.read(joinpath(spec.result_directory, "mu_search_evaluations.csv"), DataFrame)
        @test nrow(trace) == result["evaluated_mu_count"]
        # N3 sector dimensions are 42, 30, 15, 19: k=20 yields 74 states.
        @test all(trace.state_count .== 74)
        final = FastED.load_spec(joinpath(spec.result_directory, "evaluated_profile.toml"))
        @test final.cache_id == spec.cache_id && final.settings_id == spec.settings_id
        @test FastED.validate_collection(final; require_interior=true).best_mu == result["best_mu"]
        @test_throws ArgumentError FastED.release_cache(final)
        # Replay the same search: complete sector CSVs are reused, apart from the
        # intentional cold repeat at the winner.
        probe_mu = first(final.mus)
        probe_file = FastED.sector_result_path(spec, probe_mu, 1, 1)
        before = mtime(probe_file)
        repeated = FastMuSearch.run(spec)
        @test repeated["audit_passed"]
        @test mtime(probe_file) == before
        @test repeated["best_mu"] ≈ result["best_mu"] atol=2e-5
    end
end
