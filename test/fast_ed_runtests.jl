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

        # A resident sector must give the unchanged spectrum at a new mu, and reject
        # a truncated CSV instead of silently treating a partial result as complete.
        resident = FastED.load_sector_solver(spec, 1)
        resident_path = FastED.solve_sector(spec, 0.031, 1; resident=resident)
        resident_rows = FastED.CSV.read(resident_path, FastED.DataFrame)
        reference_states = solve_spectrum(direct_cache, 0.031)
        reference_energies = [s.energy for s in reference_states
            if s.sector.z == resident.z && s.sector.r == resident.r]
        @test resident_rows.energy ≈ reference_energies atol=1e-11
        @test FastED.result_is_current(resident_path, spec, 0.031, resident.z, resident.r)
        MottJainED.atomic_csv(resident_path, resident_rows[1:end-1, :])
        @test !FastED.result_is_current(resident_path, spec, 0.031, resident.z, resident.r)
        FastED.solve_sector(spec, 0.031, 1; resident=resident)
        @test FastED.result_is_current(resident_path, spec, 0.031, resident.z, resident.r)
        @test_throws ArgumentError FastED.solve_sector(spec, 0.032, 2; resident=resident)

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
        validation = FastED.validate_collection(spec)
        @test validation.best_mu == 0.03
        @test_throws ArgumentError FastED.release_cache(spec)

        config["fast_ed"]["allow_cache_release"] = true
        open(config_path, "w") do io
            TOML.print(io, config; sorted=true)
        end
        releasable = FastED.load_spec(config_path)
        released = FastED.release_cache(releasable)
        @test released.cache_id == spec.cache_id
        @test released.bytes > 0
        @test !isdir(spec.cache_directory)
        @test isfile(joinpath(spec.result_directory, "best_summary.csv"))
    end
end

include(joinpath(@__DIR__, "..", "experimental", "FastMuSearch.jl"))

@testset "Cached N7 mu search finds competing branches" begin
    profile = joinpath(@__DIR__, "..", "config", "fast_ed", "n7_uf0_165_search.toml")
    spec = FastED.load_spec(profile)
    opt = FastMuSearch.options(spec)
    @test spec.nm1 == 7 && spec.solver.k == 20 && !spec.solver.warm_start
    @test spec.couplings.Uf0 == 1.65 && spec.couplings.Vf0 == 0.55 && spec.couplings.V0 == 0.34
    @test !spec.config["fast_ed"]["allow_cache_release"]
    test_opt = merge(opt, (budget=300, tolerance=1e-7))
    calls = Dict{Float64,Int}()
    objective(mu) = begin
        calls[mu] = get(calls, mu, 0) + 1
        (valid=true, objective=(mu-0.1447)^2+0.05)
    end
    result = FastMuSearch.search_mu(objective, test_opt)
    @test result.best_mu ≈ 0.1447 atol=2e-6
    @test isempty(result.issues)
    @test all(==(1), values(calls)) # shared mu evaluations do not redo ED
    @test all(MottJainED._mu_key(mu) in keys(calls) for mu in
              range(opt.lower, opt.upper; length=opt.wide_count))
    @test all(MottJainED._mu_key(mu) in keys(calls) for mu in
              range(opt.guard_lower, opt.guard_upper; length=opt.guard_count))

    # The historical failure: the seed and whole-range Brent agree near .142,
    # while an independently sampled narrow branch near .10 is lower.
    hidden(mu) = (valid=true, objective=abs(mu-0.101) < 0.005 ?
                 0.01 + 50*(mu-0.101)^2 : 0.1 + (mu-0.142)^2)
    rival = FastMuSearch.search_mu(hidden, test_opt)
    @test rival.best_mu ≈ 0.101 atol=2e-6
    @test rival.score.objective < rival.local_result.score.objective
    @test "local_and_independent_grid_disagree" in rival.issues

    # Missing score near the candidate cannot silently pass; an outer missing score
    # is retained as an explicit coverage limit, never a proof of a global optimum.
    hole(mu) = (valid=abs(mu-0.12)>0.0001, objective=(mu-0.1447)^2+0.05)
    with_hole = FastMuSearch.search_mu(hole, test_opt)
    @test "unscored_point_inside_dense_guard" in with_hole.issues
    @test !isempty(with_hole.invalid)
    boundary = FastMuSearch.search_mu(mu -> (valid=true, objective=mu+1), test_opt)
    @test "winner_outside_dense_guard" in boundary.issues
    @test !boundary.neighbor_converged
    @test_throws ErrorException FastMuSearch.search_mu(objective, merge(test_opt, (budget=3,)))
end

@testset "N=7 fine-scan profile" begin
    config_root = joinpath(@__DIR__, "..", "config", "fast_ed")
    scout = FastED.load_spec(joinpath(config_root, "n7_retained_k20.toml"))
    refine = FastED.load_spec(joinpath(config_root, "n7_retained_k20_refine.toml"))

    @test refine.nm1 == 7
    @test refine.solver.k == 20
    @test refine.cache_id == scout.cache_id
    @test refine.run_name == "n7_retained_refine"
    @test length(refine.mus) == 7
    @test first(refine.mus) == 0.14450
    @test last(refine.mus) == 0.14600
    @test all(isapprox.(diff(refine.mus), 0.00025; atol=1.0e-14, rtol=0.0))
end

@testset "Stage-12 guided local-FSS profiles" begin
    profile_root = joinpath(@__DIR__, "..", "config", "two_size_tuning")
    reference = TOML.parsefile(joinpath(profile_root, "s_stage12_vf_retained.toml"))
    search_keys = filter(key -> startswith(key, "optimize_") || startswith(key, "mu_") ||
        startswith(key, "score_") || key == "k", keys(reference["two_size_tuning"]))
    hamiltonian_points = Set{Tuple}()
    cases = (
        ("n56_retained_local_uf0.toml", "Uf0", [1.65, 1.834, 2.00]),
        ("n56_retained_local_vf0.toml", "Vf0", [0.45, 0.65]),
    )
    for (name, parameter, values) in cases
        config = TOML.parsefile(joinpath(profile_root, name))
        tuning = config["two_size_tuning"]
        @test tuning["guide_nm_values"] == [3, 4]
        @test tuning["match_nm_values"] == [5, 6]
        @test tuning["scan_parameter"] == parameter
        @test tuning["scan_values"] == values
        @test config["hamiltonian"] == reference["hamiltonian"]
        for key in search_keys
            @test tuning[key] == reference["two_size_tuning"][key]
        end
        for value in values
            point = parameter == "Uf0" ? (value, 0.55) : (1.834, value)
            @test !(point in hamiltonian_points)
            push!(hamiltonian_points, point)
        end
    end
    @test length(hamiltonian_points) == 5
    @test (1.834, 0.55) in hamiltonian_points
end
