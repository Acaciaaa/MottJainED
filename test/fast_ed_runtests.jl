using Test
using TOML
using SHA
using MottJainED

include(joinpath(@__DIR__, "..", "experimental", "FastED.jl"))
using .FastED
include(joinpath(@__DIR__, "..", "experimental", "FastEDPipeline.jl"))
using .FastEDPipeline

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
        rebuilt = FastED.prepare_caches(spec; force=true)
        @test rebuilt["cache_id"] == manifest["cache_id"]
        @test all(mtime(joinpath(spec.cache_directory, file)) > timestamp
                  for (file, timestamp) in mtimes)

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
        completed_at = mtime(resident_path)
        FastED.solve_sector(spec, 0.031, 1; force=true)
        @test mtime(resident_path) > completed_at
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

@testset "Bounded N7 automatic pipeline decisions and resources" begin
    root = joinpath(@__DIR__, "..")
    pilot_path = joinpath(root, "config", "fast_ed", "n7_uf0_200_auto.toml")
    retirement_path = joinpath(root, "config", "fast_ed", "n7_uf0_165_cache_retirement.toml")
    completed_retirement_path = joinpath(root, "config", "fast_ed", "n7_uf0_200_cache_retirement.toml")
    next_path = joinpath(root, "config", "fast_ed", "n7_vf0_045_auto.toml")
    vf0_retirement_path = joinpath(root, "config", "fast_ed", "n7_vf0_045_cache_retirement.toml")
    final_vf0_path = joinpath(root, "config", "fast_ed", "n7_vf0_065_auto.toml")
    final_vf0_retirement_path = joinpath(root, "config", "fast_ed",
                                         "n7_vf0_065_cache_retirement.toml")
    original_path = joinpath(root, "config", "fast_ed", "n7_original_audit_auto.toml")
    original_uf165_path = joinpath(root, "config", "fast_ed",
                                   "n7_original_uf0_165_auto.toml")
    original_uf165_retirement_path = joinpath(root, "config", "fast_ed",
                                              "n7_original_uf0_165_cache_retirement.toml")
    original_vf052_path = joinpath(root, "config", "fast_ed",
                                   "n7_original_vf0_052_auto.toml")
    pilot = FastED.load_spec(pilot_path)
    opt = FastEDPipeline.pipeline_options(pilot)

    damped = opt.n6_mu + 0.5*(opt.n6_mu-opt.n5_mu)
    @test pilot.nm1 == 7 && pilot.solver.k == 20 && !pilot.solver.warm_start
    @test pilot.couplings.Uf0 == 2.0 && pilot.couplings.Vf0 == 0.55
    @test pilot.mus ≈ damped .+ [-0.01, -0.005, 0.0, 0.005, 0.01] atol=1e-14
    @test pilot.terms == [:ds_s, :j, :curlj, :dj_rank1, :t_rank1]
    @test pilot.metric == :q && FastED.plan(pilot).tasks == 20
    @test FastEDPipeline.resource_fields(pilot_path) == (8, 8, 8, 8, 4, 1)
    @test opt.max_total_mus == 16 && opt.max_adaptive_rounds == 3 && !opt.auto_release
    @test opt.scan_parameter == "Uf0" && opt.scan_value == 2.0
    @test opt.n5_delta_s == 1.3736311582780456
    @test opt.n6_delta_o == 3.0744735460523773
    retirement = FastED.load_spec(retirement_path)
    retirement_config = TOML.parsefile(retirement_path)
    prior_scout = FastED.load_spec(joinpath(root, "config", "fast_ed",
                                            "n7_uf0_165_k20_scout_restart.toml"))
    @test retirement_config["cache_retirement"]["expected_cache_id"] == "aaa8ccb4a8dd9fef"
    @test retirement.cache_id == prior_scout.cache_id
    @test retirement.cache_id != pilot.cache_id

    completed_retirement = FastED.load_spec(completed_retirement_path)
    completed_retirement_config = TOML.parsefile(completed_retirement_path)
    @test completed_retirement.cache_id == pilot.cache_id
    @test completed_retirement_config["cache_retirement"]["expected_cache_id"] ==
          "ce74ad933acda136"
    @test completed_retirement_config["cache_retirement"]["verified_local_archive_sha256"] ==
          "bedb6a2f9d444be294fff1d6f3753ba42b7dbf326df3bcd269e6f211f718ae36"

    next = FastED.load_spec(next_path)
    next_opt = FastEDPipeline.pipeline_options(next)
    next_damped = next_opt.n6_mu + 0.5*(next_opt.n6_mu-next_opt.n5_mu)
    @test next.nm1 == 7 && next.solver.k == 20 && !next.solver.warm_start
    @test next.couplings.Uf0 == 1.834 && next.couplings.Vf0 == 0.45
    @test next.mus ≈ next_damped .+ [-0.01, -0.005, 0.0, 0.005, 0.01] atol=1e-14
    @test next.terms == pilot.terms && next.metric == :q && FastED.plan(next).tasks == 20
    @test FastEDPipeline.resource_fields(next_path) == (8, 8, 8, 8, 4, 1)
    @test next_opt.scan_parameter == "Vf0" && next_opt.scan_value == 0.45
    @test next_opt.n5_q == 0.21819656826690698
    @test next_opt.n6_delta_o == 3.0438948900444003
    @test next.cache_id != pilot.cache_id
    next_config = TOML.parsefile(next_path)
    @test next_config["pipeline"]["retire_before_start"] == [
        "config/fast_ed/n7_uf0_165_cache_retirement.toml",
        "config/fast_ed/n7_uf0_200_cache_retirement.toml",
    ]

    vf0_retirement = FastED.load_spec(vf0_retirement_path)
    vf0_retirement_config = TOML.parsefile(vf0_retirement_path)
    @test vf0_retirement.cache_id == next.cache_id
    @test vf0_retirement_config["cache_retirement"]["expected_cache_id"] ==
          "342e15745d81a91a"
    @test vf0_retirement_config["cache_retirement"]["verified_local_archive_sha256"] ==
          "f4f75f0e9248244e802b944f99d4f6a1bb3548739d57ebf4abc6a4d51d661372"

    final_vf0 = FastED.load_spec(final_vf0_path)
    final_opt = FastEDPipeline.pipeline_options(final_vf0)
    final_damped = final_opt.n6_mu + 0.5*(final_opt.n6_mu-final_opt.n5_mu)
    final_main_scout = final_damped .+ [-0.01, -0.005, 0.0, 0.005, 0.01]
    @test final_vf0.nm1 == 7 && final_vf0.solver.k == 20 && !final_vf0.solver.warm_start
    @test final_vf0.couplings.Uf0 == 1.834 && final_vf0.couplings.Vf0 == 0.65
    expected_final_scout = vcat(final_main_scout, final_opt.guard_mus)
    @test length(final_vf0.mus) == length(expected_final_scout)
    @test all(any(isapprox(expected, actual; atol=1e-14, rtol=0)
                  for actual in final_vf0.mus) for expected in expected_final_scout)
    @test final_opt.guard_mus == [0.12771774983325]
    @test final_vf0.terms == pilot.terms && final_vf0.metric == :q
    @test FastED.plan(final_vf0).tasks == 24
    @test FastEDPipeline.resource_fields(final_vf0_path) == (8, 8, 8, 8, 4, 1)
    @test final_opt.scan_parameter == "Vf0" && final_opt.scan_value == 0.65
    @test final_opt.n5_q == 0.15935728084504894
    @test final_opt.n6_q == 0.11603076339853687
    @test final_opt.n6_delta_o == 2.9429301053671573
    @test final_opt.max_total_mus == 16 && final_opt.mu_min == 0.11
    @test final_vf0.cache_id != next.cache_id
    final_config = TOML.parsefile(final_vf0_path)
    @test final_config["pipeline"]["retire_before_start"] == [
        "config/fast_ed/n7_uf0_165_cache_retirement.toml",
        "config/fast_ed/n7_uf0_200_cache_retirement.toml",
        "config/fast_ed/n7_vf0_045_cache_retirement.toml",
    ]

    final_vf0_retirement = FastED.load_spec(final_vf0_retirement_path)
    final_vf0_retirement_config = TOML.parsefile(final_vf0_retirement_path)
    @test final_vf0_retirement.cache_id == final_vf0.cache_id
    @test final_vf0_retirement_config["cache_retirement"]["expected_cache_id"] ==
          "acae800d795e294b"
    @test final_vf0_retirement_config["cache_retirement"]["verified_local_archive_sha256"] ==
          "e2f3e2f6545e8a7026ff402e8fb06271d3a8650e990751b5e6a27f08244f1d77"

    original = FastED.load_spec(original_path)
    original_opt = FastEDPipeline.pipeline_options(original)
    original_damped = original_opt.n6_mu + 0.5*(original_opt.n6_mu-original_opt.n5_mu)
    @test original.nm1 == 7 && original.solver.k == 20 && !original.solver.warm_start
    @test original.couplings == FastED.Couplings(
        Uf=0.46, Uf0=1.834, U0=4.14, Vf=0.0, Vf0=0.41, V0=0.525,
        t=0.5, mu=original_damped,
    )
    @test original.mus ≈ original_damped .+ [-0.01, -0.005, 0.0, 0.005, 0.01] atol=1e-14
    @test original.terms == pilot.terms && original.metric == :q
    @test FastED.plan(original).tasks == 20
    @test FastEDPipeline.resource_fields(original_path) == (8, 8, 8, 8, 4, 1)
    @test original_opt.scan_parameter == "Vf0" && original_opt.scan_value == 0.41
    @test original_opt.n5_q == 0.11755568086489525
    @test original_opt.n6_q == 0.06517470973203797
    @test original_opt.n5_delta_s == 1.4611807307943308
    @test original_opt.n6_delta_o == 3.0224153568005727
    @test isempty(original_opt.guard_mus)
    @test original_opt.max_total_mus == 16
    @test original_opt.mu_min == 0.08 && original_opt.mu_max == 0.17
    @test original.cache_id != final_vf0.cache_id
    original_config = TOML.parsefile(original_path)
    @test original_config["pipeline"]["retire_before_start"] == [
        "config/fast_ed/n7_uf0_165_cache_retirement.toml",
        "config/fast_ed/n7_uf0_200_cache_retirement.toml",
        "config/fast_ed/n7_vf0_045_cache_retirement.toml",
        "config/fast_ed/n7_vf0_065_cache_retirement.toml",
    ]

    original_uf165 = FastED.load_spec(original_uf165_path)
    original_uf165_opt = FastEDPipeline.pipeline_options(original_uf165)
    original_uf165_damped = original_uf165_opt.n6_mu +
                            0.5*(original_uf165_opt.n6_mu-original_uf165_opt.n5_mu)
    @test original_uf165.nm1 == 7 && original_uf165.solver.k == 20
    @test !original_uf165.solver.warm_start
    @test original_uf165.couplings == FastED.Couplings(
        Uf=0.46, Uf0=1.65, U0=4.14, Vf=0.0, Vf0=0.41, V0=0.525,
        t=0.5, mu=0.12139298275993,
    )
    @test original_uf165.couplings.mu ≈ original_uf165_damped atol=1e-14
    @test original_uf165.mus ≈
          original_uf165_damped .+ [-0.01, -0.005, 0.0, 0.005, 0.01] atol=1e-14
    @test original_uf165.terms == original.terms && original_uf165.metric == :q
    @test FastED.plan(original_uf165).tasks == 20
    @test FastEDPipeline.resource_fields(original_uf165_path) == (8, 8, 8, 8, 4, 1)
    @test original_uf165_opt.scan_parameter == "Uf0"
    @test original_uf165_opt.scan_value == 1.65
    @test original_uf165_opt.n5_mu == 0.12126476635669
    @test original_uf165_opt.n6_mu == 0.12135024395885
    @test original_uf165_opt.n5_q == 0.0999568257332489
    @test original_uf165_opt.n6_q == 0.06254640447215823
    @test original_uf165_opt.n5_factor == 0.03790095452441143
    @test original_uf165_opt.n6_factor == 0.03412172635806117
    @test original_uf165_opt.n5_delta_s == 1.5793229163256086
    @test original_uf165_opt.n6_delta_o == 2.9532304784496852
    @test isempty(original_uf165_opt.guard_mus)
    @test original_uf165_opt.max_total_mus == 16
    @test original_uf165_opt.mu_min == 0.08 && original_uf165_opt.mu_max == 0.17
    @test original_uf165.cache_id != original.cache_id
    @test original_uf165.cache_id != prior_scout.cache_id
    original_uf165_config = TOML.parsefile(original_uf165_path)
    @test original_uf165_config["fast_ed"]["run_name"] ==
          "n7_original_uf0_165_auto_scout"
    @test !original_uf165_config["fast_ed"]["allow_cache_release"]
    @test isempty(original_uf165_config["pipeline"]["retire_before_start"])
    @test !original_uf165_config["pipeline"]["auto_release"]

    original_uf165_retirement = FastED.load_spec(original_uf165_retirement_path)
    original_uf165_retirement_config = TOML.parsefile(original_uf165_retirement_path)
    @test original_uf165_retirement.couplings == FastED.Couplings(
        Uf=0.46, Uf0=1.65, U0=4.14, Vf=0.0, Vf0=0.41, V0=0.525,
        t=0.5, mu=0.12125,
    )
    @test original_uf165_retirement_config["cache_retirement"]["expected_cache_id"] ==
          "3bfe378c1c294c11"
    @test original_uf165_retirement_config["cache_retirement"]["verified_local_archive_sha256"] ==
          "0a0479121496c7e05e4011f3e1835a273a4495af420e713b34c112396a56a185"

    original_vf052 = FastED.load_spec(original_vf052_path)
    original_vf052_opt = FastEDPipeline.pipeline_options(original_vf052)
    original_vf052_damped = original_vf052_opt.n6_mu +
                            0.5*(original_vf052_opt.n6_mu-original_vf052_opt.n5_mu)
    @test original_vf052.nm1 == 7 && original_vf052.solver.k == 20
    @test !original_vf052.solver.warm_start
    @test original_vf052.couplings == FastED.Couplings(
        Uf=0.46, Uf0=1.834, U0=4.14, Vf=0.0, Vf0=0.52, V0=0.525,
        t=0.5, mu=0.112517687151475,
    )
    @test original_vf052.couplings.mu ≈ original_vf052_damped atol=1e-14
    @test original_vf052.mus ≈
          original_vf052_damped .+ [-0.01, -0.005, 0.0, 0.005, 0.01] atol=1e-14
    @test original_vf052.terms == original.terms && original_vf052.metric == :q
    @test FastED.plan(original_vf052).tasks == 20
    @test FastEDPipeline.resource_fields(original_vf052_path) == (8, 8, 8, 8, 4, 1)
    @test original_vf052_opt.scan_parameter == "Vf0"
    @test original_vf052_opt.scan_value == 0.52
    @test original_vf052_opt.n5_mu == 0.11307786788428
    @test original_vf052_opt.n6_mu == 0.11270441406241
    @test original_vf052_opt.n5_q == 0.1084034224823743
    @test original_vf052_opt.n6_q == 0.0725442138872437
    @test original_vf052_opt.n5_factor == 0.0381956042627588
    @test original_vf052_opt.n6_factor == 0.0341222647875728
    @test original_vf052_opt.n5_delta_s == 1.5077691100575907
    @test original_vf052_opt.n6_delta_o == 2.9758958354340583
    @test isempty(original_vf052_opt.guard_mus)
    @test original_vf052_opt.max_total_mus == 16
    @test original_vf052_opt.mu_min == 0.07 && original_vf052_opt.mu_max == 0.16
    @test original_vf052.cache_id != original.cache_id
    @test original_vf052.cache_id != original_uf165.cache_id
    original_vf052_config = TOML.parsefile(original_vf052_path)
    @test original_vf052_config["fast_ed"]["run_name"] ==
          "n7_original_vf0_052_auto_scout"
    @test !original_vf052_config["fast_ed"]["allow_cache_release"]
    @test original_vf052_config["pipeline"]["retire_before_start"] == [
        "config/fast_ed/n7_original_uf0_165_cache_retirement.toml",
    ]
    @test !original_vf052_config["pipeline"]["auto_release"]

    mktempdir() do directory
        missing_guard = deepcopy(final_config)
        filter!(mu -> mu != final_opt.guard_mus[1], missing_guard["fast_ed"]["mus"])
        path = joinpath(directory, "missing_guard.toml")
        open(path, "w") do io
            TOML.print(io, missing_guard; sorted=true)
        end
        @test_throws ArgumentError FastEDPipeline.pipeline_options(FastED.load_spec(path))
    end

    pending = Dict{String,Any}(
        "action" => "bundle", "complete" => false, "accepted" => true,
        "stage" => "followup", "next_mus" => [0.145125],
        "next_profile" => "round_2_followup.toml", "next_tasks" => 4,
        "bundle_list" => "bundle_file_list.txt",
    )
    packaged = FastEDPipeline.completed_state(pending, "/tmp/final.tar.gz")
    @test packaged["action"] == "complete" && packaged["complete"]
    @test packaged["accepted"] && packaged["stage"] == "final"
    @test all(!haskey(packaged, key) for key in
              ("next_mus", "next_profile", "next_tasks", "bundle_list"))
    @test pending["action"] == "bundle" && !pending["complete"]

    make_row(mu, q; valid=true, factor=0.03, gaps=Float64[0.04, 0.06, 0.09, 0.09, 0.09]) =
        (mu=Float64(mu), q=Float64(q), valid=valid, factor=Float64(factor),
         raw_gaps=copy(gaps), delta_s=1.5, delta_o=2.9)
    qcurve(mu; center=0.1455, q0=0.09) = sqrt(q0^2 + 400*(mu-center)^2)

    guarded_scout = [make_row(mu, mu == final_opt.guard_mus[1] ? 0.24 :
                                   qcurve(mu; center=final_damped, q0=0.10))
                     for mu in final_vf0.mus]
    guarded_refine = FastEDPipeline.decide_next(
        guarded_scout, final_vf0.mus, "scout", 0, final_opt,
    )
    @test guarded_refine.action == "solve" && guarded_refine.kind == "refine"
    @test all(abs(mu-final_damped) < 0.002 for mu in guarded_refine.mus)

    guard_wins = [make_row(mu, mu == final_opt.guard_mus[1] ? 0.05 : 0.20)
                  for mu in final_vf0.mus]
    guard_extension = FastEDPipeline.decide_next(
        guard_wins, final_vf0.mus, "scout", 0, final_opt,
    )
    @test guard_extension.action == "solve" && guard_extension.kind == "scout_extension"
    @test all(mu < final_opt.guard_mus[1] for mu in guard_extension.mus)

    scout_mus = pilot.mus
    scout = [make_row(mu, qcurve(mu)) for mu in scout_mus]
    refine_decision = FastEDPipeline.decide_next(scout, scout_mus, "scout", 0, opt)
    @test refine_decision.action == "solve" && refine_decision.kind == "refine"
    @test length(refine_decision.mus) == 7
    @test all(isapprox.(diff(refine_decision.mus), 0.00025; atol=1e-14, rtol=0))

    refined = vcat(scout, [make_row(mu, qcurve(mu)) for mu in refine_decision.mus])
    accepted = FastEDPipeline.decide_next(
        refined, refine_decision.mus, "refine", 1, opt,
    )
    @test accepted.action == "accept"
    @test accepted.best.mu in refine_decision.mus
    @test accepted.reason == "measured_minimum_bracketed"

    edge_mus = collect(0.144:0.00025:0.1455)
    edge_rows = [make_row(mu, sqrt(0.09^2 + 100*(mu-0.1458)^2)) for mu in edge_mus]
    followup = FastEDPipeline.decide_next(edge_rows, edge_mus, "refine", 1, opt)
    @test followup.action == "solve" && followup.kind == "followup"
    @test 1 <= length(followup.mus) <= 2
    @test all(mu -> mu > maximum(edge_mus), followup.mus)

    invalid = copy(scout)
    invalid[3] = make_row(invalid[3].mu, invalid[3].q; valid=false)
    @test FastEDPipeline.decide_next(invalid, scout_mus, "scout", 0, opt).action == "review"
    @test FastEDPipeline.decide_next(edge_rows, edge_mus, "followup", 3, opt).reason ==
          "adaptive_round_limit"

    # Regression for the Vf0=.52 review: the fitted vertex lay between valid,
    # rising neighbors, but the configured snap mapped it onto the sampled
    # minimum.  Preserve the strict acceptance threshold and solve the actual
    # fitted vertex instead of stopping with no_new_followup_points.
    stalled_mus = [0.112375, 0.1125, 0.112517687151475]
    stalled_curve(mu) = sqrt(0.0534^2 + 1000*(mu-0.11245)^2)
    stalled_rows = [make_row(mu, stalled_curve(mu)) for mu in stalled_mus]
    stalled = FastEDPipeline.decide_next(stalled_rows, [0.1125], "followup", 2,
                                         original_vf052_opt)
    @test stalled.action == "solve" && stalled.kind == "followup"
    @test stalled.reason == "sample_unsnapped_fit_vertex"
    @test length(stalled.mus) == 1
    @test stalled.mus[1] ≈ 0.11245 atol=1e-13 rtol=0

    mktempdir() do directory
        config = TOML.parsefile(pilot_path)
        config["pipeline"]["name"] = "state_test"
        config["pipeline"]["state_root"] = joinpath(directory, "states")
        config["pipeline"]["archive_root"] = joinpath(directory, "archives")
        config["fast_ed"]["cache_root"] = joinpath(directory, "cache")
        config["fast_ed"]["result_root"] = joinpath(directory, "runs")
        path = joinpath(directory, "pilot.toml")
        open(path, "w") do io
            TOML.print(io, config; sorted=true)
        end
        state = FastEDPipeline.initialize(path)
        @test state["action"] == "ready" && state["stage"] == "scout"
        @test isfile(FastEDPipeline.state_path(FastEDPipeline.pipeline_options(FastED.load_spec(path))))
        @test FastEDPipeline.claim_launch(path)["action"] == "launching"
        @test_throws ArgumentError FastEDPipeline.claim_launch(path)
        @test FastEDPipeline.reset_launch(path)["action"] == "ready"
        FastEDPipeline.claim_launch(path)
        submitted = FastEDPipeline.record_submission(
            path; stage="scout", prepare_job="101", solve_job="102",
            controller_job="103",
        )
        @test submitted["action"] == "awaiting_results"
        @test length(submitted["submission_log"]) == 1
    end

    mktempdir(root) do directory
        config = TOML.parsefile(pilot_path)
        config["pipeline"]["name"] = "bundle_state_test"
        config["pipeline"]["state_root"] = joinpath(directory, "states")
        config["pipeline"]["archive_root"] = joinpath(directory, "archives")
        config["fast_ed"]["cache_root"] = joinpath(directory, "cache")
        config["fast_ed"]["result_root"] = joinpath(directory, "runs")
        path = joinpath(directory, "pilot.toml")
        open(path, "w") do io
            TOML.print(io, config; sorted=true)
        end
        spec = FastED.load_spec(path)
        local_opt = FastEDPipeline.pipeline_options(spec)
        mkpath(spec.cache_directory)
        MottJainED.atomic_toml(FastED.cache_manifest_path(spec), Dict(
            "cache_id" => spec.cache_id, "complete" => true,
        ))
        spectrum = joinpath(directory, "measured_spectrum.csv")
        write(spectrum, "energy\n0.0\n")
        row(mu, q) = (
            cache_id=spec.cache_id, settings_id=spec.settings_id, nm1=7,
            mu=mu, k=20, score_valid=true, definition="custom",
            score_terms="ds_s,j,curlj,dj_rank1,t_rank1", metric="q",
            objective=q, q=q, cost=q^2, factor=0.03, delta_s=1.5,
            delta_o=2.9, state_count=80, reason="", valid=true,
            raw_gaps=[0.03, 0.06, 0.09, 0.09, 0.09], profile=path,
        )
        rows = [row(0.145, 0.08), row(0.145125, 0.07), row(0.14525, 0.08)]
        relation = (cache_id=spec.cache_id, settings_id=spec.settings_id,
                    nm1=7, mu=0.145125, term="ds_s", label="dS-S",
                    raw_gap=0.03, target_gap=1.0, scaled_gap=1.0, residual=0.0)
        selected = (mu=0.145125, label="S", l2=0, c2=0, raw_rank=2,
                    z=1, r=1, sector_rank=2, energy=0.03, gap=0.03, dimension=1.0)
        analysis = (
            base=spec, rows=rows, relations=[relation], selected=[selected],
            files=[(file=relpath(path, FastED.PROJECT_ROOT), sha256="abc", bytes=1)],
            max_quantum_error=0.0, max_copy_split=0.0, max_selected_rank=2,
            issues=String[],
            spectrum_paths=Dict(FastEDPipeline.key(0.145125) => spectrum),
        )
        decision = (
            best=rows[2], triplet=(xs=[0.145, 0.145125, 0.14525],
                                   qs=[0.08, 0.07, 0.08], vertex=0.145125,
                                   predicted_q=0.07, improvement=0.0),
            factor_jump=0.0, gap_jump=0.0,
        )
        state = FastEDPipeline.initialize(path)
        state["action"] = "awaiting_results"
        state["stage"] = "followup"
        state["next_mus"] = [0.145125]
        state["next_profile"] = "round_2_followup.toml"
        state["next_tasks"] = 4
        FastEDPipeline.finalise!(state, local_opt, analysis, decision)
        live = TOML.parsefile(FastEDPipeline.state_path(local_opt))
        snapshot_path = joinpath(FastEDPipeline.final_directory(local_opt), "pipeline_state.toml")
        snapshot = TOML.parsefile(snapshot_path)
        bundle_paths = readlines(live["bundle_list"])
        @test live["action"] == "bundle" && !live["complete"] && live["stage"] == "final"
        @test snapshot["action"] == "complete" && snapshot["complete"]
        @test snapshot["stage"] == "final" && snapshot["packaged_snapshot"]
        @test !haskey(snapshot, "next_mus") && !haskey(snapshot, "bundle_list")
        @test relpath(FastEDPipeline.state_path(local_opt), FastED.PROJECT_ROOT) ∉ bundle_paths
        @test relpath(FastEDPipeline.final_directory(local_opt), FastED.PROJECT_ROOT) ∈ bundle_paths
        archive = live["bundle_path"]
        mkpath(dirname(archive))
        write(archive, "verified archive bytes")
        digest = bytes2hex(sha256(read(archive)))
        finished = FastEDPipeline.mark_bundled(path, archive, digest)
        @test finished["action"] == "complete" && finished["complete"]
        @test finished["stage"] == "final" && finished["bundle_sha256"] == digest
        @test finished["bundle_bytes"] == filesize(archive)
        @test all(!haskey(finished, key) for key in
                  ("next_mus", "next_profile", "next_tasks", "bundle_list"))
    end

    launcher = read(joinpath(root, "scripts", "submit_fast_ed_pipeline.sh"), String)
    bootstrap = read(joinpath(root, "slurm", "fast_ed_pipeline_bootstrap.sbatch"), String)
    controller = read(joinpath(root, "slurm", "fast_ed_pipeline_controller.sbatch"), String)
    recovery = read(joinpath(root, "slurm", "fast_ed_pipeline_review_recovery.sbatch"), String)
    cli = read(joinpath(root, "scripts", "fast_ed_pipeline.jl"), String)
    @test occursin(raw"%${max_concurrent}", launcher)
    @test occursin(raw"afterany:$solve_job", launcher)
    @test occursin("fast_ed_pipeline_bootstrap.sbatch", launcher)
    @test occursin("retire-configured-caches", bootstrap)
    @test occursin("JULIA_NUM_PRECOMPILE_TASKS=1", bootstrap)
    @test occursin("retire-configured-caches", cli)
    @test occursin(raw"%${max_concurrent}", controller)
    @test occursin(raw"afterany:$solve_job", controller)
    @test occursin("resume-stalled-fit", recovery)
    @test occursin("fast_ed_pipeline_controller.sbatch", recovery)
    @test occursin("resume-stalled-fit", cli)
    @test !occursin("sleep ", launcher*bootstrap*controller*recovery)
end

@testset "Fresh N7 five-point scout uses the validated array workflow" begin
    config_root = joinpath(@__DIR__, "..", "config", "fast_ed")
    scout = FastED.load_spec(joinpath(config_root, "n7_uf0_165_k20_scout_restart.toml"))
    retired = FastED.load_spec(joinpath(config_root, "n7_uf0_165_search.toml"))
    retained = FastED.load_spec(joinpath(config_root, "n7_retained_k20.toml"))
    center = 2*0.14320460558017-0.14128848931940
    @test scout.mus ≈ center .+ [-0.01, -0.005, 0.0, 0.005, 0.01] atol=1e-14
    @test FastED.plan(scout).tasks == 20
    @test scout.nm1 == 7 && scout.solver.k == 20 && !scout.solver.warm_start
    @test scout.cache_id == retired.cache_id && scout.settings_id == retired.settings_id
    @test scout.result_directory != retired.result_directory
    @test scout.cache_id != retained.cache_id
    @test !haskey(scout.config, "mu_search")
    @test !scout.config["fast_ed"]["allow_cache_release"]
    @test scout.config["scout"]["requires_review"] && !scout.config["scout"]["muc_confirmed"]
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

@testset "Uf0=1.65 refinement preserves the completed scout cache and solver" begin
    root = joinpath(@__DIR__, "..", "config", "fast_ed")
    scout = FastED.load_spec(joinpath(root, "n7_uf0_165_k20_scout_restart.toml"))
    refine = FastED.load_spec(joinpath(root, "n7_uf0_165_k20_refine.toml"))
    @test refine.nm1 == 7 && refine.solver.k == 20
    @test refine.cache_id == scout.cache_id
    @test refine.cache_directory == scout.cache_directory
    @test refine.settings_id == scout.settings_id
    @test FastED.solver_dict(refine.solver) == FastED.solver_dict(scout.solver)
    @test FastED.coupling_dict(refine.couplings; include_mu=false) ==
          FastED.coupling_dict(scout.couplings; include_mu=false)
    @test refine.terms == scout.terms && refine.metric == scout.metric
    @test refine.run_name == "n7_uf0_165_refine" && refine.result_directory != scout.result_directory
    @test refine.mus == [0.14425, 0.14450, 0.14475, 0.14500, 0.14525, 0.14550, 0.14575]
    @test all(isapprox.(diff(refine.mus), 0.00025; atol=1e-14, rtol=0))
    @test FastED.plan(refine).tasks == 28 && isempty(intersect(refine.mus, scout.mus))
    @test all(first(refine.mus) < mu < last(refine.mus) for mu in
              [0.14493829304945727, 0.1448838386844746, 0.14512072184094])
    @test !haskey(refine.config, "mu_search") && !refine.config["fast_ed"]["allow_cache_release"]
    @test refine.config["refinement"]["source_array_job"] == "548422" &&
          refine.config["refinement"]["requires_review"] && !refine.config["refinement"]["muc_confirmed"]
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

@testset "Original-audit local FSS stops at N=6" begin
    root = joinpath(@__DIR__, "..")
    profile_root = joinpath(root, "config", "two_size_tuning")
    reference = TOML.parsefile(joinpath(profile_root, "s_stage12_vf_retained.toml"))
    search_keys = filter(key -> startswith(key, "optimize_") || startswith(key, "mu_") ||
        startswith(key, "score_") || key == "k", keys(reference["two_size_tuning"]))
    fixed = Dict(
        "Uf" => 0.46, "U0" => 4.14, "Vf" => 0.0, "V0" => 0.525, "t" => 0.5,
    )
    cases = (
        ("n56_original_uf0_165.toml", "Uf0", 1.65, 1.65, 0.41),
        ("n56_original_uf0_200.toml", "Uf0", 2.00, 2.00, 0.41),
        ("n56_original_vf0_030.toml", "Vf0", 0.30, 1.834, 0.30),
        ("n56_original_vf0_052.toml", "Vf0", 0.52, 1.834, 0.52),
    )
    points = Set{Tuple{Float64,Float64}}()
    for (name, parameter, value, uf0, vf0) in cases
        config = TOML.parsefile(joinpath(profile_root, name))
        hamiltonian = config["hamiltonian"]
        tuning = config["two_size_tuning"]
        @test tuning["guide_nm_values"] == [3, 4]
        @test tuning["match_nm_values"] == [5, 6]
        @test maximum(vcat(tuning["guide_nm_values"], tuning["match_nm_values"])) == 6
        @test tuning["scan_parameter"] == parameter
        @test tuning["scan_values"] == [value]
        @test hamiltonian["Uf0"] == uf0 && hamiltonian["Vf0"] == vf0
        for (key, expected) in fixed
            @test hamiltonian[key] == expected
        end
        for key in search_keys
            @test tuning[key] == reference["two_size_tuning"][key]
        end
        push!(points, (uf0, vf0))
    end
    @test length(points) == 4
    @test !((1.834, 0.41) in points)

    submission = read(joinpath(root, "slurm", "fss_original_local_n56.sbatch"), String)
    @test occursin("#SBATCH --array=0-3%4", submission)
    @test count(name -> occursin(name, submission), first.(cases)) == 4
    @test !occursin("scripts/fast_ed.jl", submission)
end
