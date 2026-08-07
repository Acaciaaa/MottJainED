using Test
using MottJainED
using DataFrames

@testset "Portable FuzzifiED dependency" begin
    root = dirname(@__DIR__)
    project = MottJainED.TOML.parsefile(joinpath(root, "Project.toml"))
    source = project["sources"]["FuzzifiED"]
    @test source["url"] == "https://github.com/FuzzifiED/FuzzifiED.jl.git"
    @test source["rev"] == "29a0cc9e06bcb5b30d3cf9f6db6416917f8a573f"
    @test !haskey(source, "path")

    manifest = MottJainED.TOML.parsefile(joinpath(root, "Manifest.toml"))
    fuzzified = only(manifest["deps"]["FuzzifiED"])
    @test fuzzified["git-tree-sha1"] == "fe4f9de48a7b76014281b87a385088dea0733aac"
    @test fuzzified["repo-rev"] == source["rev"]
    @test !haskey(fuzzified, "path")
end

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

@testset "Quantum-number labels" begin
    boundary = state(1.0, 20.0 + 3e-14, 3.0 - 2e-14)
    rejected = state(2.0, 20.2, 3.0)
    @test MottJainED._state_quantum_labels(boundary, 2e-3) == (20, 3)
    @test MottJainED._state_quantum_labels(rejected, 2e-3) === nothing
    selected = MottJainED._states_in_sector(
        [rejected, boundary], 20, 3; quantum_tol=2e-3,
    )
    @test selected == [boundary]

    table = MottJainED.spectrum_dataframe(
        [boundary, rejected]; mu=0.1, nm1=6, quantum_tol=2e-3,
    )
    @test table.l2[1] == 20
    @test table.c2[1] == 3
    @test ismissing(table.l2[2])
    @test table.l2_raw[1] == boundary.l2
end

@testset "Legacy raw-rank selection" begin
    states = SpectrumState[
        state(0.0, 1e-14, -1e-14),
        state(0.75, 2.0 + 1e-14, 3.0 - 1e-14),
        state(1.25, 2e-14, 1e-14),
        state(2.0, 2.0, 0.0),
    ]
    @test MottJainED._legacy_singlet_gap(states, 2e-3) ≈ 1.25
    @test MottJainED._legacy_j_gap(states, 2e-3) ≈ 0.75
end

@testset "Gap plot output" begin
    data = DataFrame(
        status=["ok", "ok"], nm1=[4, 4], mu=[0.0, 0.1],
        scaled_scalar_gap=[0.2, 0.3], scaled_j_gap=[0.4, 0.5],
    )
    mktempdir() do directory
        path = joinpath(directory, "gap.png")
        figure = MottJainED._plot_gap_curves(
            data, :scaled_j_gap, "test gap", "test", path,
        )
        @test figure !== nothing
        @test isfile(path)
    end
end

@testset "Density ground state and plot" begin
    model = build_model(nm1=2)
    settings = SolverSettings(k=3)
    cache = prepare_spectrum(model, Couplings(), settings)
    ground = MottJainED.solve_ground_state(cache, 0.05)
    reference = first(solve_spectrum(cache, 0.05; keep_vectors=true))
    @test ground.energy ≈ reference.energy
    @test ground.vector !== nothing

    data = DataFrame(
        status=["ok", "ok"], mu=[0.0, 0.1], nf_per_orbital=[1.2, 1.3],
    )
    mktempdir() do directory
        path = joinpath(directory, "density.png")
        figure = MottJainED._plot_density_curve(data, path)
        @test figure !== nothing
        @test isfile(path)
    end
end

@testset "CFT tower classification" begin
    states = SpectrumState[
        state(0.0, 0, 0), state(1.2, 0, 0), state(3.2, 0, 0),
        state(2.2, 2, 0), state(2.7, 6, 0), state(3.0, 6, 0),
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
    @test score.labels == ["dS-S", "J", "curlJ", "dJ(rank1)", "T(rank2)"]
    @test score.raw_gaps ≈ [1.0, 2.0, 3.0, 3.0, 3.0]
    @test score.delta_s ≈ 1.2
    @test score.delta_o ≈ 1.6
end

@testset "Selectable legacy CFT scores" begin
    states = SpectrumState[
        state(0.0, 0, 0), state(1.2, 0, 0), state(3.2, 0, 0),
        state(2.2, 2, 0), state(3.0, 6, 0), state(3.2, 6, 0),
        state(2.0, 2, 3; z=1), state(2.0, 2, 3; z=-1), state(3.0, 2, 3),
        state(3.0, 6, 3; z=1), state(3.0, 6, 3; z=-1), state(3.0, 6, 3),
        state(2.0, 2, 6), state(3.0, 6, 6), state(1.6, 0, 3),
    ]
    seven = cft_score(states; definition=:fss7)
    @test seven.valid
    @test seven.definition == :fss7
    @test seven.metric == :cost
    @test seven.objective ≈ 0.0 atol=1e-12
    @test seven.raw_gaps ≈ [1, 1, 2, 2, 3, 3, 3]

    eight = cft_score(states; definition=:optimization8)
    @test eight.valid
    @test eight.metric == :cost
    @test eight.objective ≈ 0.0 atol=1e-12
    @test eight.raw_gaps ≈ [1, 1, 1, 2, 2, 3, 3, 3]

    seven_q = cft_score(states; definition=:fss7, metric=:q)
    @test seven_q.metric == :q
    @test seven_q.objective == seven_q.q

    custom = cft_score(states; terms=[:j, :t_rank1], metric=:q)
    @test custom.valid
    @test custom.definition == :custom
    @test custom.terms == [:j, :t_rank1]
    @test custom.raw_gaps ≈ [2, 3]
    @test custom.target_gaps ≈ [2, 3]

    # 未选择的关系不应要求其 sector/rank 存在。
    minimal = SpectrumState[state(0.0, 0, 0), state(2.0, 2, 3)]
    only_j = cft_score(minimal; terms=["J"], metric=:q)
    @test only_j.valid
    @test only_j.terms == [:j]
end

@testset "Configuration" begin
    config = load_config(joinpath(dirname(@__DIR__), "config", "default.toml"))
    @test config["model"]["nm1"] == 5
    @test config["density"]["k"] == 3
    @test config["critical"]["k"] == 10
    @test config["critical"]["mu_count"] == 9
    @test !haskey(config["critical"], "coarse_points")
    @test config["fss"]["mu_count"] == 9
    @test !haskey(config["critical"], "score_terms")
    @test !haskey(config["fss"], "score_terms")
    @test config["fss"]["methods"] == ["grid", "optimize"]
    @test !haskey(config["optimization"], "score_terms")
    @test !haskey(config["optimization"], "free")
    @test !haskey(config["optimization"], "values")
    @test !haskey(config["optimization"], "bounds")
    @test config["optimization"]["u0_over_uf"] == 9.0
    @test config["optimization"]["algorithm"] == "auto"
    @test config["generator"]["config_root"] == "config/generator"
    @test !haskey(config["generator"], "data_root")
    fit_config = MottJainED.TOML.parsefile(joinpath(
        dirname(@__DIR__), "config", "generator", "templates", "generator_fit.toml",
    ))
    tower_config = MottJainED.TOML.parsefile(joinpath(
        dirname(@__DIR__), "config", "generator", "templates", "tower.toml",
    ))
    @test fit_config["fit"]["source"] == "S"
    @test fit_config["states"]["dS"]["rank"] == 1
    @test !haskey(tower_config, "fit")
    critical_case = load_config(
        joinpath(dirname(@__DIR__), "config", "default.toml");
        override=joinpath(dirname(@__DIR__), "config", "critical_profiles", "critical5.toml"),
    )
    fss_case = load_config(
        joinpath(dirname(@__DIR__), "config", "default.toml");
        override=joinpath(dirname(@__DIR__), "config", "fss_profiles", "fss7.toml"),
    )
    fss5_case = load_config(
        joinpath(dirname(@__DIR__), "config", "default.toml");
        override=joinpath(dirname(@__DIR__), "config", "fss_profiles", "fss5.toml"),
    )
    @test length(critical_case["critical"]["score_terms"]) == 5
    @test critical_case["output"]["run_name"] == "critical5"
    @test length(fss_case["fss"]["score_terms"]) == 7
    @test fss_case["fss"]["k"] == 30
    @test fss_case["output"]["run_name"] == "fss7"
    @test fss5_case["fss"]["score_terms"] ==
          ["ds_s", "j", "curlj", "dj_rank1", "t_rank1"]
    @test fss5_case["fss"]["score_metric"] == "q"
    @test fss5_case["fss"]["k"] == 10
    @test fss5_case["output"]["run_name"] == "fss5"
    mktempdir() do directory
        first_case = MottJainED._ordinary_task_output(
            critical_case, :critical, directory,
        )
        @test basename(first_case) == "critical5_01"
        unrelated_change = deepcopy(critical_case)
        unrelated_change["gap"]["mu_count"] = 99
        @test MottJainED._ordinary_task_output(
            unrelated_change, :critical, directory,
        ) == first_case
        relevant_change = deepcopy(critical_case)
        relevant_change["critical"]["mu_count"] = 11
        @test basename(MottJainED._ordinary_task_output(
            relevant_change, :critical, directory,
        )) == "critical5_02"
    end
    optimization_couplings = MottJainED._optimization_couplings(
        config, config["optimization"],
    )
    @test optimization_couplings.Uf == 0.5
    @test optimization_couplings.mu == 0.05
    @test MottJainED._fss_methods(config["fss"], Dict{String,String}()) == [:grid, :optimize]
    @test MottJainED._fss_methods(config["fss"], Dict("method" => "grid")) == [:grid]
    @test MottJainED._fss_methods(config["fss"], Dict("method" => "optimize")) == [:optimize]
    @test MottJainED._fss_methods(config["fss"], Dict("method" => "both")) == [:grid, :optimize]

    profile_path = joinpath(
        dirname(@__DIR__), "config", "optimization_profiles", "mu_uf0_v0.toml",
    )
    profile = load_config(
        joinpath(dirname(@__DIR__), "config", "default.toml"); override=profile_path,
    )
    @test profile["optimization"]["free"] == ["Uf0", "V0", "mu"]
    @test length(profile["optimization"]["score_terms"]) == 8
    @test profile["optimization"]["values"]["Vf"] == 0.0
    @test profile["optimization"]["values"]["Vf0"] == 0.4
    MottJainED._apply_cli_config_overrides!(
        profile, "plan", Dict("nm1" => "7", "k" => "12"), profile_path,
    )
    @test profile["model"]["nm1"] == 7
    @test profile["optimization"]["k"] == 12
    @test !haskey(profile["output"], "run_name")

    mktempdir() do directory
        first_output = MottJainED._optimization_output(
            profile, profile_path, directory,
        )
        @test basename(first_output) == "mu_uf0_v0_nm7_k12_01"
        @test MottJainED._optimization_output(
            profile, profile_path, directory,
        ) == first_output
        changed_profile = deepcopy(profile)
        changed_profile["optimization"]["values"]["Uf0"] = 4.0
        second_output = MottJainED._optimization_output(
            changed_profile, profile_path, directory,
        )
        @test basename(second_output) == "mu_uf0_v0_nm7_k12_02"
    end
    @test endswith(
        MottJainED._output_root(config, :spectrum),
        joinpath("output", "spectrum"),
    )

    expected_profiles = Dict(
        "mu_only.toml" => ["mu"],
        "mu_uf0.toml" => ["Uf0", "mu"],
        "mu_uf0_v0.toml" => ["Uf0", "V0", "mu"],
        "mu_uf_uf0_vf0.toml" => ["Uf", "Uf0", "Vf0", "mu"],
        "mu_uf_uf0_vf0_v0.toml" => ["Uf", "Uf0", "Vf0", "V0", "mu"],
    )
    for (filename, expected_free) in expected_profiles
        candidate = load_config(
            joinpath(dirname(@__DIR__), "config", "default.toml");
            override=joinpath(dirname(profile_path), filename),
        )
        @test candidate["optimization"]["free"] == expected_free
        @test candidate["optimization"]["values"]["Vf"] == 0.0
        @test haskey(candidate["optimization"], "score_terms")
    end
    mktempdir() do directory
        resolved_path = MottJainED.write_resolved_config(
            directory, profile; base_config="base.toml", override_config=profile_path,
        )
        resolved = MottJainED.TOML.parsefile(resolved_path)
        @test resolved["model"]["nm1"] == 7
        @test resolved["optimization"]["k"] == 12
        @test endswith(resolved["resolved_sources"]["override_config"], "mu_uf0_v0.toml")
    end
    mktempdir() do directory
        state = MottJainED.start_task_logging(directory)
        @info "file-only test progress" point=3
        MottJainED.stop_task_logging(state)
        contents = read(joinpath(directory, "run.log"), String)
        @test occursin("file-only test progress", contents)
        @test occursin("point = 3", contents)
    end
end

@testset "Generator point registry and ED snapshot identity" begin
    mktempdir() do directory
        registry = joinpath(directory, "generator_points.csv")
        MottJainED.ensure_generator_registry(registry)
        best = joinpath(directory, "best.csv")
        MottJainED.atomic_csv(best, DataFrame([(
            nm1=6, Uf=0.8, Uf0=3.1, U0=7.2, Vf=0.0, Vf0=0.6,
            V0=1.0, t=0.5, mu_initial=0.22, factor=0.04,
            objective=0.1, q=0.1, cost=0.2, delta_s=1.3, delta_o=0.6,
            score_definition="custom-q-a,b",
        )]))
        registered = MottJainED.register_optimization_point(
            registry, "nm6_good_01", best; notes="checked",
        )
        @test registered.point_id == "nm6_good_01"
        @test registered.nm1 == 6
        @test registered.couplings.mu == 0.22
        @test registered.couplings.V0 == 1.0
        @test registered.factor == 0.04
        @test registered.metadata["notes"] == "checked"
        @test_throws ArgumentError MottJainED.register_optimization_point(
            registry, "nm6_good_01", best,
        )

        settings = SolverSettings(k=40)
        provenance = Dict(
            "package_source_signature" => "package-test",
            "fuzzified_source_signature" => "fuzzified-test",
        )
        id1 = MottJainED.generator_snapshot_id(
            registered, settings; include_adjoint=true, adjoint_k=30,
            provenance=provenance,
        )
        id2 = MottJainED.generator_snapshot_id(
            registered, settings; include_adjoint=true, adjoint_k=31,
            provenance=provenance,
        )
        @test id1 != id2
        changed = GeneratorPoint(
            point_id=registered.point_id, nm1=registered.nm1,
            couplings=MottJainED.with_coupling(registered.couplings, :V0, 1.1),
            factor=registered.factor, metadata=registered.metadata,
        )
        @test id1 != MottJainED.generator_snapshot_id(
            changed, settings; include_adjoint=true, adjoint_k=30,
            provenance=provenance,
        )
        @test MottJainED.generator_snapshot_directory(
            directory, registered, id1,
        ) == joinpath(directory, registered.point_id)

        case_directory = MottJainED.ensure_generator_case_config(
            registered.point_id,
            joinpath(directory, "config", "generator"),
            joinpath(dirname(@__DIR__), "config", "generator", "templates"),
        )
        @test isfile(joinpath(case_directory, "generator_fit.toml"))
        @test isfile(joinpath(case_directory, "tower.toml"))

        first_tower = MottJainED._tower_output_directory(directory, "tower-signature-a")
        @test first_tower.analysis_id == "tower_01"
        @test isfile(joinpath(first_tower.directory, ".tower_identity.toml"))
        @test MottJainED._tower_output_directory(
            directory, "tower-signature-a",
        ).directory == first_tower.directory
        second_tower = MottJainED._tower_output_directory(directory, "tower-signature-b")
        @test second_tower.analysis_id == "tower_02"

        raw_states = SpectrumState[
            SpectrumState(0.0, 0.0, 0.0, SectorKey(1, 1), 1, [1.0, 0.0], :basis),
            SpectrumState(1.0, 2.0, 0.0, SectorKey(1, 1), 2, [0.0, 1.0], :basis),
        ]
        packed = MottJainED._pack_spectrum(raw_states, :standard)
        snapshot = GeneratorEDSnapshot(
            snapshot_id=id1, created_at="test", point=registered,
            settings=settings, include_adjoint=false, adjoint_k=30,
            sectors=packed, model_summary=Dict("nm1" => 6),
            provenance=provenance,
        )
        path = joinpath(directory, "snapshot.jld2")
        MottJainED.atomic_jldsave(path; snapshot=snapshot)
        loaded = MottJainED.load_generator_snapshot(path)
        restored = snapshot_states(loaded)
        @test length(restored) == 2
        @test restored[2].vector == [0.0, 1.0]
        @test restored[2].basis == :basis

        # 旧快照曾把 Project.toml 的路径文字算进 ED 身份。只要快照旧签名自洽、
        # Hamiltonian/solver/FuzzifiED 都相同，就迁移一次；之后 sidecar 严格锁定新身份。
        legacy_provenance = Dict{String,Any}(
            "ed_source_signature" => "legacy-project-path-hash",
            "fuzzified_source_signature" => "same-fuzzified-source",
        )
        legacy_signature = MottJainED.generator_snapshot_id(
            registered, settings; include_adjoint=false, adjoint_k=30,
            provenance=legacy_provenance,
        )
        legacy_provenance["ed_identity_signature"] = legacy_signature
        legacy_snapshot = GeneratorEDSnapshot(
            snapshot_id=registered.point_id, created_at="test", point=registered,
            settings=settings, include_adjoint=false, adjoint_k=30,
            sectors=packed, model_summary=Dict("nm1" => 6),
            provenance=legacy_provenance,
        )
        current_provenance = Dict{String,Any}(
            "ed_source_signature" => "physics-only-hash",
            "fuzzified_source_signature" => "same-fuzzified-source",
        )
        current_signature = MottJainED.generator_snapshot_id(
            registered, settings; include_adjoint=false, adjoint_k=30,
            provenance=current_provenance,
        )
        migration_directory = joinpath(directory, "legacy_migration")
        mkpath(migration_directory)
        @test MottJainED._validate_snapshot_identity(
            legacy_snapshot, registered, settings, current_provenance,
            current_signature, migration_directory;
            include_adjoint=false, adjoint_k=30,
        )
        @test isfile(joinpath(migration_directory, ".ed_identity.toml"))
        changed_provenance = copy(current_provenance)
        changed_provenance["ed_source_signature"] = "changed-physics-code"
        changed_signature = MottJainED.generator_snapshot_id(
            registered, settings; include_adjoint=false, adjoint_k=30,
            provenance=changed_provenance,
        )
        @test !MottJainED._validate_snapshot_identity(
            legacy_snapshot, registered, settings, changed_provenance,
            changed_signature, migration_directory;
            include_adjoint=false, adjoint_k=30,
        )

        tower_table = DataFrame(
            status=["ok", "ok", "skipped"],
            relation=["curlJ_same_L", "curlJ_same_L", "optional"],
            mode=["same_angular", "same_angular", "regular"],
            input=["curlJ_other", "curlJ_other", "missing"],
            target=["J_other", "boxJ_other", "unknown"],
            target_l=Union{Missing,Int}[1, 1, missing],
            overlap=Union{Missing,Float64}[0.7, 0.2, missing],
            total_overlap=Union{Missing,Float64}[0.9, 0.9, missing],
            delta_energy=Union{Missing,Float64}[-1.0, 1.0, missing],
            scaled_delta=Union{Missing,Float64}[-0.5, 0.5, missing],
            message=["", "", "target not found"],
        )
        io = IOBuffer()
        MottJainED._print_tower_summary(
            tower_table; point_id="nm6_good_01", analysis_id="tower_01",
            state_names=Dict(
                "curlJ_other" => "curl J", "J_other" => "J",
                "boxJ_other" => "box J",
            ), io=io,
        )
        printed = String(take!(io))
        @test occursin("| Input  | l' | Target 1 Ovlp", printed)
        @test occursin("| curl J | 1  | J(-0.50)", printed)
        @test occursin("box J(0.50)", printed)
        @test occursin("0.9000 |", printed)
        @test occursin("Skipped optional: target not found", printed)
    end
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

    grid_result = scan_mu(cache, [0.02, 0.03])
    @test grid_result.mus == [0.02, 0.03]
    @test grid_result.evaluations == 2
    @test length(grid_result.scores) == 2
    @test grid_result.completed

    family = MottJainED._prepare_linear_family(
        model, couplings, [:Uf], settings; u0_over_uf=9.0,
    )
    tied_states = MottJainED._solve_linear(family, [0.7])
    tied_couplings = MottJainED.with_coupling(
        MottJainED.with_coupling(couplings, :Uf, 0.7), :U0, 6.3,
    )
    direct_cache = prepare_spectrum(model, tied_couplings, settings)
    direct_states = solve_spectrum(direct_cache, tied_couplings.mu)
    @test first(tied_states).energy ≈ first(direct_states).energy atol=1e-11

    fixed_uf_family = MottJainED._prepare_linear_family(
        model, Couplings(Uf=0.5, U0=99.0), [:mu], settings; u0_over_uf=9.0,
    )
    @test fixed_uf_family.base.U0 == 4.5
    @test MottJainED._couplings_from_values(fixed_uf_family, [0.2]).U0 == 4.5
    @test_throws ArgumentError MottJainED._prepare_linear_family(
        model, couplings, [:U0], settings; u0_over_uf=9.0,
    )
end
