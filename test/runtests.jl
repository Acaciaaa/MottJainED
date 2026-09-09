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

    manifest_path = joinpath(
        root, "Manifest-v$(VERSION.major).$(VERSION.minor).toml",
    )
    @test isfile(manifest_path)
    manifest = MottJainED.TOML.parsefile(manifest_path)
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

@testset "Filtered rescaled spectrum" begin
    states = SpectrumState[
        state(0.0, 0, 0),
        state(0.6, 2, 0),
        state(0.8, 6, 0),
        state(0.9, 0, 3),
        state(1.0, 2, 3; z=1),
        state(1.0 + 1e-10, 2, 3; z=-1), # 同一物理能级的离散 sector 副本
        state(1.5, 2, 3),               # 不同能量必须保留为另一行
        state(2.0, 2, 3),
        state(1.2, 6, 3),
        state(0.9, 12, 8),              # 未请求的量子数必须排除
    ]
    result = spectrum_level_table(
        states; l2_values=[0, 2, 6], c2_values=[0, 3],
        levels_per_block=1, factor=0.5,
    )
    table = result.data
    @test names(table) == [
        "L2=0 C2=0", "L2=2 C2=0", "L2=6 C2=0",
        "L2=0 C2=3", "L2=2 C2=3", "L2=6 C2=3",
    ]
    @test nrow(table) == 1
    @test collect(table[1, :]) ≈ [0.0, 1.2, 1.6, 1.8, 2.0, 2.4]
    @test length(result.catalog[(2, 3)]) == 3
    @test result.catalog[(2, 3)][1].multiplicity == 2

    pair = spectrum_level_table(
        states; l2_values=[2], c2_values=[3], levels_per_block=2, factor=0.5,
    ).data
    @test pair[!, "L2=2 C2=3"] ≈ [2.0, 3.0]
    @test_throws ArgumentError spectrum_level_table(
        states; l2_values=[0], c2_values=[0], levels_per_block=2, factor=0.5,
    )
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

@testset "FSS anchor and size-continuation search" begin
    evaluations = NamedTuple[]
    function multiwell(mu)
        left = (mu - 0.25)^2 + 0.08
        narrow_global = 25.0 * (mu - 0.77)^2 + 0.01
        return (valid=true, objective=min(left, narrow_global))
    end
    result = MottJainED._global_grid_refine(
        multiwell; mu_min=0.0, mu_max=1.0, mu_count=21,
        abs_tol=1.0e-7, max_iterations=100,
        on_evaluation=(mu, score, source) -> push!(
            evaluations, (mu=mu, objective=score.objective, source=source),
        ),
    )
    @test result.mu ≈ 0.77 atol=5e-5
    @test result.score.objective ≈ 0.01 atol=1e-8
    @test result.candidate_count == 2
    @test result.refined_count == 2
    @test result.refinements_converged == 2
    @test result.grid_best_mu == 0.75
    @test result.grid_best_objective >= result.score.objective
    @test startswith(result.best_source, "anchor_refine_")
    @test result.evaluations == length(evaluations)
    @test count(row -> row.source == "anchor_grid", evaluations) == 21

    # 大尺寸从前一 size 的结果附近续接；远处 Brent 只是候选，不能覆盖更低的局部支路。
    continuation = MottJainED._continuation_grid_refine(
        multiwell; center=0.76, mu_min=0.0, mu_max=1.0,
        local_half_width=0.04, local_count=9, max_expansions=2,
        abs_tol=1.0e-7, max_iterations=100,
    )
    @test continuation.mu ≈ 0.77 atol=5e-5
    @test continuation.search_mode == "continuation"
    @test continuation.center_mu == 0.76
    @test continuation.wide_brent_converged
    @test continuation.score.objective <= continuation.wide_brent_objective

    # 回归 N=6 型故障：宽 Brent 只看见 0.142 的宽谷，前一 size 的 seed 能抓住
    # 0.096 附近更低但很窄的谷；最后必须保留后者。
    function hidden_narrow_valley(mu)
        if abs(mu-0.096) < 0.006
            return (valid=true, objective=0.29 + 80.0*(mu-0.096)^2)
        end
        return (valid=true, objective=0.57 + (mu-0.142)^2)
    end
    guarded = MottJainED._continuation_grid_refine(
        hidden_narrow_valley; center=0.09584,
        mu_min=-0.095, mu_max=0.305,
        local_half_width=0.02, local_count=9, max_expansions=3,
        abs_tol=1.0e-7, max_iterations=100,
    )
    @test guarded.mu ≈ 0.096 atol=5e-5
    @test guarded.wide_brent_mu ≈ 0.142 atol=5e-5
    @test guarded.score.objective < guarded.wide_brent_objective
    @test startswith(guarded.best_source, "local_")

    # 自适应模式在平稳 continuation 上省掉 wide challenger；审计点仍强制比较。
    adaptive = MottJainED._continuation_grid_refine(
        mu -> (valid=true, objective=(mu-0.51)^2);
        center=0.50, mu_min=0.0, mu_max=1.0,
        local_half_width=0.05, local_count=9, max_expansions=2,
        abs_tol=1.0e-7, max_iterations=100, wide_mode=:adaptive,
    )
    @test adaptive.mu ≈ 0.51 atol=5e-5
    @test !adaptive.wide_brent_ran
    @test isempty(adaptive.wide_guard_reason)
    @test adaptive.completed

    audited = MottJainED._continuation_grid_refine(
        hidden_narrow_valley; center=0.09584,
        mu_min=-0.095, mu_max=0.305,
        local_half_width=0.02, local_count=9, max_expansions=3,
        abs_tol=1.0e-7, max_iterations=100, wide_mode=:adaptive,
        force_wide=true, force_wide_reason="first_point_audit",
    )
    @test audited.wide_brent_ran
    @test audited.wide_guard_reason == "first_point_audit"
    @test audited.wide_disagreement
    @test occursin("muc", audited.wide_disagreement_reason)
    @test occursin("objective", audited.wide_disagreement_reason)

    # 前一 size 的 seed 偏了一点时，窗口最低落在边界会触发扩大，而不是直接接受边界。
    expanded = MottJainED._continuation_grid_refine(
        mu -> (valid=true, objective=(mu-0.56)^2);
        center=0.50, mu_min=0.0, mu_max=1.0,
        local_half_width=0.02, local_count=5, max_expansions=3,
        abs_tol=1.0e-7, max_iterations=100,
    )
    @test expanded.mu ≈ 0.56 atol=5e-5
    @test expanded.expansions >= 2
    @test expanded.search_lower <= 0.56 <= expanded.search_upper

    jump_guarded = MottJainED._continuation_grid_refine(
        mu -> (valid=true, objective=(mu-0.56)^2);
        center=0.50, mu_min=0.0, mu_max=1.0,
        local_half_width=0.02, local_count=5, max_expansions=3,
        abs_tol=1.0e-7, max_iterations=100, wide_mode=:adaptive,
        wide_jump_tol=0.03,
    )
    @test jump_guarded.wide_brent_ran
    @test occursin("large_muc_jump", jump_guarded.wide_guard_reason)

    # 无效区间不能胜出；若最低点确实在总边界，则保留边界而不强行做 Brent。
    endpoint = MottJainED._global_grid_refine(
        mu -> (valid=mu >= 0.2, objective=1.0-mu);
        mu_min=0.0, mu_max=1.0, mu_count=6,
    )
    @test endpoint.mu == 1.0
    @test endpoint.at_boundary
    @test endpoint.invalid == 1
    @test endpoint.refined_count == 0

    @test MottJainED._grid_local_minimum_indices([
        (valid=true, objective=2.0),
        (valid=true, objective=1.0),
        (valid=true, objective=1.0),
        (valid=true, objective=2.0),
    ]) == [2]
    @test_throws ArgumentError MottJainED._global_grid_refine(
        multiwell; mu_min=0.0, mu_max=1.0, mu_count=2,
    )
    @test_throws ArgumentError MottJainED._continuation_grid_refine(
        multiwell; center=0.5, mu_min=0.0, mu_max=1.0,
        local_half_width=0.02, local_count=5, max_expansions=1,
        wide_mode=:unknown,
    )
end

@testset "Configuration" begin
    config = load_config(joinpath(dirname(@__DIR__), "config", "default.toml"))
    @test config["model"]["nm1"] == 5
    @test config["spectrum"]["k"] == 200
    @test config["spectrum"]["l2_values"] == [0, 2, 6]
    @test config["spectrum"]["c2_values"] == [0, 3]
    @test config["spectrum"]["levels_per_block"] == 7
    @test config["spectrum"]["factor"] == 1.0
    @test !haskey(config["spectrum"], "mu_min")
    @test config["density"]["k"] == 3
    @test config["critical"]["k"] == 10
    @test config["critical"]["mu_count"] == 9
    @test !haskey(config["critical"], "coarse_points")
    @test config["fss"]["mu_count"] == 9
    @test config["fss"]["optimize_strategy"] == "size_continuation"
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
    @test fss5_case["fss"]["score_metric"] == "cost"
    @test fss5_case["fss"]["k"] == 15
    @test fss5_case["output"]["run_name"] == "fss5"
    current_config_path = joinpath(dirname(@__DIR__), "config", "my_run.toml")
    current_config = load_config(current_config_path)
    @test current_config["model"]["nm_values"] == [3, 4, 5, 6]
    current_uf0_fss = load_config(
        current_config_path;
        override=joinpath(dirname(@__DIR__), "config", "fss_profiles", "no_w_uf0.toml"),
    )
    current_vf0_fss = load_config(
        current_config_path;
        override=joinpath(dirname(@__DIR__), "config", "fss_profiles", "no_w_vf0.toml"),
    )
    @test current_uf0_fss["fss"]["scan_parameter"] == "Uf0"
    @test current_uf0_fss["fss"]["scan_values"] == [1.2, 1.5, 1.834, 2.2, 2.6]
    @test current_uf0_fss["output"]["run_name"] == "no_w_uf0"
    @test current_vf0_fss["fss"]["scan_parameter"] == "Vf0"
    @test current_vf0_fss["fss"]["scan_values"] == [0.0, 0.2, 0.41, 0.7, 1.0]
    @test current_vf0_fss["output"]["run_name"] == "no_w_vf0"
    @test current_uf0_fss["fss"]["methods"] == ["optimize"]
    @test current_uf0_fss["fss"]["optimize_strategy"] == "size_continuation"
    @test current_uf0_fss["fss"]["optimize_anchor_nm"] == 4
    @test current_uf0_fss["fss"]["optimize_local_count"] == 9
    @test current_uf0_fss["fss"]["mu_min"] == -0.095
    @test current_uf0_fss["fss"]["mu_max"] == 0.305
    mktemp() do _, plan_io
        redirect_stdout(plan_io) do
            MottJainED._plan(fss5_case)
        end
        flush(plan_io)
        seekstart(plan_io)
        plan_text = read(plan_io, String)
        @test occursin("FSS 每 sector 的 k", plan_text)
        @test occursin("= 15", plan_text)
        @test occursin("FSS 外层参数", plan_text)
        @test occursin("FSS 参数点", plan_text)
        @test occursin("FSS optimize 搜索策略", plan_text)
        @test occursin("发现网格 ED", plan_text)
    end
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
