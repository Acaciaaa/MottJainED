using Test
using TOML
using MottJainED

include(joinpath(@__DIR__, "..", "experimental", "EndpointEntanglement.jl"))
using .EndpointEntanglement

@testset "Targeted endpoint entanglement" begin
    @test fiqh_root_lz2(6, 3) == -27
    @test laughlin_root_lz2(16, 3) == -27

    # N=2 is only a smoke test of the full ED -> RSES/OES -> plot pipeline.
    # The production N=6 calculation is deliberately never launched by tests.
    mktempdir() do directory
        output = joinpath(directory, "output")
        config_path = joinpath(directory, "endpoint_n2.toml")
        config = Dict{String,Any}(
            "model" => Dict("nm1" => 2),
            "hamiltonian" => Dict(
                "Uf" => 0.46, "Uf0" => 1.834, "U0" => 4.14,
                "Vf" => 0.0, "Vf0" => 0.41, "V0" => 0.525,
                "t" => 0.5, "mu" => 0.0,
            ),
            "solver" => Dict(
                "k" => 2, "eig_tol" => 1e-10, "dense_cutoff" => 128,
                "ncv_extra" => 4, "warm_start" => false,
            ),
            "endpoint_entanglement" => Dict(
                "mu_left" => 0.02, "mu_right" => 0.22,
                "qa" => 3, "f3a" => 0, "f8a" => 0,
                "lz2_values" => [-3],
                "left_expected_counts" => [1],
                "right_expected_counts" => [1],
                "cut_x" => 0.5, "nm1_a" => 1, "nm0_a" => 2,
                "lambda_plot_min" => 1e-14, "plot_xi_max" => 20.0,
                "phase_fraction_min" => 0.0, "output" => output,
            ),
        )
        open(config_path, "w") do io
            TOML.print(io, config; sorted=true)
        end

        spec = load_spec(config_path)
        @test occursin("N=2", plan(spec))
        result = run_endpoint_entanglement(spec)
        @test size(result.ground, 1) == 2
        @test size(result.left, 1) > 0
        @test size(result.right, 1) > 0
        @test result.diagnostics.expected_count == [1, 1]
        @test all(result.diagnostics.dim_a .> 0)
        @test all(result.diagnostics.dim_b .> 0)
        @test isfile(joinpath(output, "endpoint_entanglement_comparison.png"))
        @test isfile(joinpath(output, "left_fiqh_rses_spectrum.csv"))
        @test isfile(joinpath(output, "right_laughlin_oes_spectrum.csv"))
        @test isfile(joinpath(output, "ground_states.jld2"))
        @test isfile(joinpath(output, "completed.toml"))

        # A restart must reuse both endpoint vectors and the per-sector SVD files.
        restarted = run_endpoint_entanglement(spec)
        @test restarted.left.lambda ≈ result.left.lambda atol=0 rtol=0
        @test restarted.right.lambda ≈ result.right.lambda atol=0 rtol=0
    end
end
