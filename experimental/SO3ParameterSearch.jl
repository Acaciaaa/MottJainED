module SO3ParameterSearch

using LinearAlgebra
using MottJainED
using ..SO3lverED

const TRACKED_STATE_SPECS = (
    (label=:G, representation=:singlet, ell=0, rank=1),
    (label=:S, representation=:singlet, ell=0, rank=2),
    (label=:boxS, representation=:singlet, ell=0, rank=3),
    (label=:dS, representation=:singlet, ell=1, rank=1),
    (label=:T, representation=:singlet, ell=2, rank=1),
    (label=:ddS, representation=:singlet, ell=2, rank=2),
    (label=:O, representation=:adjoint, ell=0, rank=1),
    (label=:J, representation=:adjoint, ell=1, rank=1),
    (label=:curlJ, representation=:adjoint, ell=1, rank=2),
    (label=:dJ, representation=:adjoint, ell=2, rank=1),
)

function build_cft_problem(nm1::Int, couplings::Couplings; disp_std::Bool=true)
    workspaces = Dict{Symbol,SO3Workspace}()
    hamiltonians = Dict{Tuple{Symbol,Int},SO3Hamiltonian}()
    workspace_seconds = Dict{String,Float64}()
    operator_seconds = Dict{String,Float64}()
    shared_heavy_space = nothing
    for representation in (:singlet, :adjoint)
        started = time()
        model = build_so3_model(; nm1, representation)
        workspace = build_workspace(
            model; heavy_space=shared_heavy_space, disp_std,
        )
        shared_heavy_space = workspace.heavy_space
        workspaces[representation] = workspace
        workspace_seconds[String(representation)] = time() - started
        for ell in 0:2
            started = time()
            hamiltonians[(representation, ell)] = build_hamiltonian(
                workspace, ell, couplings; disp_std,
            )
            operator_seconds["$(representation)_L$(ell)"] = time() - started
        end
    end
    return (
        workspaces=workspaces,
        hamiltonians=hamiltonians,
        workspace_seconds=workspace_seconds,
        operator_seconds=operator_seconds,
    )
end

function solve_cft_blocks!(
    hamiltonians,
    couplings::Couplings;
    k::Int=6,
    tol::Float64=1.0e-8,
    ncv::Int=max(2k, k + 10),
    warm_vectors::AbstractDict=Dict{Tuple{Symbol,Int},Vector{Float64}}(),
    disp_std::Bool=true,
)
    energies = Dict{Tuple{Symbol,Int},Vector{Float64}}()
    vectors = Dict{Tuple{Symbol,Int},Matrix{Float64}}()
    solve_seconds = Dict{String,Float64}()
    for key in CFT_BLOCK_KEYS
        hamiltonian = hamiltonians[key]
        retune!(hamiltonian, couplings)
        started = time()
        block_energies, block_vectors = solve(
            hamiltonian; k, tol, ncv, vectors=true,
            initvec=get(warm_vectors, key, nothing), disp_std,
        )
        solve_seconds["$(key[1])_L$(key[2])"] = time() - started
        energies[key] = block_energies
        vectors[key] = block_vectors
        isempty(block_energies) || (warm_vectors[key] = copy(block_vectors[:, 1]))
    end
    return (energies=energies, vectors=vectors, solve_seconds=solve_seconds)
end

function track_reference_states(
    reference_vectors,
    current_vectors;
    specifications=TRACKED_STATE_SPECS,
    minimum_overlap::Real=0.70,
    require_same_rank::Bool=true,
)
    0 <= minimum_overlap <= 1 || throw(ArgumentError(
        "minimum_overlap must lie in [0,1]",
    ))
    rows = NamedTuple[]
    for specification in specifications
        key = (specification.representation, specification.ell)
        haskey(reference_vectors, key) || throw(ArgumentError(
            "reference vectors are missing block $key",
        ))
        haskey(current_vectors, key) || throw(ArgumentError(
            "current vectors are missing block $key",
        ))
        reference = reference_vectors[key]
        current = current_vectors[key]
        size(reference, 1) == size(current, 1) || throw(DimensionMismatch(
            "block $key changed dimension between reference and current point",
        ))
        specification.rank <= size(reference, 2) || throw(ArgumentError(
            "reference block $key has only $(size(reference, 2)) vectors; " *
            "cannot track $(specification.label) rank $(specification.rank)",
        ))
        specification.rank <= size(current, 2) || throw(ArgumentError(
            "current block $key has only $(size(current, 2)) vectors; " *
            "cannot track $(specification.label) rank $(specification.rank)",
        ))
        overlaps = abs2.(current' * view(reference, :, specification.rank))
        best_rank = argmax(overlaps)
        expected_overlap = overlaps[specification.rank]
        best_overlap = overlaps[best_rank]
        same_rank = best_rank == specification.rank
        passed = expected_overlap >= minimum_overlap &&
                 (!require_same_rank || same_rank)
        push!(rows, (
            label=specification.label,
            representation=specification.representation,
            ell=specification.ell,
            expected_rank=specification.rank,
            best_rank=best_rank,
            expected_overlap=Float64(expected_overlap),
            best_overlap=Float64(best_overlap),
            same_rank=same_rank,
            passed=passed,
        ))
    end
    return (passed=all(row -> row.passed, rows), rows=rows)
end

function _conformal_state(label::Symbol, name::String, level, l2::Int, c2::Int, rank::Int)
    member = first(level.members)
    member.vector === nothing && error("generator audit requires retained eigenvectors")
    return ConformalState(
        label=label, name=name, l2=l2, c2=c2, rank=rank,
        state=member.vector, basis=member.basis, energy=member.energy,
    )
end

function _low_scalar_states(states, l2::Int, count::Int, prefix::String)
    catalog, _ = level_catalog(states)
    levels = get(catalog, (l2, 0), nothing)
    isnothing(levels) && error("generator audit found no (L2,C2)=($l2,0) levels")
    length(levels) >= count || error(
        "generator audit needs $count physical (L2,C2)=($l2,0) levels; " *
        "found $(length(levels))",
    )
    return [
        _conformal_state(Symbol("$(prefix)$(rank)"), "$(prefix) rank $rank",
                         levels[rank], l2, 0, rank)
        for rank in 1:count
    ]
end

"""
Use the conventional Fock-basis generator once at the N=6 anchor to test the
provisional scalar ranks.  The fitted generator uses only S -> dS.  The gate
then asks whether boxS (scalar L=0 rank 3) and ddS (scalar L=2 rank 2) are the
largest non-parent components and whether their expected two-state subspaces
capture the generated vectors.
"""
function audit_scalar_generator(
    nm1::Int,
    couplings::Couplings;
    k::Int=20,
    eig_tol::Float64=1.0e-8,
    low_level_count::Int=6,
    minimum_fit_fidelity::Real=0.95,
    minimum_expected_subspace_overlap::Real=0.75,
    minimum_expected_subspace_fraction::Real=0.80,
)
    low_level_count >= 4 || throw(ArgumentError("low_level_count must be at least 4"))
    model = build_model(; nm1)
    settings = SolverSettings(
        k=k, eig_tol=eig_tol, warm_start=false,
        ncv_extra=max(12, k),
    )
    cache = prepare_spectrum(model, couplings, settings)
    states = solve_spectrum(cache, couplings.mu; keep_vectors=true)
    l0 = _low_scalar_states(states, 0, low_level_count, "scalarL0_")
    l1 = _low_scalar_states(states, 2, 1, "scalarL1_")
    l2 = _low_scalar_states(states, 6, low_level_count, "scalarL2_")

    # Operational labels used by the established tower convention.
    s = l0[2]
    s.label = :S
    ds = l1[1]
    ds.label = :dS
    candidates = generator_candidates(model)
    fit = fit_generator(s, ds, candidates)
    l0_overlap = generator_overlap(
        ds, l0, fit.terms; target_l=0, l2_terms=model.l2,
    )
    l2_overlap = generator_overlap(
        ds, l2, fit.terms; target_l=2, l2_terms=model.l2,
    )

    l0_values = [Float64(l0_overlap.overlaps[state.label]) for state in l0]
    l2_values = [Float64(l2_overlap.overlaps[state.label]) for state in l2]
    l0_expected = l0_values[2] + l0_values[3]
    l2_expected = l2_values[1] + l2_values[2]
    l0_total = sum(l0_values)
    l2_total = sum(l2_values)
    l0_fraction = l0_total > 0 ? l0_expected / l0_total : 0.0
    l2_fraction = l2_total > 0 ? l2_expected / l2_total : 0.0
    boxs_leading = l0_values[3] == maximum(l0_values[[1; 3:low_level_count]])
    dds_leading = l2_values[2] == maximum(l2_values[2:low_level_count])
    passed = fit.fidelity >= minimum_fit_fidelity &&
             l0_expected >= minimum_expected_subspace_overlap &&
             l2_expected >= minimum_expected_subspace_overlap &&
             l0_fraction >= minimum_expected_subspace_fraction &&
             l2_fraction >= minimum_expected_subspace_fraction &&
             boxs_leading && dds_leading
    return (
        passed=passed,
        fit_fidelity=Float64(fit.fidelity),
        numerical_rank=fit.numerical_rank,
        l0_overlaps=l0_values,
        l2_overlaps=l2_values,
        l0_expected_subspace_overlap=l0_expected,
        l2_expected_subspace_overlap=l2_expected,
        l0_expected_subspace_fraction=l0_fraction,
        l2_expected_subspace_fraction=l2_fraction,
        boxs_is_leading_non_s=boxs_leading,
        dds_is_leading_non_t=dds_leading,
        low_level_count=low_level_count,
    )
end

function normalized_residual_jacobian(plus_scores, minus_scores, parameters)
    isempty(parameters) && throw(ArgumentError("parameters must not be empty"))
    first_score = plus_scores[first(parameters)]
    relation_count = length(first_score.residuals)
    jacobian = zeros(Float64, relation_count, length(parameters))
    for (column, parameter) in enumerate(parameters)
        plus = plus_scores[parameter]
        minus = minus_scores[parameter]
        plus.terms == first_score.terms || throw(ArgumentError(
            "plus score for $parameter uses different relations",
        ))
        minus.terms == first_score.terms || throw(ArgumentError(
            "minus score for $parameter uses different relations",
        ))
        jacobian[:, column] .= (plus.residuals .- minus.residuals) ./ 2
    end
    singular_values = svdvals(jacobian)
    cutoff = isempty(singular_values) ? Inf : maximum(singular_values) * 1.0e-7
    numerical_rank = count(>(cutoff), singular_values)
    condition_number = if isempty(singular_values) || last(singular_values) <= cutoff
        Inf
    else
        first(singular_values) / last(singular_values)
    end
    return (
        matrix=jacobian,
        singular_values=singular_values,
        numerical_rank=numerical_rank,
        condition_number=condition_number,
    )
end

export TRACKED_STATE_SPECS, build_cft_problem, solve_cft_blocks!,
       track_reference_states, audit_scalar_generator,
       normalized_residual_jacobian

end
