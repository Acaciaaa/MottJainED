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

function parameter_search_to_physical(values, center, scales)
    length(values) == length(center) == length(scales) ||
        throw(DimensionMismatch("values, center, and scales must have equal length"))
    all(>(0), scales) || throw(ArgumentError("all parameter scales must be positive"))
    return Float64.(center) .+ Float64.(values) .* Float64.(scales)
end

function physical_to_parameter_search(values, center, scales)
    length(values) == length(center) == length(scales) ||
        throw(DimensionMismatch("values, center, and scales must have equal length"))
    all(>(0), scales) || throw(ArgumentError("all parameter scales must be positive"))
    return (Float64.(values) .- Float64.(center)) ./ Float64.(scales)
end

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

function _low_scalar_states(
    catalog, l2::Int, maximum_count::Int, minimum_count::Int, prefix::String,
)
    levels = get(catalog, (l2, 0), nothing)
    isnothing(levels) && error("generator audit found no (L2,C2)=($l2,0) levels")
    length(levels) >= minimum_count || error(
        "generator audit needs at least $minimum_count physical " *
        "(L2,C2)=($l2,0) levels; " *
        "found $(length(levels))",
    )
    count = min(maximum_count, length(levels))
    selected = [
        _conformal_state(Symbol("$(prefix)$(rank)"), "$(prefix) rank $rank",
                         levels[rank], l2, 0, rank)
        for rank in 1:count
    ]
    return (states=selected, available_count=length(levels))
end

function assess_scalar_generator_overlaps(
    l0_values::AbstractVector{<:Real},
    l2_values::AbstractVector{<:Real};
    minimum_expected_subspace_overlap::Real=0.75,
    minimum_expected_subspace_fraction::Real=0.80,
)
    length(l0_values) >= 4 || throw(ArgumentError(
        "at least four L=0 scalar overlaps are required",
    ))
    length(l2_values) >= 3 || throw(ArgumentError(
        "at least three L=2 scalar overlaps are required",
    ))
    all(>=(0), l0_values) && all(>=(0), l2_values) || throw(ArgumentError(
        "generator overlaps must be nonnegative",
    ))
    l0 = Float64.(l0_values)
    l2 = Float64.(l2_values)
    l0_total = sum(l0)
    l2_total = sum(l2)
    l0_unresolved = max(0.0, 1.0 - l0_total)
    l2_unresolved = max(0.0, 1.0 - l2_total)
    l0_expected = l0[2] + l0[3]
    l2_expected = l2[1] + l2[2]
    l0_fraction = l0_total > 0 ? l0_expected / l0_total : 0.0
    l2_fraction = l2_total > 0 ? l2_expected / l2_total : 0.0

    # The unresolved total is an upper bound on every omitted state's overlap.
    # Including it as one adversarial competitor prevents a truncated catalog
    # from making boxS or ddS look artificially dominant.
    boxs_competitor = maximum([l0[1]; l0[4:end]; l0_unresolved])
    dds_competitor = maximum([l2[3:end]; l2_unresolved])
    boxs_margin = l0[3] - boxs_competitor
    dds_margin = l2[2] - dds_competitor
    boxs_leading = boxs_margin >= -1.0e-12
    dds_leading = dds_margin >= -1.0e-12
    passed = l0_expected >= minimum_expected_subspace_overlap &&
             l2_expected >= minimum_expected_subspace_overlap &&
             l0_fraction >= minimum_expected_subspace_fraction &&
             l2_fraction >= minimum_expected_subspace_fraction &&
             boxs_leading && dds_leading
    return (
        passed=passed,
        l0_expected_subspace_overlap=l0_expected,
        l2_expected_subspace_overlap=l2_expected,
        l0_expected_subspace_fraction=l0_fraction,
        l2_expected_subspace_fraction=l2_fraction,
        l0_resolved_overlap=l0_total,
        l2_resolved_overlap=l2_total,
        l0_unresolved_overlap_upper_bound=l0_unresolved,
        l2_unresolved_overlap_upper_bound=l2_unresolved,
        boxs_competitor_overlap_upper_bound=boxs_competitor,
        dds_competitor_overlap_upper_bound=dds_competitor,
        boxs_leading_margin=boxs_margin,
        dds_leading_margin=dds_margin,
        boxs_is_leading_non_s=boxs_leading,
        dds_is_leading_non_t=dds_leading,
    )
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
    minimum_l0_level_count::Int=4,
    minimum_l2_level_count::Int=3,
    minimum_fit_fidelity::Real=0.95,
    minimum_expected_subspace_overlap::Real=0.75,
    minimum_expected_subspace_fraction::Real=0.80,
)
    low_level_count >= minimum_l0_level_count >= 4 || throw(ArgumentError(
        "require low_level_count >= minimum_l0_level_count >= 4",
    ))
    low_level_count >= minimum_l2_level_count >= 3 || throw(ArgumentError(
        "require low_level_count >= minimum_l2_level_count >= 3",
    ))
    model = build_model(; nm1)
    settings = SolverSettings(
        k=k, eig_tol=eig_tol, warm_start=false,
        ncv_extra=max(12, k),
    )
    cache = prepare_spectrum(model, couplings, settings)
    states = solve_spectrum(cache, couplings.mu; keep_vectors=true)
    catalog, _ = level_catalog(states)
    l0_selection = _low_scalar_states(
        catalog, 0, low_level_count, minimum_l0_level_count, "scalarL0_",
    )
    l1_selection = _low_scalar_states(catalog, 2, 1, 1, "scalarL1_")
    l2_selection = _low_scalar_states(
        catalog, 6, low_level_count, minimum_l2_level_count, "scalarL2_",
    )
    l0 = l0_selection.states
    l1 = l1_selection.states
    l2 = l2_selection.states

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
    assessment = assess_scalar_generator_overlaps(
        l0_values, l2_values;
        minimum_expected_subspace_overlap,
        minimum_expected_subspace_fraction,
    )
    passed = fit.fidelity >= minimum_fit_fidelity && assessment.passed
    return (
        passed=passed,
        fit_fidelity=Float64(fit.fidelity),
        numerical_rank=fit.numerical_rank,
        l0_overlaps=l0_values,
        l2_overlaps=l2_values,
        l0_expected_subspace_overlap=assessment.l0_expected_subspace_overlap,
        l2_expected_subspace_overlap=assessment.l2_expected_subspace_overlap,
        l0_expected_subspace_fraction=assessment.l0_expected_subspace_fraction,
        l2_expected_subspace_fraction=assessment.l2_expected_subspace_fraction,
        l0_resolved_overlap=assessment.l0_resolved_overlap,
        l2_resolved_overlap=assessment.l2_resolved_overlap,
        l0_unresolved_overlap_upper_bound=assessment.l0_unresolved_overlap_upper_bound,
        l2_unresolved_overlap_upper_bound=assessment.l2_unresolved_overlap_upper_bound,
        boxs_competitor_overlap_upper_bound=assessment.boxs_competitor_overlap_upper_bound,
        dds_competitor_overlap_upper_bound=assessment.dds_competitor_overlap_upper_bound,
        boxs_leading_margin=assessment.boxs_leading_margin,
        dds_leading_margin=assessment.dds_leading_margin,
        boxs_is_leading_non_s=assessment.boxs_is_leading_non_s,
        dds_is_leading_non_t=assessment.dds_is_leading_non_t,
        low_level_count=low_level_count,
        minimum_l0_level_count=minimum_l0_level_count,
        minimum_l2_level_count=minimum_l2_level_count,
        l0_level_count=length(l0),
        l2_level_count=length(l2),
        l0_available_level_count=l0_selection.available_count,
        l2_available_level_count=l2_selection.available_count,
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
       assess_scalar_generator_overlaps, normalized_residual_jacobian,
       parameter_search_to_physical,
       physical_to_parameter_search

end
