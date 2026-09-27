# Included inside SO3lverED. All intermediate states remain in the complete
# projected SO(3) blocks. Only the small candidate-coefficient space is reduced.

function _generator_coordinates(problem, operators; rtol=1e-10)
    !isnothing(problem.generator_coordinates[]) && return problem.generator_coordinates[]
    n = length(problem.candidates[:singlet].names)
    gram = zeros(n, n)
    # Deterministic full-block probes identify operator redundancies, not just
    # directions invisible on the chosen primaries. Two independent probes per
    # domain retain improvements that act only after the first generator.
    for (rep, initial) in sort!(collect(keys(problem.hamiltonians)))
        dim = problem.hamiltonians[(rep, initial)].space.dim
        dim == 0 && continue
        for final in _vector_target_ells(initial)
            haskey(problem.hamiltonians, (rep, final)) || continue
            for probe in 1:2
                v = [sin((j + probe) * sqrt(2.0 + probe)) +
                     cos(j * sqrt(5.0 + probe)) for j in 1:dim]
                v ./= norm(v)
                X = _candidate_action_matrix(v, operators(rep, initial, final))
                gram .+= X' * X
            end
        end
    end
    scales = sqrt.(max.(diag(gram), 0.0))
    active = findall(>(sqrt(eps(Float64)) * maximum(scales)), scales)
    isempty(active) && error("empty generator operator space")
    normalizer = zeros(n, length(active))
    for (j, i) in enumerate(active)
        normalizer[i, j] = inv(scales[i])
    end
    decomposition = eigen(Symmetric(normalizer' * gram * normalizer))
    kept = findall(>(rtol * maximum(decomposition.values)), decomposition.values)
    W = normalizer * decomposition.vectors[:, kept] *
        Diagonal(inv.(sqrt.(decomposition.values[kept])))
    problem.generator_coordinates[] = W
    return W
end

"""Reduced coefficient of [A^(1) × B^(1)]^(k), in actual CG conventions."""
function _coupled_vector_product_coefficient(initial, middle, final, rank)
    abs(initial-rank) <= final <= initial+rank || return 0.0
    for m in -initial:initial, Q in -rank:rank
        abs(m+Q) <= final || continue
        denominator = Float64(clebschgordan(initial,m,rank,Q,final,m+Q))
        abs(denominator) > 1e-10 || continue
        value = 0.0
        for q in -1:1, r in -1:1
            q+r == Q || continue
            value += Float64(clebschgordan(1,q,1,r,rank,Q)) *
                _wigner_eckart_coefficient(initial,m,r,middle) *
                _wigner_eckart_coefficient(middle,m+r,q,final)
        end
        return value / denominator
    end
    return 0.0
end

function _symmetric_pair_columns(matrix, pairs, count)
    hcat((a == b ? view(matrix,:,a+(b-1)*count) :
        view(matrix,:,a+(b-1)*count) + view(matrix,:,b+(a-1)*count)
        for (a,b) in pairs)...)
end

function _compact_commutator_grams(source, rep, ell, energy, ground,
        hamiltonians, operators, W)
    r = size(W, 2)
    pairs = [(a,b) for b in 1:r for a in 1:b]
    features = 3length(pairs)
    gm = zeros(features,features)
    gp = zeros(features,features)
    gk = zeros(features,features)
    cross_d = zeros(features)
    cross_l = zeros(features)
    products = Dict{Tuple{Int,Int},NTuple{4,Matrix{Float64}}}()
    for middle in _vector_target_ells(ell)
        hi = hamiltonians[(rep,middle)]
        X = _candidate_action_matrix(source,operators(rep,ell,middle)) * W
        Y = _apply_hamiltonian_columns(hi,X) - energy*X
        Z = _apply_hamiltonian_columns(hi,Y) - energy*Y
        for final in _vector_target_ells(middle)
            hf = hamiltonians[(rep,final)]
            ops = operators(rep,middle,final)
            A = hcat((_candidate_action_matrix(view(X,:,j),ops)*W for j in 1:r)...)
            B = hcat((_candidate_action_matrix(view(Y,:,j),ops)*W for j in 1:r)...)
            C = hcat((_candidate_action_matrix(view(Z,:,j),ops)*W for j in 1:r)...)
            F10 = _apply_hamiltonian_columns(hf,A) - energy*A - B
            F11 = _apply_hamiltonian_columns(hf,B) - energy*B - C
            products[(middle,final)] = Tuple(
                _symmetric_pair_columns(M,pairs,r) for M in (A,B,F10,F11)
            )
        end
    end
    finals = sort!(unique(last.(collect(keys(products)))))
    for rank in 0:2, final in finals
        abs(ell-rank) <= final <= ell+rank || continue
        dim = hamiltonians[(rep,final)].space.dim
        A,B,C,E = (zeros(dim,length(pairs)) for _ in 1:4)
        for middle in _vector_target_ells(ell)
            haskey(products,(middle,final)) || continue
            alpha = _coupled_vector_product_coefficient(ell,middle,final,rank)
            parts = products[(middle,final)]
            A .+= alpha .* parts[1]; B .+= alpha .* parts[2]
            C .+= alpha .* parts[3]; E .+= alpha .* parts[4]
        end
        zero = zeros(dim,length(pairs))
        # Exchanging the two rank-one tensor indices contributes (-1)^rank.
        mixed = iseven(rank) ? hcat(zero,(B-C)/2,zero) : hcat(A/2,zero,-E/2)
        weight = (2final+1)/(2ell+1)
        gm .+= weight .* (mixed' * mixed)
        if rank == 1
            pp = hcat(A/2,(B+C)/2,E/2)
            kk = hcat(A/2,-(B+C)/2,E/2)
            # The public objective sums i<j, half the complete tensor norm.
            gp .+= (weight/2) .* (pp' * pp)
            gk .+= (weight/2) .* (kk' * kk)
        end
        if final == ell && rank == 0
            cross_d .+= weight .* (mixed' * (-2sqrt(3)*(energy-ground)*source))
        elseif final == ell && rank == 1
            cross_l .+= weight .* (mixed' * (2sqrt(2ell*(ell+1))*source))
        end
    end
    return (mixed=gm, pp=gp, kk=gk, cross_d=cross_d, cross_l=cross_l,
        rhs_d=12(energy-ground)^2, rhs_l=8ell*(ell+1), pairs=pairs)
end

function _compact_first_grams(channels, W, ell)
    r = size(W,2)
    names = (:lambda,:cross,:comm,:second,:lambda_second,
             :leak_lambda,:leak_cross,:leak_comm,
             :null_lambda,:null_cross,:null_comm)
    grams = Dict(name=>zeros(r,r) for name in names)
    for channel in channels
        X = channel.lambda_design * W
        Y = channel.commutator_numerator * W
        Z = channel.second_commutator_numerator * W
        LX = channel.lambda_leakage_design * W
        LY = channel.commutator_leakage_numerator * W
        weight = channel.weight
        for (name, matrix) in ((:lambda,X'*X),(:cross,X'*Y),(:comm,Y'*Y),
                (:second,Z'*Z),(:lambda_second,X'*Z),(:leak_lambda,LX'*LX),
                (:leak_cross,LX'*LY),(:leak_comm,LY'*LY))
            grams[name] .+= weight .* matrix
        end
        if channel.target_ell == ell-1
            grams[:null_lambda] .+= weight .* (X'*X)
            grams[:null_cross] .+= weight .* (X'*Y)
            grams[:null_comm] .+= weight .* (Y'*Y)
        end
    end
    return grams
end

function _compact_source_values(first, second, c, factor)
    q(name) = dot(c,first[name]*c)
    lambda = q(:lambda)
    p = (lambda + 2q(:cross)/factor + q(:comm)/factor^2)/4
    k = (lambda - 2q(:cross)/factor + q(:comm)/factor^2)/4
    d = (q(:second)/factor^4 - 2q(:lambda_second)/factor^2 + lambda)/4
    leakage = (q(:leak_lambda)+2q(:leak_cross)/factor+q(:leak_comm)/factor^2)/4
    shortening = (q(:null_lambda)+2q(:null_cross)/factor+q(:null_comm)/factor^2)/4
    v = [c[a]*c[b] for (a,b) in second.pairs]
    features = vcat(v,v/factor,v/factor^2)
    rhs = second.rhs_d/factor^2 + second.rhs_l
    mixed = dot(features,second.mixed*features) -
        2dot(features,second.cross_d/factor+second.cross_l) + rhs
    if lambda <= 1e-10 || p <= 1e-10 || rhs <= eps(Float64)
        return (lambda_norm2=lambda,primary_k=Inf,dilatation=Inf,
            low_energy_leakage=Inf,shortening=Inf,mixed=Inf,
            p_commutator=Inf,k_commutator=Inf)
    end
    safe = eps(Float64)
    return (lambda_norm2=lambda,
        primary_k=max(0.,k)/max(lambda,safe),
        dilatation=max(0.,d)/max(p,safe),
        low_energy_leakage=max(0.,leakage)/max(p,safe),
        shortening=max(0.,shortening)/max(p,safe),
        mixed=max(0.,mixed)/max(rhs,safe),
        p_commutator=max(0.,dot(features,second.pp*features))/max(rhs,safe),
        k_commutator=max(0.,dot(features,second.kk*features))/max(rhs,safe))
end

function _compact_conformal_score(cache, c, factor, weights, worst_weight)
    values = Float64[]; term_weights = Float64[]
    lambda = 0.0
    for source in cache.sources
        row = _compact_source_values(source.first,source.second,c,factor)
        lambda += row.lambda_norm2
        for term in (:primary_k,:dilatation,:mixed,:p_commutator,:k_commutator,
                     :shortening,:low_energy_leakage)
            term == :primary_k && !source.is_primary && continue
            term == :shortening && !source.is_conserved && continue
            weight = get(weights,term,0.)
            weight > 0 || continue
            push!(values,getproperty(row,term)); push!(term_weights,weight)
        end
    end
    if get(weights,:vacuum,0.) > 0
        push!(values,max(0.,dot(c,cache.vacuum*c))/max(lambda,eps(Float64)))
        push!(term_weights,weights[:vacuum])
    end
    isempty(values) && error("empty full-algebra training objective")
    return dot(values,term_weights)/sum(term_weights)+worst_weight*maximum(values)
end

function _refine_full_conformal_fit(problem, channels, vacuum, states, energies,
        specs, primary_labels, conserved_labels, operators, initial_fit;
        training_labels, term_weights, worst_weight, iterations, starts, gram_rtol)
    iterations > 0 && starts > 0 || throw(ArgumentError("inner fit limits must be positive"))
    W = _generator_coordinates(problem,operators;rtol=gram_rtol)
    selected = filter(spec->spec.label in training_labels,specs)
    length(selected)==length(Set(training_labels)) || error("unknown full-fit training source")
    ground = energies[(:singlet,0)][1]
    sources = map(selected) do spec
        key = (spec.representation,spec.ell)
        source = states[key][:,spec.rank]
        (label=spec.label,is_primary=spec.label in primary_labels,
         is_conserved=spec.label in conserved_labels,
         first=_compact_first_grams(filter(ch->ch.label==spec.label,channels),W,spec.ell),
         second=_compact_commutator_grams(source,spec.representation,spec.ell,
             energies[key][spec.rank],ground,problem.hamiltonians,operators,W))
    end
    V = vacuum*W
    cache = (sources=sources,vacuum=V'*V)
    c0 = W \ initial_fit.coefficients
    lower,upper = initial_fit.factor_bounds
    initial = vcat(c0,log(initial_fit.factor))
    function objective(x)
        log(lower) <= x[end] <= log(upper) || return 1e6 +
            min(abs(x[end]-log(lower)),abs(x[end]-log(upper)))^2
        norm(view(x,1:length(x)-1)) > 1e-10 || return 1e6
        value = _compact_conformal_score(cache,view(x,1:length(x)-1),exp(x[end]),
            term_weights,worst_weight)
        return isfinite(value) ? value : 1e6
    end
    baseline = objective(initial)
    best = copy(initial); best_value = baseline; converged = false
    total_iterations = 0
    for start in 1:starts
        seed = copy(initial)
        start > 1 && (seed[1:end-1] .*= 0.75^(start-1))
        # Same algebra objective, including amplitude and scale. The old
        # generalized fit supplies only deterministic initial directions.
        result = optimize(objective,seed,
            NelderMead(initial_simplex=Optim.AffineSimplexer(a=0.02,b=0.10)),
            Optim.Options(iterations=iterations,g_tol=1e-7,
                f_reltol=1e-8,x_abstol=1e-7,show_trace=false))
        total_iterations += Optim.iterations(result)
        if Optim.minimum(result) <= best_value
            best = Optim.minimizer(result); best_value = Optim.minimum(result)
            converged = Optim.converged(result)
        end
    end
    coefficients = W*best[1:end-1]
    factor = exp(best[end])
    distance = min(log(factor/lower),log(upper/factor))
    return (coefficients=coefficients,factor=factor,factor_at_boundary=distance<1e-5,
        commutator_normalization_valid=all(isfinite,coefficients) && norm(coefficients)>1e-10,
        full_algebra_objective=best_value,full_algebra_initial_objective=baseline,
        full_algebra_converged=converged,full_algebra_iterations=total_iterations,
        full_algebra_coordinate_rank=size(W,2),
        full_algebra_cache=cache,full_algebra_coordinates=best[1:end-1],
        full_algebra_basis=W)
end
