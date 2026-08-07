# CFT score 由“候选关系池 + 每个 workflow 自己选择的 terms”组成。
# critical5/fss7/optimization8 只保留为旧文件的快捷预设，不再限制用户只能三选一。

const CFT_RELATION_SPECS = Dict{Symbol,Any}(
    # upper/lower = (L², C₂, 原始 rank)；lower=nothing 表示减全局基态 E0。
    :ds_s => (
        label="dS-S", upper=(2, 0, 1), lower=(0, 0, 2), target=1.0,
    ),
    :dds_ds => (
        label="ddS-dS", upper=(6, 0, 2), lower=(2, 0, 1), target=1.0,
    ),
    :c2_6 => (
        label="C2=6", upper=(6, 6, 1), lower=(2, 6, 1), target=1.0,
    ),
    :boxs_s => (
        label="boxS-S", upper=(0, 0, 3), lower=(0, 0, 2), target=2.0,
    ),
    # 这一项在旧 optimization.jl 中存在但被注释；现在可以按需单独启用。
    :boxo_o => (
        label="boxO-O", upper=(0, 3, 2), lower=(0, 3, 1), target=2.0,
    ),
    :j => (
        label="J", upper=(2, 3, 1), lower=nothing, target=2.0,
    ),
    :curlj => (
        label="curlJ", upper=(2, 3, 3), lower=nothing, target=3.0,
    ),
    :dj_rank1 => (
        label="dJ(rank1)", upper=(6, 3, 1), lower=nothing, target=3.0,
    ),
    :dj_rank3 => (
        label="dJ(rank3)", upper=(6, 3, 3), lower=nothing, target=3.0,
    ),
    :t_rank1 => (
        label="T(rank1)", upper=(6, 0, 1), lower=nothing, target=3.0,
    ),
    :t_rank2 => (
        label="T(rank2)", upper=(6, 0, 2), lower=nothing, target=3.0,
    ),
)

const CFT_SCORE_PRESETS = Dict(
    :critical5 => [:ds_s, :j, :curlj, :dj_rank1, :t_rank2],
    :fss7 => [:ds_s, :dds_ds, :boxs_s, :j, :curlj, :dj_rank1, :t_rank1],
    :optimization8 => [
        :ds_s, :dds_ds, :c2_6, :boxs_s, :j, :curlj, :dj_rank3, :t_rank1,
    ],
)

function normalize_score_definition(value)
    name = Symbol(lowercase(String(value)))
    aliases = Dict(
        :five => :critical5, :critical => :critical5, :critical5 => :critical5,
        :seven => :fss7, :fss => :fss7, :fss7 => :fss7,
        :eight => :optimization8, :optimization => :optimization8,
        :optimize => :optimization8, :optimization8 => :optimization8,
        :custom => :custom,
    )
    haskey(aliases, name) || throw(ArgumentError(
        "Unknown CFT score '$value'; choose critical5, fss7, optimization8, or custom",
    ))
    return aliases[name]
end

function normalize_score_term(value)
    raw = lowercase(strip(String(value)))
    token = strip(replace(raw, r"[^a-z0-9]+" => "_"), '_')
    aliases = Dict(
        "ds_s" => :ds_s, "ds_minus_s" => :ds_s,
        "dds_ds" => :dds_ds, "dds_minus_ds" => :dds_ds,
        "c2_6" => :c2_6, "c26" => :c2_6,
        "boxs_s" => :boxs_s, "box_s_s" => :boxs_s, "boxs_minus_s" => :boxs_s,
        "boxo_o" => :boxo_o, "box_o_o" => :boxo_o, "boxo_minus_o" => :boxo_o,
        "j" => :j, "curlj" => :curlj, "curl_j" => :curlj,
        "dj_rank1" => :dj_rank1, "dj1" => :dj_rank1,
        "dj_rank3" => :dj_rank3, "dj3" => :dj_rank3,
        "t_rank1" => :t_rank1, "t1" => :t_rank1,
        "t_rank2" => :t_rank2, "t2" => :t_rank2,
    )
    haskey(aliases, token) || throw(ArgumentError(
        "Unknown CFT score term '$value'. Available terms: " *
        join(sort!(String.(collect(keys(CFT_RELATION_SPECS)))), ", "),
    ))
    return aliases[token]
end

function normalize_score_terms(values)
    values === nothing && return nothing
    raw = values isa AbstractVector ? collect(values) : split(String(values), ',')
    isempty(raw) && throw(ArgumentError("score_terms must contain at least one relation"))
    terms = normalize_score_term.(raw)
    length(unique(terms)) == length(terms) || throw(ArgumentError(
        "score_terms must not contain duplicate relations",
    ))
    return terms
end

function resolve_score_terms(definition, terms=nothing)
    definition = normalize_score_definition(definition)
    normalized = normalize_score_terms(terms)
    if normalized !== nothing
        return normalized
    end
    definition == :custom && throw(ArgumentError(
        "score='custom' requires a nonempty score_terms list",
    ))
    return copy(CFT_SCORE_PRESETS[definition])
end

function score_identity(terms)
    selected = normalize_score_terms(terms)
    selected === nothing && throw(ArgumentError("score terms are required"))
    for definition in (:critical5, :fss7, :optimization8)
        preset = CFT_SCORE_PRESETS[definition]
        length(selected) == length(preset) && Set(selected) == Set(preset) &&
            return definition
    end
    return :custom
end

default_score_metric(definition) =
    normalize_score_definition(definition) == :critical5 ? :q : :cost

function normalize_score_metric(value, definition)
    value === nothing && return default_score_metric(definition)
    metric = Symbol(lowercase(String(value)))
    metric in (:q, :cost) || throw(ArgumentError("score_metric must be 'q' or 'cost'"))
    return metric
end

function score_tag(definition, metric, terms=nothing)
    selected = resolve_score_terms(definition, terms)
    identity = score_identity(selected)
    metric = normalize_score_metric(metric, definition)
    if identity == :custom
        digest = bytes2hex(sha1(join(sort(String.(selected)), ",")))[1:10]
        return "custom-$metric-$digest-v2"
    end
    return "$identity-$metric-v2"
end

function _invalid_score(definition, metric, reason; terms=nothing)
    selected = resolve_score_terms(definition, terms)
    return CFTScore(
        definition=score_identity(selected), terms=selected, metric=metric,
        reason=String(reason),
    )
end

function _raw_sector!(cache, states, l2, c2, settings)
    return get!(cache, (Int(l2), Int(c2))) do
        _states_in_sector(states, l2, c2; quantum_tol=settings.quantum_tol)
    end
end

function _relation_energy!(cache, states, reference, ground, settings)
    reference === nothing && return ground, nothing
    l2, c2, rank = reference
    selected = _raw_sector!(cache, states, l2, c2, settings)
    length(selected) >= rank || return NaN, "(L2,C2)=($l2,$c2) rank $rank"
    return selected[rank].energy, nothing
end

"""
计算用户选择的 CFT tower relations。

`definition` 是旧预设快捷名；传入 `terms` 后，以 `terms` 为准，可从候选池任意
选取和组合。每个 term 都保留原始 `(Z,R)` 等能副本和原始 rank。
"""
function cft_score(
    states::Vector{SpectrumState};
    settings::SolverSettings=SolverSettings(),
    definition=:critical5,
    terms=nothing,
    metric=nothing,
)
    definition = normalize_score_definition(definition)
    selected_terms = resolve_score_terms(definition, terms)
    metric = normalize_score_metric(metric, definition)
    identity = score_identity(selected_terms)
    isempty(states) && return _invalid_score(
        definition, metric, "empty spectrum"; terms=selected_terms,
    )

    ground = minimum(state.energy for state in states)
    sector_cache = Dict{Tuple{Int,Int},Vector{SpectrumState}}()
    raw = Float64[]
    targets = Float64[]
    labels = String[]
    missing = String[]
    for term in selected_terms
        spec = CFT_RELATION_SPECS[term]
        upper, upper_missing = _relation_energy!(
            sector_cache, states, spec.upper, ground, settings,
        )
        lower, lower_missing = _relation_energy!(
            sector_cache, states, spec.lower, ground, settings,
        )
        if upper_missing !== nothing || lower_missing !== nothing
            details = filter(value -> value !== nothing, (upper_missing, lower_missing))
            push!(missing, "$(String(term)): $(join(details, " / "))")
            continue
        end
        push!(raw, upper - lower)
        push!(targets, spec.target)
        push!(labels, spec.label)
    end
    if !isempty(missing)
        rejected = count(
            state -> isnothing(_state_quantum_labels(state, settings.quantum_tol)), states,
        )
        detail = rejected == 0 ? "" : "; $rejected states had non-integer quantum numbers"
        return _invalid_score(
            definition, metric, "missing levels: $(join(missing, "; "))$detail";
            terms=selected_terms,
        )
    end

    all(isfinite, raw) || return _invalid_score(
        definition, metric, "non-finite tower gap"; terms=selected_terms,
    )
    raw_norm2 = sum(abs2, raw)
    raw_norm2 > 1.0e-14 || return _invalid_score(
        definition, metric, "zero tower-gap norm"; terms=selected_terms,
    )
    target_norm2 = sum(abs2, targets)
    projection = dot(raw, targets)
    factor = projection / target_norm2
    abs(factor) > eps(Float64) || return _invalid_score(
        definition, metric, "zero fitted factor"; terms=selected_terms,
    )
    residual = raw ./ factor .- targets
    q = sqrt(mean(abs2, residual))
    cost = max(0.0, target_norm2 - projection^2 / raw_norm2)
    objective = metric == :q ? q : cost

    s00 = _raw_sector!(sector_cache, states, 0, 0, settings)
    delta_s = length(s00) >= 2 ? (s00[2].energy - ground) / factor : NaN
    a03 = _raw_sector!(sector_cache, states, 0, 3, settings)
    delta_o = isempty(a03) ? NaN : (a03[1].energy - ground) / factor
    return CFTScore(
        valid=true, definition=identity, terms=selected_terms, metric=metric,
        objective=objective, q=q, cost=cost, factor=factor,
        delta_s=delta_s, delta_o=delta_o, raw_gaps=raw,
        target_gaps=targets, labels=labels,
    )
end

"""只在给定 μ 网格上计算所选 score，然后选 objective 最小点。"""
function scan_mu(
    cache::ModelCache,
    mus;
    definition=:critical5,
    terms=nothing,
    metric=nothing,
)
    definition = normalize_score_definition(definition)
    selected_terms = resolve_score_terms(definition, terms)
    metric = normalize_score_metric(metric, definition)
    grid = Float64.(collect(mus))
    isempty(grid) && throw(ArgumentError("mu grid must not be empty"))
    all(isfinite, grid) || throw(ArgumentError("all mu grid values must be finite"))
    scores = CFTScore[]
    errors = String[]
    total = length(grid)

    for (index, mu) in enumerate(grid)
        @info "CFT mu grid" index total mu definition metric terms=selected_terms
        score = try
            cft_score(
                solve_spectrum(cache, mu); settings=cache.settings,
                definition=definition, terms=selected_terms, metric=metric,
            )
        catch err
            push!(errors, "mu=$mu: $(sprint(showerror, err))")
            _invalid_score(
                definition, metric, sprint(showerror, err); terms=selected_terms,
            )
        end
        push!(scores, score)
        @info "CFT mu result" index total mu valid=score.valid objective=score.objective q=score.q cost=score.cost reason=score.reason
    end

    valid = findall(score -> score.valid && isfinite(score.objective), scores)
    if isempty(valid)
        return (
            mu=NaN,
            score=_invalid_score(
                definition, metric, "no valid score on requested mu grid";
                terms=selected_terms,
            ),
            mus=grid, scores=scores, completed=true, at_boundary=false,
            evaluations=total, invalid=total, errors=errors,
        )
    end
    best_index = valid[argmin([scores[index].objective for index in valid])]
    return (
        mu=grid[best_index], score=scores[best_index], mus=grid, scores=scores,
        completed=true, at_boundary=best_index == 1 || best_index == total,
        evaluations=total, invalid=count(score -> !score.valid, scores), errors=errors,
    )
end

function score_dataframe(score::CFTScore)
    return DataFrame(
        definition=fill(String(score.definition), length(score.labels)),
        metric=fill(String(score.metric), length(score.labels)),
        term=String.(score.terms), label=score.labels,
        raw_gap=score.raw_gaps, target=score.target_gaps,
        scaled_gap=score.raw_gaps ./ score.factor,
        residual=score.raw_gaps ./ score.factor .- score.target_gaps,
    )
end
