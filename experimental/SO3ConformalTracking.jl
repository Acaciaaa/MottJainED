# Deterministic fixed-rank continuation along the anchor-to-point line. Global
# anchor overlap is a diagnostic; only local continuity and level isolation gate.
function _conformal_ground_ok(snapshot; tolerance=1e-8)
    ground = snapshot.energies[(:singlet,0)][1]
    all(isempty(es) || es[1] >= ground-tolerance for es in values(snapshot.energies))
end

function _conformal_identity_step(before, after, specs;
        minimum_overlap=0.90, gap_tolerance=1e-8)
    rows = map(specs) do spec
        key = (spec.representation,spec.ell)
        rank = spec.rank
        overlaps = abs2.(after.states[key]' * before.states[key][:,rank])
        expected = overlaps[rank]
        rival = maximum((overlaps[i] for i in eachindex(overlaps) if i!=rank);init=0.)
        energies = after.energies[key]
        gap = minimum((abs(energies[i]-energies[rank]) for i in eachindex(energies)
                       if i!=rank);init=Inf)
        (label=spec.label,overlap=expected,gap=gap,
         passed=expected>=minimum_overlap && expected>rival+1e-10 && gap>gap_tolerance)
    end
    return (passed=all(row.passed for row in rows) && _conformal_ground_ok(after),rows=rows)
end

function _continue_conformal_identity(anchor_values, values, anchor, target,
        snapshot_at, specs; scales, minimum_overlap=0.90, maximum_step=0.10,
        maximum_depth=10, maximum_solves=64, gap_tolerance=1e-8)
    scales = Float64.(scales)
    all(scales .> 0) || error("tracking scales must be positive")
    counts = Ref(0)
    minima = Dict(spec.label=>1.0 for spec in specs)
    reason = Ref("ok")
    function follow(left_values,right_values,left,right,depth)
        _conformal_ground_ok(right) || (reason[]="ground_branch"; return false)
        step = _conformal_identity_step(left,right,specs;minimum_overlap,gap_tolerance)
        distance = maximum(abs.((right_values-left_values)./scales))
        if distance<=maximum_step && step.passed
            for row in step.rows
                minima[row.label] = min(minima[row.label],row.overlap)
            end
            return true
        end
        if depth>=maximum_depth || counts[]>=maximum_solves
            reason[] = "continuation_unresolved"
            return false
        end
        midpoint = (left_values+right_values)/2
        counts[] += 1
        middle = snapshot_at(midpoint)
        follow(left_values,midpoint,left,middle,depth+1) || return false
        return follow(midpoint,right_values,middle,right,depth+1)
    end
    passed = follow(Float64.(anchor_values),Float64.(values),anchor,target,0)
    rows = map(specs) do spec
        key = (spec.representation,spec.ell)
        overlap = abs2(dot(anchor.states[key][:,spec.rank],target.states[key][:,spec.rank]))
        (label=spec.label,representation=spec.representation,ell=spec.ell,rank=spec.rank,
         overlap=overlap,minimum_step_overlap=minima[spec.label],passed=passed)
    end
    return (passed=passed,reason=reason[],solves=counts[],rows=rows,
        minimum=passed ? minimum(Base.values(minima)) : 0.0,
        minimum_anchor=minimum(row.overlap for row in rows))
end
