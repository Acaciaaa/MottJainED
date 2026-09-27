using Test, LinearAlgebra, Random, FuzzifiED, MottJainED
using FuzzifiED.SO3lver
FuzzifiED.NumThreads=1
BLAS.set_num_threads(1)
include(joinpath(@__DIR__,"..","experimental","SO3lverED.jl"))
using .SO3lverED

struct ExactConformalBlock
    ell::Int
end
struct ExactConformalTensor
    initial::Int
    final::Int
end
function SO3lverED.apply_so3_conformal_generator(v::AbstractVector,
        op::ExactConformalTensor,c::AbstractVector,
        hi::ExactConformalBlock,hf::ExactConformalBlock;factor,generator)
    l,lp=hi.ell,hf.ell
    p=lp==l+1 ? sqrt((l+1)*(2l+1)) : 0.
    k=lp==l-1 ? -sqrt(l*(2l+1)) : 0.
    (generator==:p ? p : generator==:k ? k : p+k).*v
end
@testset "Exact free scalar conformal module, including spin" begin
    blocks=Dict((:singlet,l)=>ExactConformalBlock(l) for l in 0:4)
    for l in 0:2
        audit=SO3lverED._mixed_commutator_audit([1.],:singlet,l,l+.5,0.,blocks,
            (rep,i,f)->ExactConformalTensor(i,f),[1.],1.)
        @test audit.fractional_residual<1e-27
        @test audit.p_commutator_fraction<1e-27
        @test audit.k_commutator_fraction<1e-27
    end
end

@testset "Fixed rank continuation distinguishes deformation and crossings" begin
    specs=((label=:S,representation=:singlet,ell=0,rank=2),)
    function snapshot(x)
        angle=1.2x[1]
        st=[1. 0. 0.;0. cos(angle) -sin(angle);0. sin(angle) cos(angle)]
        (states=Dict((:singlet,0)=>st),energies=Dict((:singlet,0)=>[0.,1.,2.]))
    end
    track=SO3lverED._continue_conformal_identity([0.],[1.],snapshot([0.]),
        snapshot([1.]),snapshot,specs;scales=[1.])
    @test track.passed
    @test track.minimum_anchor<0.70
    @test track.minimum>0.90
    function crossing(x)
        swapped=x[1]>.5
        st=swapped ? [1. 0. 0.;0. 0. 1.;0. 1. 0.] : Matrix{Float64}(I,3,3)
        (states=Dict((:singlet,0)=>st),
         energies=Dict((:singlet,0)=>[0.0,1.0 - abs(x[1]-.5),1.0 + abs(x[1]-.5)]))
    end
    rejected=SO3lverED._continue_conformal_identity([0.],[1.],crossing([0.]),
        crossing([1.]),crossing,specs;scales=[1.],maximum_solves=32)
    @test !rejected.passed
    @test rejected.solves<=32
    bad=snapshot([0.]);bad.energies[(:adjoint,0)]=[-.1]
    @test !SO3lverED._conformal_ground_ok(bad)
end

c=Couplings(Uf=.48,U0=4.32,Uf0=1.68,Vf=0.,Vf0=.3549772487682664,V0=0.,t=.5,mu=.22881105175318384)
s=build_workspace(build_so3_model(nm1=4,representation=:singlet);heavy_space_mode=:laughlin13,disp_std=false)
a=build_workspace(build_so3_model(nm1=4,representation=:adjoint);heavy_space_mode=:laughlin13,heavy_space=s.heavy_space,disp_std=false)
descendants=DEFAULT_CONFORMAL_DESCENDANT_SPECS
problem=build_so3_conformal_problem(s,a,c;generator_basis=:local_density,descendant_specs=descendants,disp_std=false)
labels=[:S,:O,:J,:dS,:curlJ]
weights=Dict(term=>1.0 for term in CONFORMAL_OBJECTIVE_TERMS)
weights[:low_energy_leakage]=.25
result=analyze_so3_conformal_algebra(s,a,c;generator_basis=:local_density,
    descendant_specs=descendants,prepared_problem=problem,full_algebra_fit=true,
    fit_primary_labels=[:S,:O,:J],training_labels=labels,term_weights=weights,
    factor_bounds=(.005,.25),inner_iterations=2000,inner_starts=2,disp_std=false)
function operators(rep,i,f)
    get!(problem.operator_cache,(rep,i,f)) do
        build_generator_operators(problem.hamiltonians[(rep,i)].space,
            problem.hamiltonians[(rep,f)].space,problem.candidates[rep];disp_std=false)
    end
end

@testset "Complete density moments in projected SO3 space" begin
    @test length(result.generator_candidate_names)==15
    for rep in (:singlet,:adjoint), (i,f) in ((0,1),(1,1),(1,2),(2,2))
        for (forward,reverse) in zip(operators(rep,i,f).operators,operators(rep,f,i).operators)
            A=Matrix(forward;disp_std=false);B=Matrix(reverse;disp_std=false)
            @test norm(B-(-1.)^(f-i)*sqrt((2*f+1)/(2*i+1))*A')<
                1e-10*max(1.,norm(A),norm(B))
        end
    end
    light=GetDensityObs(4,3;norm_r2=4.)
    heavy=GetDensityObs(10,1;norm_r2=4.)
    cpd=ContactCouple([Laplacian(light;norm_r2=4.),heavy],zeros(Int64,4,2),2)
    cpd=filter(ch->abs(ch.coeff)>1e-13 && all(p->!isempty(GetComponent(
        ch.amd[p],ch.ch[1,p]/2,ch.ch[1,p]/2)),1:2),cpd)
    pieces=Matrix{Float64}[]; targets=Vector{Float64}[]
    for rep in (:singlet,:adjoint), (i,f) in ((0,1),(1,1),(1,2),(2,2))
        hi=problem.hamiltonians[(rep,i)].space;hf=problem.hamiltonians[(rep,f)].space
        sg=BuildSegOperators(hi.sgsp,hf.sgsp,cpd;disp_std=false)
        actual=BuildCompOperator(hi,hf,cpd,sg,2;disp_std=false)
        push!(targets,vec(Matrix(actual;disp_std=false)))
        push!(pieces,hcat((vec(Matrix(o;disp_std=false)) for o in operators(rep,i,f).operators)...))
    end
    B=vcat(pieces...);target=vcat(targets...)
    @test norm(B*(pinv(B;rtol=1e-10)*target)-target)/norm(target)<1e-10
    old=B[:,1:8]
    @test norm(old*(pinv(old;rtol=1e-10)*target)-target)/norm(target)>.70
end

@testset "Identical inner and outer algebra objectives" begin
    fit=result.fit
    direct=score_so3_conformal_algebra(result;labels,term_weights=weights)
    @test direct.objective≈fit.full_algebra_objective atol=1e-8
    @test fit.full_algebra_objective<fit.full_algebra_initial_objective
    @test !(:T in getproperty.(fit.full_algebra_cache.sources,:label))
    @test !any(row->row.label in (:dS,:curlJ) && row.term==:primary_k,direct.rows)
    @test all(row.used_for_fit for row in result.primary_rows if row.label in (:dS,:curlJ))
    @test !only(filter(row->row.label==:T,result.primary_rows)).used_for_fit
    @test fit.scalar_commutator_residuals≈[
        row.kp_commutator_lhs-row.kp_commutator_target for row in result.primary_rows
        if row.ell==0 && row.label in fit.fit_primary_labels]
    @test any(row->row.label==:J && row.term==:shortening,direct.rows)
    holdout=score_so3_conformal_algebra(result;labels=[:T],term_weights=weights)
    @test any(row->row.label==:T && row.term==:shortening,holdout.rows)
    @test all(row->row.label in (:T,:vacuum),holdout.rows)
    rng=MersenneTwister(413)
    for multiplier in (.9,1.1)
        coefficients=fit.full_algebra_coordinates+.02randn(rng,length(fit.full_algebra_coordinates))
        factor=fit.factor*multiplier
        for src in fit.full_algebra_cache.sources
            spec=only(filter(x->x.label==src.label,(DEFAULT_CONFORMAL_PRIMARY_SPECS...,descendants...)))
            key=(spec.representation,spec.ell)
            audit=SO3lverED._mixed_commutator_audit(result.states[key][:,spec.rank],
                spec.representation,spec.ell,result.energies[key][spec.rank],
                result.energies[(:singlet,0)][1],problem.hamiltonians,operators,
                fit.full_algebra_basis*coefficients,factor)
            compact=SO3lverED._compact_source_values(src.first,src.second,coefficients,factor)
            @test compact.mixed≈audit.fractional_residual atol=1e-8 rtol=1e-8
            @test compact.p_commutator≈audit.p_commutator_fraction atol=1e-8 rtol=1e-8
            @test compact.k_commutator≈audit.k_commutator_fraction atol=1e-8 rtol=1e-8
        end
    end
    @test !isfinite(SO3lverED._compact_conformal_score(fit.full_algebra_cache,
        zeros(length(fit.full_algebra_coordinates)),fit.factor,weights,.25))
    @info "full fit" before=fit.full_algebra_initial_objective after=fit.full_algebra_objective converged=fit.full_algebra_converged factor=fit.factor
end
