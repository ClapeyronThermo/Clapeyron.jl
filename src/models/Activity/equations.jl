#for use in models that have activity coefficient defined.
function recombine_impl!(model::ActivityModel)
    recombine!(model.puremodel)
    return model
end

function excess_gibbs_free_energy(model::ActivityModel,p,T,z)
    γ = activity_coefficient(model,p,T,z)
    return sum(z[i]*R̄*T*log(γ[i]) for i ∈ @comps)
end

function test_excess_gibbs_free_energy(model::ActivityModel,p,T,z)
    γ = activity_coefficient(model,p,T,z)
    return sum(z[i]*R̄*T*log(γ[i]) for i ∈ @comps)
end


#for use in models that have gibbs free energy defined.
function activity_coefficient(model::ActivityModel,p,T,z)
    X = gradient_type(p,T,z)
    return exp.(ForwardDiff.gradient(x->excess_gibbs_free_energy(model,p,T,x),z)/(R̄*T))::X
end

function test_activity_coefficient(model::ActivityModel,p,T,z)
    X = gradient_type(p,T,z)
    return exp.(ForwardDiff.gradient(x->excess_gibbs_free_energy(model,p,T,x),z)/(R̄*T))::X
end

x0_sat_pure(model::ActivityModel,T) = x0_sat_pure(model.puremodel[1],T)

function saturation_pressure(model::ActivityModel,T::Real,method::SaturationMethod)
    return saturation_pressure(model.puremodel[1],T,method)
end

function eos(model::ActivityModel,V,T,z)
    Σz = sum(z)
    lnΣz = log(Σz)
    
    pures = model.puremodel
    p = sum(z[i]*pressure(pures[i],V/Σz,T) for i ∈ @comps)/Σz
    g_E = excess_gibbs_free_energy(model,p,T,z)
    g_ideal = sum(z[i]*R̄*T*(log(z[i])-lnΣz) for i ∈ @comps)
    g_pure = sum(z[i]*VT_gibbs_free_energy(pures[i],V/Σz,T) for i ∈ @comps)
    
    return g_E+g_ideal+g_pure-p*V
end

function eos_res(model::ActivityModel,V,T,z)
    Σz = sum(z)
    pures = model.puremodel
    g_pure_res = sum(z[i]*VT_gibbs_free_energy_res(pures[i],V/Σz,T) for i ∈ @comps)
    p = sum(z[i]*pressure(pures[i],V,T) for i ∈ @comps)/Σz
    g_E = excess_gibbs_free_energy(model,p,T,z)
    p_res = p - Σz*R̄*T/V
    return g_E+g_pure_res-p_res*V
end


function mixing(model::ActivityModel,p,T,z,::typeof(enthalpy))
    f(x) = excess_gibbs_free_energy(model,p,x,z)/x
    df(x) = Solvers.derivative(f,x)
    return -df(T)*T^2
end

function mixing(model::ActivityModel,p,T,z,::typeof(gibbs_free_energy))
    x = z./sum(z)
    return excess_gibbs_free_energy(model,p,T,z)+dot(z,log.(x))*R̄*T
end

function mixing(model::ActivityModel,p,T,z,::typeof(entropy))
    f(x) = excess_gibbs_free_energy(model,p,x,z)/x
    g,dg = Solvers.f∂f(f,T)
    return -dg*T-g
end

function gibbs_solvation(model::ActivityModel,T)
    z = [1.0,1e-30] 
    p,v_l,v_v = saturation_pressure(model.puremodel[1],T)
    p2,v_l2,v_v2 = saturation_pressure(model.puremodel[2],T)
    γ = activity_coefficient(model,p,T,z)
    K = v_v/v_l*γ[2]*p2/p   
    return -R̄*T*log(K)
end

function lb_volume(model::ActivityModel,z = SA[1.0])
    b = sum(lb_volume(model.puremodel[i])*z[i] for i in @comps)
    return b
end

function T_scale(model::ActivityModel,z=SA[1.0])
    prod(T_scale(model.puremodel[i])^1/z[i] for i in @comps)^(sum(z))
end

function p_scale(model::ActivityModel,z=SA[1.0])
    0.33*R̄*T_scale(model,z)/lb_volume(model,z)
end

function x0_volume_liquid(model::ActivityModel,T,z = SA[1.0])
    pures = model.puremodel
    return sum(z[i]*x0_volume_liquid(pures[i],T,SA[1.0]) for i ∈ @comps)
end

include("methods/methods.jl")
