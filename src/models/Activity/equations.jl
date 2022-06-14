#for use in models that have activity coefficient defined.
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
    g_E = excess_gibbs_free_energy(model,V,T,z)
    pures = model.puremodel
    p = sum(z[i]*pressure(pures[i],V,T) for i ∈ @comps)/Σz
    g_ideal = sum(z[i]*R̄*T*(log(z[i])-lnΣz) for i ∈ @comps)
    g_pure = sum(z[i]*VT_gibbs_free_energy(pures[i],V,T) for i ∈ @comps)
    return g_E+g_ideal+g_pure-p*V
end

function eos_res(model::ActivityModel,V,T,z)
    Σz = sum(z)
    g_E = excess_gibbs_free_energy(model,V,T,z)
    pures = model.puremodel
    g_pure_res = sum(VT_chemical_potential_res(pures[i],V,T)[1]*z[i] for i ∈ @comps)
    p = sum(z[i]*pressure(pures[i],V,T) for i ∈ @comps)/Σz
    p_res = p - Σz*R̄*T/V
    return g_E+g_pure_res-p_res*V
end

function bubble_pressure(model::ActivityModel,T,x)
    sat = saturation_pressure.(model.puremodel,T)
    p_sat = [tup[1] for tup in sat]
    γ     = activity_coefficient(model,1e-4,T,x)
    p     = sum(x.*γ.*p_sat)
    y     = x.*γ.*p_sat ./ p
    return (p,y)
end

function bubble_temperature(model::ActivityModel,p,x;T0=nothing)
    f(z) = Obj_bubble_temperature(model,z,p,x)
    if T0===nothing
        pure = model.puremodel
        sat = saturation_temperature.(pure,p)
        Ti   = zero(x)
        for i ∈ 1:length(x)
            if isnan(sat[i][1])
                Tc,pc,vc = crit_pure(pure[i])
                g(x) = p-pressure(pure[i],vc,x,[1.])
                Ti[i] = Roots.find_zero(g,(Tc))
            else
                Ti[i] = sat[i][1]
            end
        end
        T = Roots.find_zero(f,(minimum(Ti)*0.9,maximum(Ti)*1.1))
    else
        T = Roots.find_zero(f,T0)
    end
    p,y = bubble_pressure(model,T,x)
    return (T,y)
end

function Obj_bubble_temperature(model::ActivityModel,T,p,x)
    sat = saturation_pressure.(model.puremodel,T)
    p_sat = [tup[1] for tup in sat]
    γ     = activity_coefficient(model,1e-4,T,x)
    y     = x.*γ.*p_sat ./ p
    return sum(y)-1
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
    b = maximum([model.puremodel[i].params.b.values[1,1] for i ∈ @comps])
    return b
end

function T_scale(model::ActivityModel,z=SA[1.0])
    n = sum(z)
    invn2 = one(n)/(n*n)
    Ωa,Ωb = ab_consts(model.puremodel[1])
    _a = [model.puremodel[i].params.a.values[1,1] for i in @comps]
    _b = [model.puremodel[i].params.b.values[1,1] for i in @comps]
    a = dot(z, _a)*invn2/Ωa
    b = dot(z, _b)*invn2/Ωb
    return a/b/R̄
end

function p_scale(model::ActivityModel,z=SA[1.0])
    n = sum(z)
    invn2 = (1/n)^2
    Ωa,Ωb = ab_consts(model.puremodel[1])
    _a = [model.puremodel[i].params.a.values[1,1] for i in @comps]
    _b = [model.puremodel[i].params.b.values[1,1] for i in @comps]
    a = dot(z, _a)*invn2/Ωa
    b = dot(z, _b)*invn2/Ωb
    return a/ (b^2) # Pc mean
end

