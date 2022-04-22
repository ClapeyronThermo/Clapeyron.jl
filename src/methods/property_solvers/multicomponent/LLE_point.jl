## LLE pressure solver
function x0_LLE_pressure(model::EoSModel,T,x)
    xx = Fractions.neg(x)
    pure = split_model(model)
    sat = saturation_pressure.(pure,T)
    V_l_sat = getindex.(sat,2)
    V0_l = dot(x,V_l_sat)/sum(x)
    V0_ll = dot(xx,V_l_sat)
    prepend!(xx,log10.((V0_l,V0_ll)))
    return xx[1:end-1]
end

"""
    LLE_pressure(model::EoSModel, T, x; v0 = x0_LLE_pressure(model,T,x))

calculates the Liquid-Liquid equilibrium pressure and properties at a given temperature.

Returns a tuple, containing:
- LLE Pressure `[Pa]`
- liquid volume of composition `x₁ = x` at LLE Point [`m³`]
- liquid volume of composition `x₂` at LLE Point  [`m³`]
- Liquid composition `x₂`
"""
function LLE_pressure(model::EoSModel, T, x; v0 =nothing)
    model_r,idx_r = index_reduction(model,x)
    if length(model_r)==1
        error("There is no LLE for a pure component")
    end
    x_r = x[idx_r]
    ts = T_scales(model_r)
    pmix = p_scale(model_r,x_r)
    if v0 === nothing
        v0 = x0_LLE_pressure(model_r,T,x_r)
    end
    
    f! = (F,z) -> Obj_bubble_pressure(model_r, F, T, exp10(z[1]), exp10(z[2]), x_r,z[3:end],ts,pmix)
    r  =Solvers.nlsolve(f!,v0,LineSearch(Newton()))
    sol = Solvers.x_sol(r)
    v_l = exp10(sol[1])
    v_ll = exp10(sol[2])
    xx_r = FractionVector(sol[3:end])
    P_sat = pressure(model_r,v_l,T,x_r)
    xx = zeros(length(model))
    xx[idx_r] = xx_r
    return (P_sat, v_l, v_ll, xx)
end

"""
    LLE_temperature(model::EoSModel, p, x; T0 = x0_LLE_temperature(model,p,x))

calculates the Liquid-Liquid equilibrium temperature and properties at a given pressure.

Returns a tuple, containing:
- LLE Pressure `[Pa]`
- liquid volume of composition `x₁ = x` at LLE Point [`m³`]
- liquid volume of composition `x₂` at LLE Point  [`m³`]
- Liquid composition `x₂`
"""
function LLE_temperature(model,p,x;T0=nothing)
    if T0===nothing
        T0 = x0_LLE_temperature(model,p,x)
    end
    TT = promote_type(typeof(p),eltype(x))
    cache = Base.RefValue{Tuple{TT,TT,TT,Vector{TT}}}()
    f(z) = Obj_LLE_temperature(model,z,p,x,cache)
    fT = Roots.ZeroProblem(f,T0)
    Roots.solve(fT,Roots.Order0())
    return cache[]
    #p,v_l,v_ll,xx = LLE_pressure(model,T,x)
    #return T,v_l,v_ll,xx
end

function Obj_LLE_temperature(model,T,p,x,cache)
    p̃,v_l,v_ll,xx = LLE_pressure(model,T,x)
    cache[] = (T,v_l,v_ll,xx)
    return p̃-p
end

function x0_LLE_temperature(model,p,x)
    return  1.5*sum(T_scales(model))/length(x)
end

function presents_LLE(model,p,T)
    pure = split_model(model)
    g_pure = gibbs_free_energy.(pure,p,T)
    
    function mixing_gibbs(x1)
        z = FractionVector(x1)
        log∑z = log(sum(z))
        g_mix = gibbs_free_energy(model,p,T,z)
        for i in 1:length(z)
            g_mix -= z[i]*(g_pure[i] + R̄*T*(log(z[i]) - log∑z))
        end
        return g_mix
    end

    f0(x1) = ForwardDiff.derivative(mixing_gibbs,x1)
    prob = Roots.ZeroProblem(f0,(0.9))
    res = Roots.solve(prob)
    return (res,mixing_gibbs(res))

end
