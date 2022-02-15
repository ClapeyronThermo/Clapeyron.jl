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
    ts = T_scales(model)
    pmix = p_scale(model,x)
    if v0 === nothing
        v0 = x0_LLE_pressure(model,T,x)
    end
    
    f! = (F,z) -> Obj_bubble_pressure(model, F, T, exp10(z[1]), exp10(z[2]), x,z[3:end],ts,pmix)
    r  =Solvers.nlsolve(f!,v0,LineSearch(Newton()))
    sol = Solvers.x_sol(r)
    v_l = exp10(sol[1])
    v_ll = exp10(sol[2])
    xx = FractionVector(sol[3:end])
    P_sat = pressure(model,v_l,T,x)
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
    cache = Ref{Tuple{TT,TT,TT,FractionVector{TT,Vector{TT}}}}()
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