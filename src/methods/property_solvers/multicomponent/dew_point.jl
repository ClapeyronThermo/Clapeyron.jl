## dew pressure solver
function x0_dew_pressure(model::EoSModel,T,y)
    pure = split_model(model)
    crit = crit_pure.(pure)
    
    T_c = [tup[1] for tup in crit]
    V_c = [tup[3] for tup in crit]
    _0 = zero(T+first(y))
    nan = _0/_0 
    sat_nan = (nan,nan,nan)
    replaceP = ifelse.(T_c .< T,true,false)
    sat = [if !replaceP[i] saturation_pressure(pure[i],T) else sat_nan end for i in 1:length(pure)]
    
    P_sat = [tup[1] for tup in sat]
    V_l_sat = [tup[2] for tup in sat]
    V_v_sat = [tup[3] for tup in sat]
    P⁻¹ = zero(T)
    V0_l = zero(T)
    V0_v = zero(T)
    Pi   = zero(y)
    for i in 1:length(y)
        if !replaceP[i]
            Pi[i] = P_sat[i][1]
            P⁻¹+=y[i]/Pi[i]
            V0_v += y[i]*V_v_sat[i]
        else 
            Pi[i] = pressure(pure[i],V_c[i],T)
            P⁻¹+=y[i]/Pi[i]
            V0_v += y[i]*V_c[i]*1.2
        end
    end
    P = 1/P⁻¹
    x = @. y*P/Pi
    xsum = 1/∑(x)
    x    = x.*xsum
    
    for i in 1:length(y)
        if !replaceP[i]
            V0_l += x[i]*V_l_sat[i]
        else
            V0_l += x[i]*V_c[i]
        end
    end
    prepend!(x,log10.([V0_l,V0_v]))
    return x
end

"""
    dew_pressure(model::EoSModel, T, y; v0 = x0_dew_pressure(model,T,y))

calculates the dew pressure and properties at a given temperature.
Returns a tuple, containing:
- Dew Pressure `[Pa]`
- liquid volume at Dew Point [`m³`]
- vapour volume at Dew Point [`m³`]
- Liquid composition at Dew Point
"""
function dew_pressure(model::EoSModel, T, y; v0 =nothing)
    model_r,idx_r = index_reduction(model,y)
    if length(model_r)==1
        (P_sat,v_l,v_v) = saturation_pressure(model_r,T)
        return (P_sat,v_l,v_v,y)
    end
    y_r = y[idx_r]
    ts = T_scales(model_r)
    pmix = p_scale(model_r,y_r)
    if v0 === nothing
        v0 = x0_dew_pressure(model_r,T,y_r)
    end
    len = length(v0[1:end-1])
    Fcache = zeros(eltype(v0[1:end-1]),len)
    f!(F,z) = Obj_dew_pressure(model_r, F, T, exp10(z[1]), exp10(z[2]), z[3:end],y_r,ts,pmix)
    r  =Solvers.nlsolve(f!,v0[1:end-1],LineSearch(Newton()))
    sol = Solvers.x_sol(r)
    v_l = exp10(sol[1])
    v_v = exp10(sol[2])
    x_r = FractionVector(sol[3:end])
    P_sat = pressure(model_r,v_v,T,y_r)
    x = zeros(length(model))
    x[idx_r] = x_r
    return (P_sat, v_l, v_v, x)
end

function Obj_dew_pressure(model::EoSModel, F, T, v_l, v_v, x, y,ts,ps)
    return μp_equality(model::EoSModel, F, T, v_l, v_v, FractionVector(x), y ,ts,ps)
end
"""
    dew_temperature(model::EoSModel, p, y; T0 = x0_dew_temperature(model,p,y))

calculates the dew temperature and properties at a given pressure.
Returns a tuple, containing:
- Dew Temperature `[K]`
- liquid volume at Dew Point [`m³`]
- vapour volume at Dew Point [`m³`]
- Liquid composition at Dew Point
"""
function dew_temperature(model::EoSModel,p,y;T0=nothing)
    TT = promote_type(typeof(p),eltype(y))
    if T0 === nothing
        T0 = x0_dew_temperature(model,p,y)
    end
    x0 = x0_dew_pressure(model,T0,y)
    x = FractionVector(x0[3:end-1])
    v_l = exp10(x0[1])
    v_v = exp10(x0[2])
    cache = Base.RefValue{Tuple{TT,TT,TT,Vector{TT}}}((T0,v_l,v_v,x))
    f(z) = Obj_dew_temperature(model,z,p,y,cache)
    fT = Roots.ZeroProblem(f,T0)
    T::TT = Roots.solve(fT,Roots.Order0())
    return cache[]
end

function Obj_dew_temperature(model,T,p,y,cache)
    last_result = cache[]
    x0 = collect(last(last_result))
    prepend!(x0,(log10(last_result[2]),log10(last_result[3])))
    p̃,v_l,v_v,x= dew_pressure(model,T,y,v0 = x0)  
    cache[] = (T,v_l,v_v,x)
    return p̃-p
end

x0_dew_temperature(model,p,y) = sat_T_equimix(model,p)