## Bubble pressure solver
function x0_bubble_pressure(model::EoSModel,T,x)
    #check each T with T_scale, if treshold is over, replace Pi with inf
    pure = split_model(model)
    crit = crit_pure.(pure)
    T_c = [tup[1] for tup in crit]
    V_c = [tup[3] for tup in crit]
    _0 = zero(T+first(x))
    nan = _0/_0 
    sat_nan = (nan,nan,nan)
    replaceP = ifelse.(T_c .< T,true,false)
    sat = [if !replaceP[i] saturation_pressure(pure[i],T) else sat_nan end for i in 1:length(pure)]
    
    P_sat = [tup[1] for tup in sat]
    V_l_sat = [tup[2] for tup in sat]
    V_v_sat = [tup[3] for tup in sat]

    P = zero(T)
    V0_l = zero(T)
    V0_v = zero(T)
    Pi   = zero(x)
    for i in 1:length(x)
        if !replaceP[i]
            Pi[i] = P_sat[i][1]
            P+=x[i]*Pi[i]
            V0_l += x[i]*V_l_sat[i]
        else 
            Pi[i] = pressure(pure[i],V_c[i],T)
            P+=x[i]*Pi[i]
            V0_l += x[i]*V_c[i]
        end
    end

    y = @. x*Pi/P
    ysum = 1/∑(y)
    y    = y.*ysum
    
    for i in 1:length(x)
        if !replaceP[i]
            V0_v += y[i]*V_v_sat[i]
        else
            V0_v += y[i]*V_c[i]*1.2
        end
    end
    
    prepend!(y,log10.([V0_l,V0_v]))
    return y
end
"""
    bubble_pressure(model::EoSModel, T, x; v0 = x0_bubble_pressure(model,T,x))

calculates the bubble pressure and properties at a given temperature.
Returns a tuple, containing:
- Bubble Pressure `[Pa]`
- liquid volume at Bubble Point [`m³`]
- vapour volume at Bubble Point [`m³`]
- Gas composition at Bubble Point
"""
function bubble_pressure(model::EoSModel, T, x; v0 =nothing)
    ts = T_scales(model)
    pmix = p_scale(model,x)
    if v0 === nothing
        v0 = x0_bubble_pressure(model,T,x)
    end
    len = length(v0[1:end-1])
    #xcache = zeros(eltype(x0),len)
    Fcache = zeros(eltype(v0[1:end-1]),len)
    f!(F,z) = Obj_bubble_pressure(model, F, T, exp10(z[1]),exp10(z[2]),x,z[3:end],ts,pmix)
    r  =Solvers.nlsolve(f!,v0[1:end-1],LineSearch(Newton()))
    sol = Solvers.x_sol(r)
    v_l = exp10(sol[1])
    v_v = exp10(sol[2])
    y = FractionVector(sol[3:end])
    P_sat = pressure(model,v_l,T,x)
    return (P_sat, v_l, v_v, y)
end

function Obj_bubble_pressure(model::EoSModel, F, T, v_l, v_v, x, y,ts,ps)
    return μp_equality(model::EoSModel, F, T, v_l, v_v, x, FractionVector(y),ts,ps)
end

function Obj_bubble_temperature(model,T,p,x)
    p̃,v_l,v_v,y = bubble_pressure(model,T,x)
    return p̃-p
end

"""
    bubble_temperature(model::EoSModel, p, x; T0 = x0_bubble_pressure(model,p,x))

calculates the bubble temperature and properties at a given pressure.
Returns a tuple, containing:
- Bubble Temperature `[K]`
- liquid volume at Bubble Point [`m³`]
- vapour volume at Bubble Point [`m³`]
- Gas composition at Bubble Point
"""
function bubble_temperature(model,p,x;T0=nothing)
    f(z) = Obj_bubble_temperature(model,z,p,x)
    if T0 === nothing
        T0 = x0_bubble_temperature(model,p,x)
    end
    fT = Roots.ZeroProblem(f,T0)
    T = Roots.solve(fT,Roots.Order0())
    p,v_l,v_v,y = bubble_pressure(model,T,x)
    return T,v_l,v_v,y
end

function x0_bubble_temperature(model,p,x)
    Ti = _sat_Ti(model,p)
    Tmin,Tmax = extrema(Ti)
    return (0.9*Tmin,1.1*Tmax)
end
