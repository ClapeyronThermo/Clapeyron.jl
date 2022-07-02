## Bubble pressure solver
function x0_bubble_pressure(model::EoSModel,T,x)
    #check each T with T_scale, if treshold is over, replace Pi with inf
    comps = length(model)
    pure = split_model(model)
    crit = crit_pure.(pure)
    T_c = first.(crit)
    V_c = last.(crit)
    _0 = zero(T+first(x))
    nan = _0/_0 
    sat_nan = (nan,nan,nan)
    replaceP = T_c .< T
    sat = fill(sat_nan,comps)
    for i in 1:comps
        if !replaceP[i]
        sat[i] = saturation_pressure(pure[i],T,ChemPotVSaturation(crit = crit[i]))
        end
    end
    P_sat = [tup[1] for tup in sat]
    V_l_sat = [tup[2] for tup in sat]
    V_v_sat = [tup[3] for tup in sat]

    P = zero(T)
    V0_l = zero(T)
    V0_v = zero(T)
    Pi   = zero(x)
    for i in 1:length(x)
        if !replaceP[i]
            Pi[i] = P_sat[i]
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
    model_r,idx_r = index_reduction(model,x)
    if length(model_r)==1
        (P_sat,v_l,v_v) = saturation_pressure(model_r,T)
        return (P_sat,v_l,v_v,x)
    end
    x_r = x[idx_r]
    ts = T_scales(model_r)
    pmix = p_scale(model_r,x_r)
    if v0 === nothing
        v0 = x0_bubble_pressure(model_r,T,x_r)
    end
    len = length(v0[1:end-1])
    #xcache = zeros(eltype(x0),len)
    Fcache = zeros(eltype(v0[1:end-1]),len)
    f!(F,z) = Obj_bubble_pressure(model_r, F, T, exp10(z[1]),exp10(z[2]),x_r,z[3:end],ts,pmix)
    r  =Solvers.nlsolve(f!,v0[1:end-1],LineSearch(Newton()))
    sol = Solvers.x_sol(r)
    v_l = exp10(sol[1])
    v_v = exp10(sol[2])
    y_r = FractionVector(sol[3:end])
    P_sat = pressure(model_r,v_l,T,x_r)
    y = zeros(length(model))
    y[idx_r] = y_r
    return (P_sat, v_l, v_v, y)
end

function Obj_bubble_pressure(model::EoSModel, F, T, v_l, v_v, x, y,ts,ps)
    return μp_equality(model::EoSModel, F, T, v_l, v_v, x, FractionVector(y),ts,ps)
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
function bubble_temperature(model::EoSModel,p,x;v0=nothing)
    model_r,idx_r = index_reduction(model,x)
    if length(model_r)==1
        (T_sat,v_l,v_v) = saturation_temperature(model_r,p,v0)
        return (T_sat,v_l,v_v,x)
    end
    x_r = x[idx_r]
    ts = T_scales(model_r)
    pmix = p_scale(model_r,x_r)
    if v0 === nothing
        v0 = x0_bubble_temperature(model_r,p,x_r)
    end
    
    len = length(v0[1:end-1])
    Fcache = zeros(eltype(v0[1:end-1]),len)
    f!(F,z) = Obj_bubble_temperature(model_r, F, p, z[1], exp10(z[2]), exp10(z[3]), x_r, z[4:end],ts,pmix)
    r  =Solvers.nlsolve(f!,v0[1:end-1],LineSearch(Newton()))
    sol = Solvers.x_sol(r)
    T   = sol[1]
    v_l = exp10(sol[2])
    v_v = exp10(sol[3])
    y_r = FractionVector(sol[4:end])
    y = zeros(length(model))
    y[idx_r] = y_r
    return T, v_l, v_v, y
end

function Obj_bubble_temperature(model::EoSModel, F, p, T, v_l, v_v, x, y,ts,ps)
    F = μp_equality(model::EoSModel, F, T, v_l, v_v, x, FractionVector(y),ts,ps)
    F[end] = (pressure(model,v_l,T,x) - p)/ps
    return F
end

function x0_bubble_temperature(model::EoSModel,p,x)
    comps = length(model)   
    pure = split_model(model)
    crit = crit_pure.(pure)
    
    p_c = [tup[2] for tup in crit]
    V_c = [tup[3] for tup in crit]
    _0 = zero(p+first(x))
    replaceP = p_c .< p
    T_sat = fill(_0,comps)
    V_l_sat = fill(_0,comps)
    V_v_sat = fill(_0,comps)
    for i in 1:comps
        crit_i = crit[i]
        Tci,Pci,Vci = crit_i
        if !replaceP[i]
            Ti,Vli,Vvi = saturation_temperature(pure[i],p,AntoineSaturation(crit = crit_i))
        else
            
            Ti,Vli,Vvi = Tci,Vci,1.2*Vci  
        end
        T_sat[i] = Ti
        V_l_sat[i] = Vli
        V_v_sat[i] = Vvi
    end    
    Tb = extrema(T_sat).*(0.9,1.1)

    V0_l = zero(p)
    V0_v = zero(p)
    f(T) = antoine_bubble(pure,T,x,crit)[1]-p
    fT = Roots.ZeroProblem(f,Tb)

    T0 = Roots.solve(fT,Roots.Order0())
    p,y = antoine_bubble(pure,T0,x,crit)
    for i in 1:length(x)
        if !replaceP[i]
            V0_v += y[i]*V_v_sat[i]
            V0_l += x[i]*V_l_sat[i]
        else 
            V0_v += y[i]*V_c[i]*1.2
            V0_l += x[i]*V_c[i]
        end
    end
    prepend!(y,log10.([V0_l,V0_v]))
    prepend!(y,T0)

    return y
end

function antoine_bubble(pure,T,x,crit)
    pᵢ = aprox_psat.(pure,T,crit)
    p = sum(x.*pᵢ)
    y = x.*pᵢ./p
    ysum = 1/∑(y)
    y    = y.*ysum
    return p,y
end

