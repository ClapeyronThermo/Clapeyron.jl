abstract type BubblePointMethod <: ThermodynamicMethod end
abstract type BubblePressureMethod <: BubblePointMethod end

function __x0_bubble_pressure(model::EoSModel,T,x)
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
    #prepend!(y,log10.([V0_l,V0_v]))
    return P,V0_l,V0_v,y
end

function x0_bubble_pressure(model,T,x)
    p,V0_l,V0_v,y = __x0_bubble_pressure(model,T,x)
    prepend!(y,log10.([V0_l,V0_v]))
    return y
end

function bubble_pressure_init(model,T,x,vol0,p0,y0)
    if !isnothing(y0)
        if !isnothing(p0)
            if !isnothing(vol0)
                vl,vv = vol0
            else
                vl = volume(model,p0,T,x,phase = :l)
                vv = volume(model,p0,T,y0,phase =:v)
            end
        else
            if !isnothing(vol0)
                vl,vv = vol0
                p0 = pressure(model,vv,T,y)
            else
                p0,_,_,_ = __x0_bubble_pressure(model,T,x)
                vl = volume(model,p0,T,x,phase = :l)
                vv = volume(model,p0,T,y0,phase =:v)
            end
        end
    else
        p00,vl0,vv0,y0 = __x0_bubble_pressure(model,T,x)
        if !isnothing(p0)
            vl = volume(model,p0,T,x,phase = :l)
            vv = volume(model,p0,T,y0,phase = :v)
        else
            vl = vl0
            vv = vv0
            p0 = p00
        end
    end
    return p0,vl,vv,y0
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
function bubble_pressure(model::EoSModel, T, x, method::ThermodynamicMethod)
    model_r,idx_r = index_reduction(model,x)
    if length(model_r)==1
        (P_sat,v_l,v_v) = saturation_pressure(model_r,T)
        n = sum(x)
        return (P_sat,n*v_l,n*v_v,x)
    end
    x_r = x[idx_r]
    (P_sat, v_l, v_v, y_r) = bubble_pressure_impl(model,T,x_r,method)
    y = index_expansion(y_r,idx_r)
    return (P_sat, v_l, v_v, y)
end

#=
"""
    bubble_temperature(model::EoSModel, p, x)

calculates the bubble temperature and properties at a given pressure.
Returns a tuple, containing:
- Bubble Temperature `[K]`
- liquid volume at Bubble Point [`m³`]
- vapour volume at Bubble Point [`m³`]
- Gas composition at Bubble Point
"""
function bubble_temperature(model::EoSModel,p,x,method::ThermodynamicMethod)
    model_r,idx_r = index_reduction(model,x)
    if length(model_r)==1
        (T_sat,v_l,v_v) = saturation_temperature(model_r,p)
        return (T_sat,v_l,v_v,x)
    end
    x_r = x[idx_r]
    (T, v_l, v_v, y_r) = bubble_temperature_impl(model,x_r,T,method)
    y = index_expansion(y_r,idx_r)
    return T, v_l, v_v, y
end
=#
include("bubble_point/bubble_chempot.jl")
include("bubble_point/bubble_fugacity.jl")
include("bubble_point/bubble_fugacity_non_volatile.jl")

#legacy
function bubble_pressure(model::EoSModel,T,x;v0 = nothing)
    if isnothing(v0)
        return bubble_pressure(model,T,x,ChemPotBubblePressure())
    else
        vl = exp10(v0[1])
        vv = exp10(v0[2])
        vol0 = (vl,vv)
        y = v0[3:end]
        bubble_pressure(model,T,x,ChemPotBubblePressure(;vol0,y))
    end
end