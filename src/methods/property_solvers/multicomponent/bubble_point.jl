
"""
    BubblePointMethod <: ThermodynamicMethod

Abstract type for `bubble_pressure` and `bubble_temperature` routines.

Should at least support passing the `y0` keyword, containing an initial vapour phase, if available.

"""
abstract type BubblePointMethod <: ThermodynamicMethod end


function index_reduction(method::BubblePointMethod,idx_r)
    if hasfield(typeof(method),:y0)
        if !isnothing(method.y0)
            method_r = deepcopy(method)
            y0_new = method.y0[idx_r]
            resize!(method_r.y0,length(y0_new))
            method_r.y0 .= y0_new
            return method_r
        end
    end
    return method
end

function initial_points_bd_T(pure,T,volatile = true)
    #try without critical point information
    if !volatile
        #the component is not volatile. its gas volume is 0
        vv = zero(one(eltype(pure))+T)
        vl = zero(vv)
        p = zero(vv)
        return p,vl,vv
    end
    sat = saturation_pressure(pure,T,crit_retry = false)
    !isnan(first(sat)) && return sat
    
    #calculate critical point, try again
    crit = crit_pure(pure)
    Tc,Pc,Vc = crit
    if T < Tc
        sat = saturation_pressure(pure,T,crit = crit)
        !isnan(first(sat)) && return sat
    end
    #create initial point from critical values (TODO: use pseudocritical_pressure instead?)
    p0 = pressure(pure,Vc,T)
    vl0 = Vc
    vv0 = 1.2Vc
    return p0,vl0,vv0
end

function fix_vli!(pure,vli,p,T,volatiles)
    for i in eachindex(vli)
        if !volatiles[i]
            vlii = volume(pure[i],p,T,phase =:l)
            vli[i] = volume(pure[i],p,T,phase =:l)
        end
    end
end

function __x0_bubble_pressure(model::EoSModel,T,x,y0 = nothing,volatiles = FillArrays.Fill(true,length(model)))
    #check each T with T_scale, if treshold is over, replace Pi with inf
    pure = split_model(model)
    pure_vals = initial_points_bd_T.(pure,T,volatiles) #saturation, or aproximation via critical point.
    p0 = first.(pure_vals)
    vli = getindex.(pure_vals,2)
    vvi = getindex.(pure_vals,3)
    xipi = p0 .* x
    p = sum(xipi)
    if isnothing(y0)
        y = xipi
        y ./= p
    else
        y = y0
    end
    fix_vli!(pure,vli,p,T,volatiles) #calculate volumes if not-volatiles present
    vl0  = dot(vli,x)
    vv0 = dot(vvi,y)
    return p,vl0,vv0,y
end

function x0_bubble_pressure(model,T,x)
    p,V0_l,V0_v,y = __x0_bubble_pressure(model,T,x)
    prepend!(y,log10.([V0_l,V0_v]))
    return y
end

function bubble_pressure_init(model,T,x,vol0,p0,y0,volatiles = FillArrays.Fill(true,length(model)))
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
                p0 = pressure(model,vv,T,y0)
            else
                p0,_,_,_ = __x0_bubble_pressure(model,T,x,y0,volatiles)
                vl = volume(model,p0,T,x,phase = :l)
                vv = volume(model,p0,T,y0,phase =:v)
            end
        end
    else
        p00,vl0,vv0,y0 = __x0_bubble_pressure(model,T,x,nothing,volatiles)
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
    bubble_pressure(model::EoSModel, T, x, method = ChemPotBubblePressure())

Calculates the bubble pressure and properties at a given temperature.
Returns a tuple, containing:
- Bubble Pressure `[Pa]`
- liquid volume at Bubble Point [`m³`]
- vapour volume at Bubble Point [`m³`]
- Gas composition at Bubble Point

By default, uses equality of chemical potentials, via [`ChemPotBubblePressure`](@ref)
"""
function bubble_pressure(model::EoSModel,T,x;kwargs...)
    if keys(kwargs) == (:v0,)
        nt_kwargs = NamedTuple(kwargs)
        v0 = nt_kwargs.v0
        vl = exp10(v0[1])
        vv = exp10(v0[2])
        vol0 = (vl,vv)
        y0 = v0[3:end]
        _kwargs = (;vol0,y0)
        method = init_preferred_method(bubble_pressure,model,_kwargs)
    else
        method = init_preferred_method(bubble_pressure,model,kwargs)
    end
    return bubble_pressure(model, T, x, method)
end

function bubble_pressure(model::EoSModel, T, x, method::BubblePointMethod)
    x = x/sum(x)
    T = float(T)
    model_r,idx_r = index_reduction(model,x)
    if length(model_r)==1
        (P_sat,v_l,v_v) = saturation_pressure(model_r,T)
        return (P_sat,v_l,v_v,x)
    end
    x_r = x[idx_r]
    (P_sat, v_l, v_v, y_r) = bubble_pressure_impl(model_r,T,x_r,index_reduction(method,idx_r))
    y = index_expansion(y_r,idx_r)
    converged = bubbledew_check(v_l,v_v,y,x)
    if converged
        return (P_sat, v_l, v_v, y)
    else
        nan = zero(v_l)/zero(v_l)
        y = y*nan
        return (nan,nan,nan,y)
    end
end

###Bubble Temperature

function __x0_bubble_temperature(model::EoSModel,p,x)
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
    if !any(replaceP) #p < min(pci), proceed with entalphy aproximation:
        dPdTsat = (VT_entropy.(pure,V_v_sat,T_sat) .- VT_entropy.(pure,V_l_sat,T_sat)) ./ (V_v_sat .- V_l_sat)
        y = copy(dPdTsat)
        ##initialization for T
        #= we solve the aproximate problem of finding T such as:
        p = sum(xi*pi(T))
        where pi ≈ p + dpdt(T-T0)
        then: p = sum(xi*pi) = sum(xi*(p + dpdti*(T-Ti)))
        p = sum(xi*p) + sum(xi*dpdti*(T-Ti)) # sum(xi*p) = p
        0 = sum(xi*dpdti*(T-Ti))
        0 = sum(xi*dpdti)*T - sum(xi*dpdti*Ti)
        0 = T* ∑T - ∑Ti
        T = ∑Ti/∑T
        =#
        ∑T = zero(p+first(dPdTsat))
        ∑Ti = zero(∑T)
        for i in eachindex(dPdTsat)
            ∑Ti += x[i]*dPdTsat[i]*(T_sat[i])
            ∑T += x[i]*dPdTsat[i]
        end        
        T0 = ∑Ti/∑T
        sat = saturation_pressure.(pure,T0)
        for i in 1:comps
            V_l_sat[i] = sat[i][2]
            V_v_sat[i] = sat[i][3]
        end
        y .= x .* first.(sat)
        y ./= sum(y) 
    else
        Tb = extrema(T_sat).*(0.9,1.1)
        f(T) = antoine_bubble(pure,T,x,crit)[1]-p
        fT = Roots.ZeroProblem(f,Tb)
        T0 = Roots.solve(fT,Roots.Order0())
        _,y = antoine_bubble(pure,T0,x,crit)
    end
    V0_l = zero(p)
    V0_v = zero(p)

    for i in 1:comps
        if !replaceP[i]
            V0_v += y[i]*V_v_sat[i]
            V0_l += x[i]*V_l_sat[i]
        else
            V0_v += y[i]*V_c[i]*1.2
            V0_l += x[i]*V_c[i]
        end
    end

    return T0,V0_l,V0_v,y
end

function antoine_bubble(pure,T,x,crit)
    pᵢ = aprox_psat.(pure,T,crit)
    p = sum(x.*pᵢ)
    y = x.*pᵢ./p
    ysum = 1/∑(y)
    y    = y.*ysum
    return p,y
end

function x0_bubble_temperature(model::EoSModel,p,x)
    T0,V0_l,V0_v,y = __x0_bubble_temperature(model,p,x)
    prepend!(y,log10.([V0_l,V0_v]))
    prepend!(y,T0)
    return y
end

function bubble_temperature_init(model,p,x,vol0,T0,y0)
    if !isnothing(y0)
        if !isnothing(T0)
            if !isnothing(vol0)
                vl,vv = vol0
            else
                vl = volume(model,p,T0,x,phase = :l)
                vv = volume(model,p,T0,y0,phase =:v)
            end
        else
            T0,vl0,vv0,_ = __x0_bubble_temperature(model,p,x)
            if !isnothing(vol0)
                vl,vv = vol0
            else
                vl = min(volume(model,p,T0,x,phase = :l),vl0)
                vv = max(volume(model,p,T0,y0,phase =:v),vv0)
            end
        end
    else
        T00,vl0,vv0,y0 = __x0_bubble_temperature(model,p,x)
        if !isnothing(T0)
            vl = min(vl0,volume(model,p,T0,x,phase = :l))
            vv = max(vv0,volume(model,p,T0,y0,phase = :v))
        else
            vl = vl0
            vv = vv0
            T0 = T00
        end
    end
    return T0,vl,vv,y0
end

"""
    bubble_temperature(model::EoSModel, p, x,method::BubblePointMethod = ChemPotBubbleTemperature())

Calculates the bubble temperature and properties at a given pressure.
Returns a tuple, containing:
- Bubble Temperature `[K]`
- liquid volume at Bubble Point [`m³`]
- vapour volume at Bubble Point [`m³`]
- Gas composition at Bubble Point

By default, uses equality of chemical potentials, via [`ChemPotBubbleTemperature`](@ref)
"""
function bubble_temperature(model::EoSModel,p,x;kwargs...)
    if keys(kwargs) == (:v0,)
        nt_kwargs = NamedTuple(kwargs)
        v0 = nt_kwargs.v0
        T0 = v0[1]
        vl = exp10(v0[2])
        vv = exp10(v0[3])
        vol0 = (vl,vv)
        y0 = v0[4:end]
        _kwargs = (;T0,vol0,y0)
        method = init_preferred_method(bubble_temperature,model,_kwargs)
    else
        method = init_preferred_method(bubble_temperature,model,kwargs)
    end
    return bubble_temperature(model,p,x,method)
end

function bubble_temperature(model::EoSModel, p , x, T0::Number)
    kwargs = (;T0)
    method = init_preferred_method(bubble_temperature,model,kwargs)
    return bubble_temperature(model,p,x,method)
end

function bubble_temperature(model::EoSModel, p , x, method::BubblePointMethod)
    x = x/sum(x)
    p = float(p)
    model_r,idx_r = index_reduction(model,x)
    if length(model_r)==1
        (T_sat,v_l,v_v) = saturation_temperature(model_r,p)
        return (T_sat,v_l,v_v,x)
    end
    x_r = x[idx_r]
    (T_sat, v_l, v_v, y_r) = bubble_temperature_impl(model_r,p,x_r,index_reduction(method,idx_r))
    y = index_expansion(y_r,idx_r)
    converged = bubbledew_check(v_l,v_v,y,x)
    if converged
        return (T_sat, v_l, v_v, y)
    else
        nan = zero(v_l)/zero(v_l)
        y = y*nan
        return (nan,nan,nan,y)
    end
end

include("bubble_point/bubble_activity.jl")
include("bubble_point/bubble_chempot.jl")
include("bubble_point/bubble_fugacity.jl")


#default initializers

function init_preferred_method(method::typeof(bubble_pressure),model::EoSModel,kwargs)
    return ChemPotBubblePressure(;kwargs...)
end

function init_preferred_method(method::typeof(bubble_temperature),model::EoSModel,kwargs)
    return ChemPotBubbleTemperature(;kwargs...) 
end
