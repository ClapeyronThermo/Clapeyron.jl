
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

function initial_points_bd_T(pure,T,_crit = nothing,volatile = true, bubble = true)
    #try without critical point information
    if !volatile
        #the component is not volatile/condensable. set volumes to 0.
        #one does not matter, the other one is recalculated later
        vv = zero(one(eltype(pure))+T)
        vl = zero(vv)
        if bubble
            p = zero(vv)
        else
            p = one(vv)/zero(vv) #inf
        end
        return p,vl,vv
    end
    if _crit === nothing #no crit point available, try calculating sat_p without it
        sat = saturation_pressure(pure,T,crit_retry = false)
        !isnan(first(sat)) && return sat
    end
    #calculate critical point, try again
    if _crit !== nothing
        crit = _crit
    else
        crit = crit_pure(pure)
    end
    Tc,Pc,Vc = crit
    if T < Tc
        sat = saturation_pressure(pure,T,crit = crit) #calculate sat_p with crit info
        !isnan(first(sat)) && return sat
    end
    #create initial point from critical values
    #we use a pseudo-saturation pressure extension,based on the slope at the critical point.
    dlnpdTinv,logp0,Tcinv = __dlnPdTinvsat(pure,sat,crit,T,volatile,false)
    lnp = logp0 + dlnpdTinv*(1/T - Tcinv)
    p0 = exp(lnp)
    vl0 = x0_volume(pure,p0,T,phase = :l)
    vv0 = max(1.2*Vc,3*Rgas(pure)*T/Pc)
    return p0,vl0,vv0
end

#this function does not do the crit calculation.
function initial_points_bd_p(pure,p,volatile = true)
    if !volatile
        #the component is not volatile/condensable. set volumes to 0.
        #one does not matter, the other one is recalculated later
        vv = zero(one(eltype(pure))+p)
        vl = zero(vv)
        T = zero(vv) #TODO
        return T,vl,vv
    end
    return saturation_temperature(pure,p,crit_retry = false)
end


function fix_vi!(pure,vli,p,T,in_media,phase)
    for i in eachindex(vli)
        if !in_media[i]
            # vli[i] = volume(pure[i],p,T,phase = phase)
            # if isnan(vli[i])
            vli[i] = x0_volume_liquid(pure[i],T,[1.])
            # end
        end
    end
end

function __crit_pure(sat0,pure,in_media)
    if isnan(first(sat0)) && in_media
        return crit_pure(pure)
    else
        p,vl,vv = sat0
        nan = zero(p)/zero(p)
        return nan,nan,nan
    end
end

function fix_sat_ti!(sat,pure,crit,p,in_media)
    for i in eachindex(pure)
        crit_i = crit[i]
        pc_i = crit_i[2]
        if in_media[i] && p <= pc_i
            pure_i = pure[i]
            sat_i = saturation_temperature(pure_i,p,crit = crit_i)
            if isnan(sat_i[1])
                throw(error("saturation temperataure for $pure_i not found at p = $p"))
            end
            sat[i] = sat_i
        end  
    end
end

function __dlnPdTinvsat(pure,sat,crit,xx,in_media = true,is_sat_temperature = true)
    if is_sat_temperature
        T,vl,vv = sat
        p = xx
        nan_check = isnan(T)
    else
        p,vl,vv = sat
        T = xx
        nan_check = isnan(p)
    end
    if in_media && !nan_check
        Sl = VT_entropy(pure,vl,T)
        Sv = VT_entropy(pure,vv,T)
        dpdT = (Sv - Sl)/(vv - vl)
        return -dpdT*T*T/p,log(p),1/T
    elseif !in_media
        return zero(vl)
    elseif !isnan(crit[1]) && (crit[1] <= T || crit[2] < p)
        Tc,Pc,Vc = crit
        _p(_T) = pressure(pure,Vc,_T)
        dpdT = Solvers.derivative(_p,Tc)
        return -dpdT*Tc*Tc/Pc,log(Pc),1/Tc
    else
        throw(error("dPdTsat: unreachable state with $pure"))
    end
end

function __x0_bubble_pressure(model::EoSModel,T,x,y0 = nothing,volatiles = FillArrays.Fill(true,length(model)),pure = split_model(model),crit = nothing)
    #check each T with T_scale, if treshold is over, replace Pi with inf
    pure_vals = initial_points_bd_T.(pure,T,crit,volatiles,true) #saturation, or aproximation via critical point.
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
    fix_vi!(pure,vli,p,T,volatiles,:l) #calculate volumes if not-volatiles present
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
                vv = volume(model,p0,T,y0,phase = :v)
            end
        else
            if !isnothing(vol0)
                vl,vv = vol0
                p0 = pressure(model,vv,T,y0)
            else
                p0,_,_,_ = __x0_bubble_pressure(model,T,x,y0,volatiles)
                vl = volume(model,p0,T,x,phase = :l)
                vv = volume(model,p0,T,y0,phase = :v)
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

function __x0_bubble_temperature(model::EoSModel,p,x,Tx0 = nothing,volatiles = FillArrays.Fill(true,length(model)),pure = split_model(model),crit = nothing)
    sat = initial_points_bd_p.(pure,p,volatiles)
    if crit === nothing
        _crit = __crit_pure.(sat,pure,volatiles)
    else
        _crit = crit
    end
    fix_sat_ti!(sat,pure,_crit,p,volatiles)
    if Tx0 !== nothing
        T0 = Tx0
    else
        dPdTsat = __dlnPdTinvsat.(pure,sat,_crit,p,volatiles)
        prob = antoine_bubble_problem(dPdTsat,p,x,volatiles)
        T0 = Roots.solve(prob)
    end
    #this is exactly like __x0_bubble_pressure, but we use T0, instead of an input T
    _,vl0,vv0,y = __x0_bubble_pressure(model,T0,x,nothing,volatiles,pure,crit)
    return T0,vl0,vv0,y
end

function antoine_bubble_problem(dpdt,p_bubble,x,volatiles)  
    function antoine_f0(T)
        p = zero(T+first(x)+first(dpdt)[1])
        for i in 1:length(dpdt)
            dlnpdTinv,logp0,T0inv = dpdt[i]
            if volatiles[i]
                pᵢ = exp(logp0 + dlnpdTinv*(1/T - T0inv))
                pᵢxᵢ = x[i]*pᵢ
                p += pᵢxᵢ
            end
        end
        return p - p_bubble
    end
    Tmin,Tmax = extrema(x -> 1/last(x),dpdt)
    return Roots.ZeroProblem(antoine_f0,(Tmin,Tmax))
end

function x0_bubble_temperature(model::EoSModel,p,x)
    T0,V0_l,V0_v,y = __x0_bubble_temperature(model,p,x)
    prepend!(y,log10.([V0_l,V0_v]))
    prepend!(y,T0)
    return y
end

function bubble_temperature_init(model,p,x,vol0,T0,y0,volatiles)
    if !isnothing(y0)
        if !isnothing(T0)
            if !isnothing(vol0)
                vl,vv = vol0
            else
                vl = volume(model,p,T0,x,phase = :l)
                vv = volume(model,p,T0,y0,phase =:v)
            end
        else
            T0,vl0,vv0,_ = __x0_bubble_temperature(model,p,x,T0,volatiles)
            if !isnothing(vol0)
                vl,vv = vol0
            else
                vl = min(volume(model,p,T0,x,phase = :l),vl0)
                vv = max(volume(model,p,T0,y0,phase =:v),vv0)
            end
        end
    else
        T00,vl0,vv0,y0 = __x0_bubble_temperature(model,p,x,T0,volatiles)
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
