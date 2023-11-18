
"""
    DewPointMethod <: ThermodynamicMethod

Abstract type for `dew_pressure` and `dew_temperature` routines.

Should at least support passing the `x0` keyword, containing an initial vapour phase, if available.

"""
abstract type DewPointMethod <: ThermodynamicMethod end


function index_reduction(method::DewPointMethod,idx_r)
    if hasfield(typeof(method),:x0)
        if !isnothing(method.x0)
            method_r = deepcopy(method)
            x0_new = method.x0[idx_r]
            resize!(method_r.x0,length(x0_new))
            method_r.x0 .= x0_new
            return method_r
        end
    end
    return method
end

function __x0_dew_pressure(model::EoSModel,T,y,x0=nothing,condensables = FillArrays.Fill(true,length(model)),pure = split_model(model), crit = nothing)
    pure_vals = initial_points_bd_T.(pure,T,crit,condensables,false) #saturation, or aproximation via critical point.
    p0 = first.(pure_vals)
    vli = getindex.(pure_vals,2)
    vvi = getindex.(pure_vals,3)
    yipi = y ./ p0
    p = 1/sum(yipi)
    if isnothing(x0)
        x = yipi
        x .*= p
    else
        x = x0
    end
    fix_vi!(pure,vvi,p,T,condensables,:v) #calculate volumes if not-condensables present
    vl0  = dot(vli,x)
    vv0 = dot(vvi,y)
    return p,vl0,vv0,x
end

function x0_dew_pressure(model::EoSModel,T,y)
    P,V0_l,V0_v,x = __x0_dew_pressure(model,T,y)
    prepend!(x,log10.([V0_l,V0_v]))
    return x
end

function dew_pressure_init(model,T,y,vol0,p0,x0,condensables)
    if !isnothing(x0)
        if !isnothing(p0)
            if !isnothing(vol0)
                vl,vv = vol0
            else
                vl = volume(model,p0,T,x0,phase = :l)
                vv = volume(model,p0,T,y,phase =:v)
            end
        else
            if !isnothing(vol0)
                vl,vv = vol0
                p0 = pressure(model,vv,T,y)
            else
                p0,vl0,vv0,_ = __x0_dew_pressure(model,T,y,x0,condensables)
                vl = min(vl0,volume(model,p0,T,x0,phase = :l))
                vv = max(vv0,volume(model,p0,T,y,phase =:v))
            end
        end
    else
        p00,vl0,vv0,x0 = __x0_dew_pressure(model,T,y,nothing,condensables)
        if !isnothing(p0)
            vl = min(vl0,volume(model,p0,T,x0,phase = :l))
            vv = max(vv0,volume(model,p0,T,y,phase = :v))
        else
            vl = vl0
            vv = vv0
            p0 = p00
        end
    end
    return p0,vl,vv,x0
end

"""
    dew_pressure(model::EoSModel, T, y,method = ChemPotDewPressure())

Calculates the dew pressure and properties at a given temperature.
Returns a tuple, containing:
- Dew Pressure `[Pa]`
- liquid volume at Dew Point [`m³`]
- vapour volume at Dew Point [`m³`]
- Liquid composition at Dew Point

By default, uses equality of chemical potentials, via [`ChemPotDewPressure`](@ref)
"""
function dew_pressure(model::EoSModel,T,x;kwargs...)
    if keys(kwargs) == (:v0,)
        nt_kwargs = NamedTuple(kwargs)
        v0 = nt_kwargs.v0
        vl = exp10(v0[1])
        vv = exp10(v0[2])
        vol0 = (vl,vv)
        x0 = v0[3:end]
        _kwargs = (;vol0,x0)
        method = init_preferred_method(dew_pressure,model,_kwargs)
    else
        method = init_preferred_method(dew_pressure,model,kwargs)
    end
    return dew_pressure(model, T, x, method)
end

function dew_pressure(model::EoSModel, T, y,method::DewPointMethod)
    y = y/sum(y)
    T = float(T)
    model_r,idx_r = index_reduction(model,y)
    if length(model_r)==1
        (P_sat,v_l,v_v) = saturation_pressure(model_r,T)
        return (P_sat,v_l,v_v,y)
    end
    y_r = y[idx_r]
    (P_sat, v_l, v_v, x_r) = dew_pressure_impl(model_r,T,y_r,index_reduction(method,idx_r))
    x = index_expansion(x_r,idx_r)
    converged = bubbledew_check(v_l,v_v,y,x)
    if converged
        return (P_sat, v_l, v_v, x)
    else
        nan = zero(v_l)/zero(v_l)
        x = x*nan
        return (nan,nan,nan,x)
    end
end

function antoine_dew(pure,T,y,crit,condensables)
    pinv = zero(T+first(x)+first(crit)[1])
    for i in 1:length(pure)
        if condensables[i] #if non-volatile, pi -> 0
            pᵢ = aprox_psat(pure,T,crit)
            pᵢyᵢ = y[i]/pᵢ
            pinv += pᵢyᵢ
        end
    end
    p = 1/pinv
    return p
end

function __x0_dew_temperature(model::EoSModel,p,y,Tx0 = nothing,condensables = FillArrays.Fill(true,length(model)))
    pure = split_model(model)
    crit = crit_pure.(pure)
    sat = initial_points_bd_p.(pure,crit,p,condensables,false)
    
    if Tx0 !== nothing
        T0 = Tx0
    elseif !any(crit_i -> crit_i[2] < p,crit) #p < min(pci), proceed with entalphy aproximation:
        dPdTsat = __dPdTsat.(pure,sat,condensables,false)
        ##initialization for T, dew form
        #= we solve the aproximate problem of finding T such as:
        p = sum(yi*pi(T))
        where pi ≈ p + dpdt(T-T0)
        for a dew specification:
        sum(xi*pi(T))/p - 1 = 0
        sum(yi/(pi/p))  - 1 = 0
        sum(yi/(1 + dpdt(T-T0)/p)) - 1 = 0 
        special case: non-condensable:
        pi -> infinity p*yi/pi -> 0 
        =#
        function f0p(T)
            resinv = zero(T+first(y)+first(crit)[1])
            for i in 1:length(pure)
                if condensables[i] #if non-volatile, pi -> 0
                    Ti = sat[i][1]
                    pᵢ = p + dPdTsat[i]*(T - Ti)
                    pᵢyᵢ = y[i]/pᵢ
                    resinv += pᵢyᵢ
                end
            end
            return 1/resinv - p
        end
        T00 = sum(y[i]*first(sat[i]) for i in 1:length(model))
        fTd = Roots.ZeroProblem(f0p,T00)
        T0 = Roots.solve(fTd,Roots.Order0())
    else
        Tb = extrema(first,sat).*(0.9,1.1)
        f(T) = antoine_dew(pure,T,y,crit,condensables) - p
        fT = Roots.ZeroProblem(f,Tb)
        T0 = Roots.solve(fT,Roots.Order0())
    end
    _,vl0,vv0,x = __x0_dew_pressure(model,T0,y,nothing,condensables,pure,crit)
    return T0,vl0,vv0,x
end

function x0_dew_temperature(model::EoSModel,p,y)
    T0,V0_l,V0_v,x = __x0_dew_temperature(model,p,y)
    prepend!(x,log10.([V0_l,V0_v]))
    prepend!(x,T0)
    return x
end

function dew_temperature_init(model,p,y,vol0,T0,x0,condensables)
    if !isnothing(x0)
        if !isnothing(T0)
            if !isnothing(vol0)
                vl,vv = vol0
            else
                vl = volume(model,p,T0,x0,phase = :l)
                vv = volume(model,p,T0,y,phase =:v)
            end
        else
            T0,vl0,vv0,_ = __x0_dew_temperature(model,p,y,T0,condensables)
            if !isnothing(vol0)
                vl,vv = vol0
            else
                vl = min(vl0,volume(model,p,T0,x0,phase = :l))
                vv = max(vv0,volume(model,p,T0,y,phase =:v))
            end
        end
    else
        T00,vl0,vv0,x0 = __x0_dew_temperature(model,p,y,T0,condensables)
        if !isnothing(T0)
            vl = min(vl0,volume(model,p,T0,x0,phase = :l))
            vv = max(vv0,volume(model,p,T0,y,phase = :v))
        else
            vl = vl0
            vv = vv0
            T0 = T00
        end
    end
    return T0,vl,vv,x0
end

"""
    dew_temperature(model::EoSModel, p, y, method = ChemPotDewTemperature())

calculates the dew temperature and properties at a given pressure.
Returns a tuple, containing:
- Dew Temperature `[K]`
- liquid volume at Dew Point [`m³`]
- vapour volume at Dew Point [`m³`]
- Liquid composition at Dew Point

By default, uses equality of chemical potentials, via [`ChemPotDewTemperature`](@ref)
"""
function dew_temperature(model::EoSModel,p,x;kwargs...)
    if keys(kwargs) == (:v0,)
        nt_kwargs = NamedTuple(kwargs)
        v0 = nt_kwargs.v0
        T0 = v0[1]
        vl = exp10(v0[2])
        vv = exp10(v0[3])
        vol0 = (vl,vv)
        x0 = v0[4:end]
        _kwargs = (;T0,vol0,x0)
        method = init_preferred_method(dew_temperature,model,_kwargs)
    else
        method = init_preferred_method(dew_temperature,model,kwargs)
    end
    return dew_temperature(model,p,x,method)
end

function dew_temperature(model::EoSModel, p , x, T0::Number)
    kwargs = (;T0)
    method = init_preferred_method(dew_temperature,model,kwargs)
    return dew_temperature(model,p,x,method)
end

function dew_temperature(model::EoSModel,p,y,method::DewPointMethod)
    y = y/sum(y)
    p = float(p)
    model_r,idx_r = index_reduction(model,y)
    if length(model_r)==1
        (T_sat,v_l,v_v) = saturation_temperature(model_r,p)
        return (T_sat,v_l,v_v,y)
    end
    y_r = y[idx_r]
    (T_sat, v_l, v_v, x_r) = dew_temperature_impl(model_r,p,y_r,index_reduction(method,idx_r))
    x = index_expansion(x_r,idx_r)
    converged = bubbledew_check(v_l,v_v,y,x)
    if converged
        return (T_sat, v_l, v_v, x)
    else
        nan = zero(v_l)/zero(v_l)
        x = x*nan
        return (nan,nan,nan,x)
    end
end

include("dew_point/dew_chempot.jl")
include("dew_point/dew_fugacity.jl")
include("dew_point/dew_activity.jl")

function init_preferred_method(method::typeof(dew_pressure),model::EoSModel,kwargs)
    return ChemPotDewPressure(;kwargs...)
end

function init_preferred_method(method::typeof(dew_temperature),model::EoSModel,kwargs)
    return ChemPotDewTemperature(;kwargs...)
end

