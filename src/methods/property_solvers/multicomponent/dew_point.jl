
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
    pure_vals = extended_saturation_pressure.(pure,T,crit,condensables,false) #saturation, or aproximation via critical point.
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

function dew_pressure(model::EoSModel, T, y,method::ThermodynamicMethod)
    y = y/sum(y)
    T = float(T)
    model_r,idx_r = index_reduction(model,y)
    if length(model_r)==1
        (P_sat,v_l,v_v) = saturation_pressure(model_r,T)
        return (P_sat,v_l,v_v,y)
    end
    y_r = y[idx_r]

    if has_a_res(model)
        dew_pressure_result_primal = dew_pressure_impl(primalval(model_r),primalval(T),primalval(y_r),index_reduction(method,idx_r))
        dew_pressure_result = dew_pressure_ad(model_r,T,y_r,dew_pressure_result_primal)
    else
        dew_pressure_result = dew_pressure_impl(model_r,T,y_r,index_reduction(method,idx_r))
    end

    (P_sat, v_l, v_v, x_r) = dew_pressure_result
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



function __x0_dew_temperature(model::EoSModel,p,y,Tx0 = nothing,condensables = FillArrays.Fill(true,length(model)),pure = split_model(model),crit = nothing)
    multi_component_check(x0_dew_temperature,model)
        
    if Tx0 !== nothing
        _crit = isnothing(crit) ?  FillArrays.fill(nothing,length(model)) : crit
        K = suggest_K(model,p,Tx0,y,pure,condensables,_crit)
        x = rr_flash_liquid(K,y,one(eltype(K)))
        zero_non_equilibria!(x,condensables)
        x ./= sum(x)
        vl0 = volume(model,p,Tx0,x,phase = :l)
        vv0 = volume(model,p,Tx0,y,phase = :v)
        #this is exactly like __x0_bubble_pressure, but we use T0, instead of an input T
        #_,vl0,vv0,y = __x0_bubble_pressure(model,T0,x,nothing,volatiles,pure,crit)
        return Tx0,vl0,vv0,x
    end
    
    
    
    sat = extended_saturation_temperature.(pure,p,crit,condensables)
    if crit === nothing
        _crit = __crit_pure.(sat,pure,condensables)
    else
        _crit = crit
    end
    fix_sat_ti!(sat,pure,_crit,p,condensables)
    if Tx0 !== nothing
        T0 = Tx0
    else
        dPdTsat = __dlnPdTinvsat.(pure,sat,_crit,p)
        prob = antoine_dew_problem(dPdTsat,p,y)
        T0 = Roots.solve(prob)
    end
    K = suggest_K(model,p,T0,y,pure,FillArrays.fill(true,length(model)),_crit)
    x = rr_flash_liquid(K,y,one(eltype(K)))
    zero_non_equilibria!(x,condensables)
    x ./= sum(x)
    vl0 = volume(model,p,T0,x,phase = :l)
    vv0 = volume(model,p,T0,y,phase = :v)
    #_,vl0,vv0,x = __x0_dew_pressure(model,T0,y,nothing,condensables,pure,crit)
    return T0,vl0,vv0,x
end

function antoine_dew_problem(dpdt,p_dew,y,condensables = FillArrays.Fill(true,length(dpdt)))  
    function antoine_f0(T)
        pinv = zero(T+first(y)+first(dpdt)[1])
        for i in 1:length(dpdt)
            dlnpdTinv,logp0,T0inv = dpdt[i]
            if condensables[i]
                pᵢ = exp(logp0 + dlnpdTinv*(1/T - T0inv))
                pᵢyᵢ = y[i]/pᵢ
                pinv += pᵢyᵢ
            end
        end
        return sum(y)/pinv - p_dew
    end
    Tmin,Tmax = extrema(x -> 1/last(x),dpdt)
    return Roots.ZeroProblem(antoine_f0,(Tmin,Tmax))
end

function x0_dew_temperature(model::EoSModel,p,y,T0 = nothing)
    T0,V0_l,V0_v,x = __x0_dew_temperature(model,p,y,T0)
    v0 = similar(x)
    v0 .= x
    return vcat(T0,log10(V0_l),log10(V0_v),v0)
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

function dew_temperature(model::EoSModel,p,y,method::ThermodynamicMethod)
    y = y/sum(y)
    p = float(p)
    model_r,idx_r = index_reduction(model,y)
    if length(model_r)==1
        (T_sat,v_l,v_v) = saturation_temperature(model_r,p)
        return (T_sat,v_l,v_v,y)
    end
    y_r = y[idx_r]
    
    if has_a_res(model)
        dew_temperature_result_primal =  dew_temperature_impl(primalval(model_r),primalval(p),primalval(y_r),index_reduction(method,idx_r))
        dew_temperature_result =  dew_temperature_ad(model_r,p,y_r,dew_temperature_result_primal)
    else
        dew_temperature_result =  dew_temperature_impl(model_r,p,y_r,index_reduction(method,idx_r))
    end

    (T_sat, v_l, v_v, x_r) = dew_temperature_result
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
