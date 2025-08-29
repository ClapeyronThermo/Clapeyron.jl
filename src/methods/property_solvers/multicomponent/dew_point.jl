
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

function __x0_dew_pressure(model::EoSModel,T,y,x0=nothing,condensables = FillArrays.Fill(true,length(model)),pure = split_pure_model(model,condensables), crit = nothing)
    sat = extended_saturation_pressure.(pure,T,crit) #saturation, or approximation via critical point.
    p0inv_r = 1. ./ first.(sat)
    p0inv = index_expansion(p0inv_r,condensables)
    yipi = y .* p0inv ./ sum(y)
    p0 = 1/sum(yipi)
    if isnothing(x0)
        xx = yipi
        xx .*= p0
    else
        xx = x0
    end
    high_conditions = __is_high_pressure_state(pure,sat,T)
    p,_,x,_,vl0,vv0 = improve_bubbledew_suggestion(model,p0,T,xx,y,FugEnum.DEW_PRESSURE,condensables,high_conditions)
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
    dew_pressure(model::EoSModel, T, y, method = ChemPotDewPressure())

Calculates the dew pressure and properties at a given temperature `T`.
Returns a tuple, containing:
- Dew Pressure `[Pa]`
- Liquid volume at Dew Point `[m³]`
- Vapour volume at Dew Point `[m³]`
- Liquid composition at Dew Point

By default, uses equality of chemical potentials, via [`ChemPotDewPressure`](@ref)
"""
function dew_pressure(model::EoSModel,T,x;kwargs...)
    moles_positivity(x)
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

function dew_pressure(model::EoSModel, T, y, method::ThermodynamicMethod)
    moles_positivity(y)
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
    converged = bubbledew_check(model,P_sat,T,v_l,v_v,x,y)
    if converged
        return (P_sat, v_l, v_v, x)
    else
        nan = zero(v_l)/zero(v_l)
        x = x*nan
        return (nan,nan,nan,x)
    end
end



function __x0_dew_temperature(model::EoSModel,p,y,Tx0 = nothing,condensables = FillArrays.Fill(true,length(model)),pure = split_pure_model(model,condensables),crit = nothing)
    y_r = @view y[condensables]

    if Tx0 !== nothing
        T0 = Tx0
        sat = extended_saturation_pressure.(pure,T0,crit)
        p0inv_r = 1.0 ./ first.(sat)
        high_conditions = __is_high_pressure_state(pure,sat,T0)
    else
        dPdTsat = extended_dpdT_temperature.(pure,p,crit)
        T0 = antoine_dew_solve(dPdTsat,p,y_r)
        p0inv_r = 1.0 ./ antoine_pressure.(dPdTsat,T0)
        high_conditions = __is_high_temperature_state(pure,dPdTsat,T0)
    end
    yipi_r = x_r = y_r .* p0inv_r ./ sum(y_r)
    p = 1/sum(yipi_r)
    x_r .*= p
    x0 = index_expansion(x_r,condensables)
    _,T,x,_,vl0,vv0 = improve_bubbledew_suggestion(model,p,T0,x0,y,FugEnum.DEW_TEMPERATURE,condensables,high_conditions)
    return T,vl0,vv0,x
end

function antoine_dew_solve(dpdt,p_dew,y)
    if length(dpdt) == 1
        #sum(y)/pinv - p_dew
        #y = 1.0, pinv  = 1/p(T)
        #p(T) = p_dew
        dlnpdTinv,logp0,T0inv = dpdt[1]
        Tinv = (log(p_dew) - logp0)/dlnpdTinv + T0inv
        return 1/Tinv
    end

    function antoine_f0(T)
        pinv = zero(T+first(y)+first(dpdt)[1])
        for i in 1:length(dpdt)
            pᵢ = antoine_pressure(dpdt[i],T)
            pᵢyᵢ = y[i]/pᵢ
            pinv += pᵢyᵢ
        end
        return sum(y)/pinv - p_dew
    end
    Tmin,Tmax = extrema(x -> 1/last(x),dpdt)
    if antoine_f0(Tmin)*antoine_f0(Tmax) < 0
        prob = Roots.ZeroProblem(antoine_f0,(Tmin,Tmax))
        return Roots.solve(prob)
    else
        prob = Roots.ZeroProblem(antoine_f0,0.5*(Tmin+Tmax))
        return Roots.solve(prob) 
    end
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

Calculates the dew temperature and properties at a given pressure `p`.
Returns a tuple, containing:
- Dew Temperature `[K]`
- Liquid volume at Dew Point `[m³]`
- Vapour volume at Dew Point `[m³]`
- Liquid composition at Dew Point

By default, uses equality of chemical potentials, via [`ChemPotDewTemperature`](@ref)
"""
function dew_temperature(model::EoSModel,p,x;kwargs...)
    moles_positivity(x)
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
    moles_positivity(x)
    kwargs = (;T0)
    method = init_preferred_method(dew_temperature,model,kwargs)
    return dew_temperature(model,p,x,method)
end

function dew_temperature(model::EoSModel,p,y,method::ThermodynamicMethod)
    moles_positivity(y)
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
    converged = bubbledew_check(model,p,T_sat,v_l,v_v,x,y)
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
