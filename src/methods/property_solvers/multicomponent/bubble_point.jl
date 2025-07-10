
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


function extended_saturation_pressure(pure,T,_crit = nothing; crit_retry = true)
    sat,crit,status = _extended_saturation_pressure(pure,T,_crit;crit_retry)
    if status == :supercritical
        Tc,Pc,Vc = crit
        #create initial point from critical values
        #we use a pseudo-saturation pressure extension,based on the slope at the critical point.
        dlnpdTinv,logp0,Tcinv = __dlnPdTinvsat(pure,sat,crit,T,false,:supercritical)
        lnp = logp0 + dlnpdTinv*(1/T - Tcinv)
        p0 = exp(lnp)
        vl0 = x0_volume(pure,p0,T,phase = :l)
        vv0 = max(1.2*Vc,3*Rgas(pure)*T/Pc)
        return (p0,vl0,vv0)
    else
        return sat
    end

end

function _extended_saturation_pressure(pure, T, _crit = nothing; crit_retry = true)
    #try without critical point information
    _0 = zero(Base.promote_eltype(pure,T))
    nan = _0/_0
    #no crit point available, try calculating sat_p without it
    if _crit === nothing
        sat = saturation_pressure(pure,T,crit_retry = false)
        if isnan(first(sat))
            if !crit_retry
                return sat,(nan,nan,nan),:fail
            end
        else
            return sat,(nan,nan,nan),:success
        end
    end
    #calculate critical point, try again
    if _crit !== nothing && !isnan(first(_crit))
        crit = _crit
    elseif crit_retry
        crit = crit_pure(pure)
    else
        nan = _0/_0
        crit = (nan,nan,nan)
    end
    Tc,Pc,Vc = crit
    if T < Tc
        sat = saturation_pressure(pure,T,crit = crit) #calculate sat_p with crit info
        !isnan(first(sat)) && (return sat,crit,:success)
    else
        nan = _0/_0
        sat = (nan,nan,nan)
        return sat, crit, :supercritical
    end
end

function extended_saturation_temperature(pure,p,_crit = nothing; crit_retry = true)
    sat,crit,status = _extended_saturation_temperature(pure,p,_crit;crit_retry)
    if status == :supercritical
        #create initial point from critical values
        #we use a pseudo-saturation pressure extension,based on the slope at the critical point.
        Tc,Pc,Vc = crit
        dlnpdTinv,logp0,Tcinv = __dlnPdTinvsat(pure,sat,crit,p,true,:supercritical)
        #lnp = logp0 + dlnpdTinv*(1/T - Tcinv)
        Tinv = (log(p) - logp0)/dlnpdTinv + Tcinv
        T0  = 1/Tinv
        vl0 = x0_volume(pure,p,T0,phase = :l)
        vv0 = max(1.2*Vc,3*Rgas(pure)*T0/Pc)
        return (T0,vl0,vv0)
    else
        return sat
    end
end

#this function does not do the crit calculation.
function _extended_saturation_temperature(pure, p, _crit = nothing; crit_retry = true)
    _0 = zero(Base.promote_eltype(pure,p))
    nan = _0/_0
    if _crit === nothing #no crit point available, try calculating sat_p without it
        sat = saturation_temperature(pure,p,crit_retry = false)
        if isnan(first(sat))
            if !crit_retry
                return sat,(nan,nan,nan),:fail #failed
            end
        else
            return sat,(nan,nan,nan),:success #sucess
        end
    end

    #calculate critical point, try again
    if _crit !== nothing
        crit = _crit
    elseif crit_retry
        crit = crit_pure(pure)
    else
        nan = _0/_0
        crit = (nan,nan,nan)
    end

    Tc,Pc,Vc = crit
    if p < Pc
        sat = saturation_temperature(pure,p,crit = crit) #calculate sat_p with crit info
        !isnan(first(sat)) && (return sat,crit,:success)
    else
        nan = _0/_0
        sat = (nan,nan,nan)
        (return sat,crit,:supercritical)
    end
    #create initial point from critical values
    #we use a pseudo-saturation pressure extension,based on the slope at the critical point.

    dlnpdTinv,logp0,Tcinv = __dlnPdTinvsat(pure,sat,crit,p,true,:supercritical)
    #lnp = logp0 + dlnpdTinv*(1/T - Tcinv)
    Tinv = (log(p) - logp0)/dlnpdTinv + Tcinv
    T0  = 1/Tinv
    vl0 = x0_volume(pure,p,T0,phase = :l)
    vv0 = max(1.2*Vc,3*Rgas(pure)*T0/Pc)
    return (T0,vl0,vv0),crit,:supercritical
end

function __is_high_pressure_state(pure,sat,T)
    p,vl,vv = sat
    B = second_virial_coefficient(pure,T)
    -2B > vv || p > -0.25*Rgas(pure)*T/B
end

function __is_high_pressure_state(pure::AbstractVector,sat,T)
    for i in 1:length(pure)
        __is_high_pressure_state(pure[i],sat[i],T) && return true
    end
    return false
end

function __is_high_temperature_state(pure,dpdT,T)
    p = antoine_pressure(dpdT,T)
    vv = volume(pure,p,T,phase = :v)
    B = second_virial_coefficient(pure,T)
    -2B > vv || p > -0.25*Rgas(pure)*T/B
end

function __is_high_temperature_state(pure::AbstractVector,dpdT,T)
    for i in 1:length(pure)
        __is_high_temperature_state(pure[i],dpdT[i],T) && return true
    end
    return false
end

function __crit_pure(sat0,pure)
    if isnan(first(sat0))
        return crit_pure(pure)
    else
        p,vl,vv = sat0
        nan = zero(p)/zero(p)
        return nan,nan,nan
    end
end

function __dlnPdTinvsat(pure,sat,crit,xx,is_sat_temperature,status)
    successful_saturation = status == :success
    yy,vl,vv = sat
    if is_sat_temperature
        p,T = xx,yy
    else
        p,T = yy,xx
    end

    if status == :success
        dpdT = dpdT_saturation(pure,vl,vv,T)
        return -dpdT*T*T/p,log(p),1/T
    elseif status === :supercritical
        Tc,Pc,Vc = crit
        _p(_T) = pressure(pure,Vc,_T)
        dpdT = Solvers.derivative(_p,Tc)
        return -dpdT*Tc*Tc/Pc,log(Pc),1/Tc
    elseif status == :fail
        return sat
    else
        throw(error("dPdTsat: invalid status: $status"))
    end
end

function extended_dpdT_pressure(pure,T,crit = nothing)
    sat,_crit,status = _extended_saturation_pressure(pure,T,crit)
    return __dlnPdTinvsat(pure,sat,_crit,T,false,status)
end

function extended_dpdT_temperature(pure,p,crit = nothing)
    sat,_crit,status = _extended_saturation_temperature(pure,p,crit)
    return  __dlnPdTinvsat(pure,sat,_crit,p,true,status)
end

function improve_bubbledew_suggestion_spinodal(model,p0,T0,x,y,method,in_media)
    #TODO: implement this
    return p0,T0
    #=
    if FugEnum.is_bubble(method)
        z = x
        w = y
        phasez,phasez0 = :liquid,:vapour
    else
        z = y
        w = x
        phasez,phasez0 = :vapour,:liquid
    end

    z_0 = similar(z)
    z_0 .= z
    zero_non_equilibria!(z_0,in_media)

    if FugEnum.is_pressure(method) && FugEnum.is_bubble(method)
        pmid,vmid,_ = eigmin_minimum_pressure(model,T0,x,volume(model,p0,T0,x,phase = :l))
        if p0 < 0

        end

    else
        return p0,T0
    end
    =#

end

function improve_bubbledew_suggestion(model,p0,T0,x,y,method,in_media,high_conditions)
    if high_conditions #inprove p/T via spinodal
        p,T = improve_bubbledew_suggestion_spinodal(model,p0,T0,x,y,method,in_media)
    else
        p,T = p0,T0
    end

    vlx = volume(model,p,T,x,phase = :l)
    μl = VT_chemical_potential_res(model,vlx,T,x)
    RT = Rgas(model) * T
    Zl = p*vlx/RT/sum(x)
    ϕl = K = similar(μl)
    ϕl .= exp.(μl ./ RT) ./ Zl
    ϕv = virial_phi(model,p,T,y) #virial fugacity coefficient, skips volume calculation
    if all(!isnan,@view(ϕv[in_media]))
        K .= ϕl ./ ϕv
    end
    K_r = @view K[in_media]
    if FugEnum.is_bubble(method)
        x_r = @view x[in_media]
        y_r = rr_flash_vapor(K_r,x_r,zero(eltype(K)))
        yy = index_expansion(y_r,in_media)
        yy ./= sum(yy)
        vv = volume(model,p,T,y,phase = :v)
        return p,T,x,yy,vlx,vv
    else
        y_r = @view y[in_media]
        x_r = rr_flash_liquid(K_r,y_r,one(eltype(K)))
        xx = index_expansion(x_r,in_media)
        xx ./= sum(xx)
        vl = volume(model,p,T,xx,phase = :l)
        vv = volume(model,p,T,y,phase = :v)
        return p,T,xx,y,vl,vv
    end
end

_virial(model,V,T,z) = second_virial_coefficient(model,T,z)

function virial_phi(model,p,T,z)
    pRT = p/(Rgas(model)*T)
    dB = VT_molar_gradient(model,zero(p),T,z,_virial)
    return exp.(dB .* pRT)
end

function __x0_bubble_pressure(model::EoSModel,T,x,y0 = nothing,volatiles = FillArrays.Fill(true,length(model)),pure = split_pure_model(model,volatiles),crit = nothing)
    #check each T with T_scale, if treshold is over, replace Pi with inf
    sat = extended_saturation_pressure.(pure,T,crit) #saturation, or aproximation via critical point.
    p0r = first.(sat)
    p0 = index_expansion(p0r,volatiles)
    xipi = p0 .* x
    p0 = sum(xipi)
    if isnothing(y0)
        yx = xipi
        yx ./= p0
    else
        yx = y0
    end

    high_conditions = __is_high_pressure_state(pure,sat,T)
    p,_,_,y,vl0,vv0 = improve_bubbledew_suggestion(model,p0,T,x,yx,FugEnum.BUBBLE_PRESSURE,volatiles,high_conditions)
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

function bubble_pressure(model::EoSModel, T, x, method::ThermodynamicMethod)
    x = x/sum(x)
    T = float(T)
    model_r,idx_r = index_reduction(model,x)
    if length(model_r)==1
        (P_sat,v_l,v_v) = saturation_pressure(model_r,T)
        return (P_sat,v_l,v_v,x)
    end
    x_r = x[idx_r]
    if has_a_res(model)
        bubble_pressure_result_primal = bubble_pressure_impl(primalval(model_r),primalval(T),primalval(x_r),index_reduction(method,idx_r))
        bubble_pressure_result = bubble_pressure_ad(model_r,T,x_r,bubble_pressure_result_primal)
    else
        bubble_pressure_result = bubble_pressure_impl(model_r,T,x_r,index_reduction(method,idx_r))
    end
    (P_sat, v_l, v_v, y_r) = bubble_pressure_result
    y = index_expansion(y_r,idx_r)
    converged = bubbledew_check(model,P_sat,T,v_v,v_l,y,x)
    if converged
        return (P_sat, v_l, v_v, y)
    else
        nan = zero(v_l)/zero(v_l)
        y = y*nan
        return (nan,nan,nan,y)
    end
end

###Bubble Temperature

function __x0_bubble_temperature(model::EoSModel,p,x,Tx0 = nothing,volatiles = FillArrays.Fill(true,length(model)),pure = split_pure_model(model,volatiles),crit = nothing)
    x_r = @view x[volatiles]

    if Tx0 !== nothing
        T0 = Tx0
        sat = extended_saturation_pressure.(pure,T0,crit)
        p_i_r = first.(sat)
        high_conditions = __is_high_pressure_state(pure,sat,T0)
    else
        dPdTsat = extended_dpdT_temperature.(pure,p,crit)
        T0 = antoine_bubble_solve(dPdTsat,p,x_r)
        p_i_r = antoine_pressure.(dPdTsat,T0)
        high_conditions = __is_high_temperature_state(pure,dPdTsat,T0)
    end
    xipi_r = y_r = p_i_r .* x_r
    p = sum(xipi_r)
    y_r ./= p
    y0 = index_expansion(y_r,volatiles)
    _,T,_,y,vl0,vv0 = improve_bubbledew_suggestion(model,p,T0,x,y0,FugEnum.BUBBLE_TEMPERATURE,volatiles,high_conditions)
    return T,vl0,vv0,y
end

function antoine_pressure(dpdT,T)
    dlnpdTinv,logp0,T0inv = dpdT
    return exp(logp0 + dlnpdTinv*(1/T - T0inv))
end

function antoine_bubble_solve(dpdt,p_bubble,x,T0 = nothing)

    if length(dpdt) == 1
        #p(T) = p_bubble = exp(logp0 + dlnpdTinv*(1/T - T0inv))
        dlnpdTinv,logp0,T0inv = dpdt[1]
        Tinv = (log(p_bubble) - logp0)/dlnpdTinv + T0inv
        return 1/Tinv
    end

    function antoine_f0(T)
        p = zero(T+first(x)+first(dpdt)[1])
        for i in 1:length(dpdt)
            pᵢ = antoine_pressure(dpdt[i],T)
            pᵢxᵢ = x[i]*pᵢ
            p += pᵢxᵢ
        end
        return p/sum(x) - p_bubble
    end


    if T0 === nothing
    Tmin,Tmax = extrema(x -> 1/last(x),dpdt)
        prob = Roots.ZeroProblem(antoine_f0,(Tmin,Tmax))
        return Roots.solve(prob)
    else
        return Roots.ZeroProblem(antoine_f0,T0)
        return Roots.solve(prob)
    end
end

function x0_bubble_temperature(model::EoSModel,p,x)
    T0,V0_l,V0_v,y = __x0_bubble_temperature(model,p,x)
    v0 = similar(x)
    return vcat(T0, log10(V0_l),log10(V0_v),v0)
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
- liquid volume at Bubble Point `[m³]`
- vapour volume at Bubble Point `[m³]`
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

function bubble_temperature(model::EoSModel, p, x, method::ThermodynamicMethod)
    x = x/sum(x)
    p = float(p)
    model_r,idx_r = index_reduction(model,x)
    if length(model_r)==1
        (T_sat,v_l,v_v) = saturation_temperature(model_r,p)
        return (T_sat,v_l,v_v,x)
    end
    x_r = x[idx_r]


    if has_a_res(model)
        bubble_temperature_result_primal =  bubble_temperature_impl(primalval(model_r),primalval(p),primalval(x_r),index_reduction(method,idx_r))
        bubble_temperature_result =  bubble_temperature_ad(model_r,p,x_r,bubble_temperature_result_primal)
    else
        bubble_temperature_result =  bubble_temperature_impl(model_r,p,x_r,index_reduction(method,idx_r))
    end

    (T_sat, v_l, v_v, y_r) = bubble_temperature_result
    y = index_expansion(y_r,idx_r)
    converged = bubbledew_check(model,p,T_sat,v_v,v_l,y,x)
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
