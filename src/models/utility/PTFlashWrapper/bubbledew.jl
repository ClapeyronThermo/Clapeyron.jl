function __x0_bubble_pressure(model::PTFlashWrapper,T,x,y0 = nothing,volatiles = FillArrays.Fill(true,length(model)),pure = nothing,crit = nothing;verbose = false)
    sat = model.sat #saturation, we do not approximate here.
    p0r = first.(sat)
    p0 = index_expansion(p0r,volatiles)
    xipi = p0 .* x ./ sum(x)
    p0 = sum(xipi)
    if isnothing(y0)
        yx = xipi
        yx ./= p0
    else
        yx = y0
    end
    p,_,_,y,vl0,vv0 = improve_bubbledew_suggestion(model,p0,T,x,yx,FugEnum.BUBBLE_PRESSURE,volatiles,false)
    return p,vl0,vv0,y
end

function __x0_dew_pressure(model::PTFlashWrapper,T,y,x0=nothing,condensables = FillArrays.Fill(true,length(model)),pure = nothing, crit = nothing;verbose = false)
    sat = model.sat #saturation, we do not approximate here.
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
    p,_,x,_,vl0,vv0 = improve_bubbledew_suggestion(model,p0,T,xx,y,FugEnum.DEW_PRESSURE,condensables,false)
    return p,vl0,vv0,x
end

function __x0_bubble_temperature(model::PTFlashWrapper,p,x,Tx0 = nothing,volatiles = FillArrays.Fill(true,length(model)),pure = nothing,crit = nothing;verbose = false)
    x_r = @view x[volatiles]
    pure = @view model.pures[volatiles]
    sat = @view model.sat[volatiles]
    if Tx0 !== nothing
        T0 = Tx0
        for i in 1:length(pure)
            sat[i] = saturation_pressure(pure[i],T0)
        end
        p_i_r = first.(sat)
    else
        dPdTsat = extended_dpdT_temperature.(pure,p,crit)
        T0 = antoine_bubble_solve(dPdTsat,p,x_r)
        p_i_r = antoine_pressure.(dPdTsat,T0)
    end
    xipi_r = y_r = p_i_r .* x_r ./ sum(x_r)
    p = sum(xipi_r)
    y_r ./= p
    y0 = index_expansion(y_r,volatiles)
    _,T,_,y,vl0,vv0 = improve_bubbledew_suggestion(model,p,T0,x,y0,FugEnum.BUBBLE_TEMPERATURE,volatiles,false)
    update_temperature!(model,T)
    return T,vl0,vv0,y
end

function __x0_dew_temperature(model::PTFlashWrapper,p,y,Tx0 = nothing,condensables = FillArrays.Fill(true,length(model)),pure = split_pure_model(model,condensables),crit = nothing;verbose = false)
    y_r = @view y[condensables]
    pure = @view model.pures[condensables]
    sat = @view model.sat[condensables]
    if Tx0 !== nothing
        T0 = Tx0
        for i in 1:length(pure)
            sat[i] = saturation_pressure(pure[i],T0)
        end
        p0inv_r = 1 ./ first.(sat)
    else
        dPdTsat = extended_dpdT_temperature.(pure,p,crit)
        T0 = antoine_bubble_solve(dPdTsat,p,y_r)
        p0inv_r = 1 ./ antoine_pressure.(dPdTsat,T0)
    end
    yipi_r = x_r = y_r .* p0inv_r ./ sum(y_r)
    p_r = 1/sum(yipi_r)
    x_r .*= p_r
    x0 = index_expansion(x_r,condensables)
    update_temperature!(model,T0)
    _,T,x,_,vl0,vv0 = improve_bubbledew_suggestion(model,p,T0,x0,y,FugEnum.DEW_TEMPERATURE,condensables,false)
    return T,vl0,vv0,x
end

function improve_bubbledew_suggestion(model::PTFlashWrapper,p0,T0,x,y,method,in_media,high_conditions)
    TT = Base.promote_eltype(model,p0,T0,x,y)
    p,T = TT(p0),TT(T0)
    vl = volume(model,p,T,x,phase = :l)/sum(x)
    vv = volume(model,p,T,y,phase = :v)/sum(y)
    return p,T,x,y,vl,vv
end

function x0_edge_pressure(wrapper::PTFlashWrapper,T,z,pure = nothing)
  sat = wrapper.sat
  n = sum(z)
  p_bubble = sum(z[i]*first(sat[i]) for i in 1:length(sat))/n
  p_dew = n/sum(z[i]/first(sat[i]) for i in 1:length(sat))
  return (p_bubble,p_dew),sat
end

function _edge_pressure(wrapper::PTFlashWrapper,T,z,v0 = nothing,crit_retry = true)
    _1 = one(Base.promote_eltype(wrapper,T,z))
    if v0 == nothing
        p00 = _1
    else
        p00 = 0.5*(v0[1] + v0[2])*_1
    end
    sat = wrapper.sat
    RT = Rgas(wrapper)*T
    #=
    ∑zlogϕi,_ = ∑zlogϕ(gas_model(model),p,T,w,phase = :v)
    gl = excess_gibbs_free_energy(__γ_unwrap(model),p,T,w)/RT
    gv = ∑zlogϕi + tpd_delta_g_vapour(wrapper,p,T,w)
    f(T) = gl(T) - gv(T)


    system of eqs:
    variables:
    - vv
    - p

    gl - ∑zlogϕ(model,V,T,z) - tpd_delta_g_vapour(wrapper,p,T,w) = 0
    pressure(wrapper,vv,T,z) = p

    for ideal gas: solution is non-iterative
    for real gas: use ideal gas as starting point
    =#
    model = wrapper.model
    nc = length(model)
    gl = excess_gibbs_free_energy(__γ_unwrap(model),pmin,T,z)/RT #should be independent of pressure
    ∑z = sum(z)
    ∑zlogps = sum(z[i]*log(first(sat[i])) for i in 1:nc)

    p0 = exp((gl + ∑zlogps)/∑z)
    vv = ∑z*RT/p0
    gasmodel = gas_model(wrapper)

    nan = zero(p0)/zero(p0)
    fail = (nan,nan,nan)

    if gas_model(wrapper) isa IdealModel
        result = p0,volume(wrapper,p,T,z,phase = :l),vv
        return result,fail,:success
    end
    if v0 == nothing
        p = p0
    else
        p = v00
    end
    p = p0
    p_lb = minimum(first,sat)
    p_ub = maximum(first,sat)
    ∑zlogϕsat = zero(p)
    ∑zZl = zero(p)
    ∑zvlRT = zero(p)
    lnϕsat = wrapper.fug
    for i in 1:nc
        psi,vli,_ = sat[i]
        zi = z[i]
        ∑zlogϕsat += zi*lnϕsat[i]
        Zli = vli*psi/RT
        ∑zZl += zi*Zli
        ∑zvlRT += zi*vli/RT
    end

    for i in 1:40
        vv_old = vv
        vv = volume(gasmodel,p,T,z,phase = :v,vol0 = vv)
        ∑zlogϕi,_ = ∑zlogϕ(gasmodel,p,T,z,phase = :v,vol = vv)
        p_old = p
        #p*vl/RT - vl*ps/RT + log(ϕsat[i]) + log(ps)
        p = exp((gl + p*∑zvlRT - ∑zZl + ∑zlogϕsat + ∑zlogps - ∑zlogϕi)/∑z)
        p < p_lb && (p = 0.5*(p_old + p_lb))
        p > p_ub && (p = 0.5*(p_old + p_ub))
        if abs(p - p_old)/p < sqrt(eps(eltype(p)))
            vl = volume(wrapper,p,T,z,phase = :l)
            return (p,vl,vv),fail,:success
        end
    end

    return fail,fail,:failure
end

function x0_edge_temperature(wrapper::PTFlashWrapper,p,z,pure = wrapper.pures)
    dPdTsat = extended_dpdT_temperature.(pure,p)
    T_bubble = antoine_bubble_solve(dPdTsat,p,z)
    T_dew = antoine_dew_solve(dPdTsat,p,z)
    return (T_bubble,T_dew),dPdTsat
end

function _edge_temperature(model::PTFlashWrapper,p,z,v0 = nothing)
    if v0 == nothing
        vv0,_ = x0_edge_temperature(model,p,z)
    else
        vv0 = (v0[1],v0[2])
    end

    Tmin,Tmax = minmax(vv0[1],vv0[2])
    tau0 = 0.5(1/Tmin + 1/Tmax)
    T0 = 1/tau0
    _0 = zero(Base.promote_eltype(model,T0,z))
    nan = oftype(_0,NaN)
    fail = (nan,nan,nan)

    γmodel = __γ_unwrap(model)
    gasmodel = gas_model(model)

    R = Rgas(model)
    n = sum(z)
    vv = Ref(n*R*T0/p)
   
    function obj(tau)
        T = 1/tau
        update_temperature!(model,T)
        vv0 = max(n*R*T/p,vv[])
        vvi = volume(gasmodel,p,T,z,phase = :v,vol0 = vv0)
        ∑zlogϕi,_ = ∑zlogϕ(gasmodel,p,T,z,phase = :v,vol = vvi)
        gl = excess_gibbs_free_energy(__γ_unwrap(model),p,T,z)/(R*T)
        gv = ∑zlogϕi + tpd_delta_g_vapour(model,p,T,z)
        vv[] = vvi
        return gl-gv
    end

    prob = Roots.ZeroProblem(obj,tau0)
    tau_edge = Roots.solve(prob,Roots.Order0())
    T_edge = 1/tau_edge

    if isfinite(T_edge)
        vv_edge = vv[]
        vl_edge = volume(model,p,T_edge,z,phase = :l)
        return (T_edge,vl_edge,vv_edge),fail,:success
    end

    return fail,fail,:failure
end

function init_preferred_method(method::typeof(bubble_pressure),model::PTFlashWrapper,kwargs)
    return ActivityBubblePressure(;kwargs...)
end

function init_preferred_method(method::typeof(bubble_temperature),model::PTFlashWrapper,kwargs)
    return FugBubbleTemperature(;kwargs...)
end

function init_preferred_method(method::typeof(dew_pressure),model::PTFlashWrapper,kwargs)
    return ActivityDewPressure(;kwargs...)
end

function init_preferred_method(method::typeof(dew_temperature),model::PTFlashWrapper,kwargs)
    return FugDewTemperature(;kwargs...)
end


function bubble_pressure_pproperty_method(model::PTFlashWrapper,p0,T,z,sat)
  y0 = z .* first.(sat)
  y0 ./= sum(y0)
  p,_,_,y,vl0,vv0 = improve_bubbledew_suggestion(model,p0,T,z,y0,FugEnum.BUBBLE_PRESSURE,FillArrays.Trues(length(z)),false)
  return FugBubblePressure(vol0 = (vl0,vv0),p0 = p,y0 = y)
end

function dew_pressure_pproperty_method(model::PTFlashWrapper,p0,T,z,sat)
  x0 = z ./ first.(sat)
  x0 ./= sum(x0)
  p,_,x,_,vl0,vv0 = improve_bubbledew_suggestion(model,p0,T,x0,z,FugEnum.DEW_PRESSURE,FillArrays.Trues(length(z)),false)
  return FugDewPressure(vol0 = (vl0,vv0),p0 = p,x0 = x)
end

function bubble_temperature_tproperty_method(model::PTFlashWrapper,p,T0,z,dPdT)
  y0 = z .* antoine_pressure.(dPdT,T0)
  y0 ./= sum(y0)
  _,T,_,y,vl0,vv0 = improve_bubbledew_suggestion(model,p,T0,z,y0,FugEnum.BUBBLE_TEMPERATURE,FillArrays.Trues(length(z)),false)
  return FugBubbleTemperature(vol0 = (vl0,vv0),T0 = T,y0 = y)
end

function dew_temperature_tproperty_method(model::PTFlashWrapper,p,T0,z,dPdT)
  x0 = z ./ antoine_pressure.(dPdT,T0)
  x0 ./= sum(x0)
  _,T,x,_,vl0,vv0 = improve_bubbledew_suggestion(model,p,T0,x0,z,FugEnum.DEW_TEMPERATURE,FillArrays.Trues(length(z)),false)
  return FugDewTemperature(vol0 = (vl0,vv0),T0 = T,x0 = x)
end
