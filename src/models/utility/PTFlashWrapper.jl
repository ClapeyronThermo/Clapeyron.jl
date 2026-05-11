#wrapper used to cache results in case of activity models and CompositeModel
#=
struct PTFlashWrapper{T,T2,R,S} <: EoSModel
    components::Vector{String}
    model::T
    pures::T2
    sat::Vector{R}
    fug::Vector{S}
    equilibrium::Symbol
end =#

Base.length(Base.@specialize(model::PTFlashWrapper)) = length(model.model)
Base.eltype(Base.@specialize(model::PTFlashWrapper)) = Base.promote_eltype(model.model,model.fug)
__γ_unwrap(Base.@specialize(model::PTFlashWrapper)) = __γ_unwrap(model.model)
@inline gas_model(Base.@specialize(model::PTFlashWrapper)) = gas_model(model.model)

function tp_flash_K0(wrapper::PTFlashWrapper,p,T,z)
    first.(wrapper.sat) ./ p
end

function tp_flash_K0!(K,wrapper::PTFlashWrapper,p,T,z)
    K .=  first.(wrapper.sat) ./ p 
end

function PTFlashWrapper{TT}(model,equilibrium,pures = split_pure_model(model)) where TT
    nc = length(model)
    sat = Vector{Tuple{TT,TT,TT}}(undef,nc)
    ϕpure = Vector{TT}(undef,nc)
    return PTFlashWrapper(component_list(model),model,pures,sat,ϕpure,equilibrium)
end

function PTFlashWrapper(model,equilibrium,pures = split_pure_model(model))
    TT = Base.promote_eltype(model,Float64)
    return PTFlashWrapper{TT}(model,equilibrium,pures)
end

function PTFlashWrapper(model,p,T,z,equilibrium)
    TT = Base.promote_eltype(model,p,T,z)
    wrapper = PTFlashWrapper{TT}(model,equilibrium)
    update_temperature!(wrapper,T)
    return wrapper
end

split_pure_model(model::PTFlashWrapper,splitter) = model.pures[splitter]
idealmodel(model::PTFlashWrapper) = idealmodel(model.model)
fluid_model(model::PTFlashWrapper) = fluid_model(model.model)
molecular_weight(model::PTFlashWrapper,z) = molecular_weight(model.model,z)
reference_state(model::PTFlashWrapper) = reference_state(model.model)

function update_temperature!(model::PTFlashWrapper,T)
    isnan(T) && return nothing
    pures = model.pures
    lnϕ = model.fug
    sats = model.sat
    TT = eltype(lnϕ)
    for i in 1:length(model)
        pure = pures[i]
        sat = saturation_pressure(pure,T)
        ps,vl,vv = sat
        sats[i] = sat
        if gas_model(pure) isa IdealModel
            lnϕ[i] = 0.0
        else
            lnϕ[i] = VT_lnϕ_pure(gas_model(pure),vv,T,ps)
        end
    end
    return nothing
end


function _update_temperature_with_view!(model1::TT,model2::TT,T,_view) where TT <: PTFlashWrapper
    n1,n2 = length(model1),length(model2)
    if n1 > n2
        update_temperature!(model1,T)
        model2.sat .= @view model1.sat[_view]
        model2.fug .= @view model1.fug[_view]
    else
        update_temperature!(model2,T)
        model1.sat .= @view model2.sat[_view]
        model1.fug .= @view model2.fug[_view]
    end
    return nothing
end

_update_temperature_with_view!(model1,model2,T,_view) = nothing

function volume_impl(model::PTFlashWrapper, p, T, z, phase, threaded, vol0)
    if is_unknown(phase) || phase == :stable
        new_phase = identify_phase(model, p, T, z)
        return volume_impl(model.model, p, T, z, new_phase, threaded, vol0)
    else
        volume_impl(model.model, p, T, z, phase, threaded, vol0)
    end
end

function Base.show(io::IO,mime::MIME"text/plain",wrapper::PTFlashWrapper)
    model = wrapper.model
    pure = wrapper.pures
    print(io,"PT-Flash Wrapper")
    length(model) == 1 && print(io, " with 1 component:")
    length(model) > 1 && print(io, " with ", length(model), " components:")
    println(io)
    show_pairs(io,wrapper.components)
    print(io,'\n',"Mixture model: ", typeof(model))
    print(io,'\n',"Pure model: ",eltype(pure))
    print(io,'\n',"Equilibrium type: :",wrapper.equilibrium)
    show_reference_state(io,model;space = true)
end

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

function modified_lnϕ(wrapper::PTFlashWrapper, p, T, z, cache; phase = :unknown, vol0 = nothing)
    if is_vapour(phase) || is_liquid(phase)
        lnϕz,vz = tpd_lnϕ_and_v!(cache,wrapper,p,T,z,vol0,false,phase,nothing)
        return lnϕz,vz
    elseif is_unknown(phase)
        lnϕz1,vzl = tpd_lnϕ_and_v!(cache,wrapper,p,T,z,vol0,false,:liquid,nothing)
        lnϕzl = copy(lnϕz1)
        logsumz = log(sum(z))
        minz = -1e100*one(eltype(z))
        lnϕz1 .+ log.(z) .- logsumz
        gl =  @sum(lnϕz1[i]*max(z[i],minz))
        lnϕz2,vzv = tpd_lnϕ_and_v!(cache,wrapper,p,T,z,vol0,false,:vapour,nothing)
        lnϕzv = copy(lnϕz2)
        lnϕz2 .+ log.(z) .- logsumz
        gv = @sum(lnϕz2[i]*max(z[i],minz))
        if gv < gl
            return lnϕzv,vzv
        else
            return lnϕzl,vzl
        end
    else
        throw(error("invalid phase specification, got $phase"))
    end
end

function K0_lle_init(wrapper::PTFlashWrapper,p,T,z)
    return K0_lle_init(__γ_unwrap(wrapper),p,T,z)
end

function modified_gibbs(wrapper::PTFlashWrapper,p,T,w,phase,vol)
    model = wrapper.model
    TT = Base.promote_eltype(wrapper,p,T,w)
    RT = Rgas(model)*T
    ∑w = sum(w)
    iszero(∑w) && return zero(TT), zero(TT)
    g_ideal = sum(xlogx,w) - xlogx(∑w)
    vl = zero(TT)
    if is_liquid(phase)
        return excess_gibbs_free_energy(__γ_unwrap(model),p,T,w)/RT + g_ideal,vl
    elseif is_vapour(phase)
        if isnan(vol)
            volw = volume(model,p,T,w,phase = phase)
        else
            volw = vol
        end
        ∑zlogϕi,vv = ∑zlogϕ(gas_model(model),p,T,w,phase = :v,vol = volw)
        return ∑zlogϕi + tpd_delta_g_vapour(wrapper,p,T,w) + g_ideal,vv
    elseif is_unknown(phase)
        ∑zlogϕi,vv = ∑zlogϕ(gas_model(model),p,T,w,phase = :v)
        gl = excess_gibbs_free_energy(__γ_unwrap(model),p,T,w)/RT + g_ideal
        gv = ∑zlogϕi + tpd_delta_g_vapour(wrapper,p,T,w) + g_ideal
        if gl < gv
            return gl,vl
        else
            return gv,vv
        end
    else
        throw(error("invalid phase specification: $phase"))
    end
end

function identify_phase(wrapper::PTFlashWrapper, p::Number, T, w=SA[1.]; vol0=nothing, vol = NaN)
    model = wrapper.model
    TT = Base.promote_eltype(wrapper,p,T,w)
    RT = Rgas(model)*T
    ∑w = sum(w)
    #g_ideal = sum(xlogx,w) - xlogx(∑w)
    vl = zero(TT)
    if isnan(vol)
        vv = volume(gas_model(model),p,T,w,phase = :v,vol0 = vol0)
    else
        vv = TT(vol)
    end
    ∑zlogϕi,_ = ∑zlogϕ(gas_model(model),p,T,w,phase = :v,vol = vv)
    gl = excess_gibbs_free_energy(__γ_unwrap(model),p,T,w)/RT #+ g_ideal
    gv = ∑zlogϕi + tpd_delta_g_vapour(wrapper,p,T,w) #+g_ideal
    if gl < gv
        return :liquid
    else
        return :vapour
    end
end

function saturation_pressure_ad2(result,model,T)
    return saturation_pressure_ad(result,(model,T),(model,primalval(T)))
end

function tpd_delta_d_vapour!(d,wrapper,p,T)
    lnϕsat,sat = wrapper.fug,wrapper.sat
    is_ideal = gas_model(wrapper) isa IdealModel
    RT = Rgas(gas_model(wrapper))*T
    for i in eachindex(d)
        ps,vl,vv = sat[i]
        Δd = log(ps/p)
        is_ideal || (Δd += vl*(p - ps)/RT + lnϕsat[i])
        d[i] = d[i] - Δd
    end
    return d
end

function tpd_∂delta_d∂P_vapour!(d,wrapper,p,T)
    sat = wrapper.sat
    pure = wrapper.pures
    is_ideal = gas_model(wrapper) isa IdealModel
    RT = Rgas(gas_model(wrapper))*T
    for i in eachindex(d)
        ps,vl,vv = sat[i]
        Δd = -1/p
        is_ideal || (Δd += vl/RT)
        d[i] = d[i] - Δd
    end
    return d
end

function tpd_∂delta_d∂T_vapouri(model,sat,p,T)
    is_ideal = gas_model(model) isa IdealModel
    function f(_T)
        RT = Rgas(model)*_T
        ps,vl,vv = saturation_pressure_ad2(sat,model,_T)
        gasmodel = gas_model(model)
        Δd = log(ps/p)
        if gasmodel isa IdealModel 
            Δd += vl*(p - ps)/RT + VT_lnϕ_pure(gas_model(model),vv,_T,ps)
        end
        return Δd
    end
    return Solvers.derivative(f,T)
end


function tpd_∂delta_d∂T_vapour!(d,wrapper,p,T)
    sat = wrapper.sat
    pure = wrapper.pures
    for i in eachindex(d)
        dΔddT = tpd_∂delta_d∂T_vapouri(pure[i],sat[i],p,T)
        d[i] = d[i] - dΔddT
    end
    return d
end

function tpd_delta_g_vapour(wrapper::PTFlashWrapper,p,T,w)
    lnϕsat,sat = wrapper.fug,wrapper.sat
    pure = wrapper.pures
    gasmodel = gas_model(wrapper.model)

    is_ideal = gasmodel isa IdealModel
    RT = Rgas(gasmodel)*T
    res = zero(Base.promote_eltype(gasmodel,p,T,w))
    for i in eachindex(w)
        pure_i = pure[i]
        ps,vl,vv = saturation_pressure_ad2(sat[i],pure_i,T)
        lnϕsat_i = if T isa ForwardDiff.Dual
            VT_lnϕ_pure(pure_i,vv,T,ps)
        else
            lnϕsat[i]*one(res)
        end
        Δd = lnϕsat_i + log(ps/p)
        is_ideal || (Δd += vl*(p - ps)/RT)
        res -= w[i]*Δd
    end
    return res
end

function ∂lnϕ∂n∂P∂T(wrapper::PTFlashWrapper, p, T, z=SA[1.],cache = ∂lnϕ_cache(wrapper,p,T,z,Val{true}());
            phase=:unknown,
            vol0=nothing,
            threaded = true,
            vol = nothing)

    if is_liquid(phase)
        result,aux,logγ,A1,x1,x2,∂lnγ∂P,hconfig = cache
        g_E,lnγ,∂lnγ∂ni,∂lnγ∂T = ∂lnγ∂n∂T(__γ_unwrap(wrapper), p, T, z,cache)
        ∂lnγ∂P .= 0
        V = zero(typeof(g_E))
        return lnγ,∂lnγ∂ni,∂lnγ∂P,∂lnγ∂T,V
    else
        if vol === nothing
            _vol = volume(gas_model(wrapper),p,T,z;phase,vol0,threaded)
        else
            _vol = vol
        end
        lnϕ, ∂lnϕ∂n, ∂lnϕ∂P, ∂lnϕ∂T, V = ∂lnϕ∂n∂P∂T(gas_model(wrapper), p, T, z,cache; vol = _vol)
        tpd_delta_d_vapour!(lnϕ,wrapper,p,T)
        tpd_∂delta_d∂P_vapour!(∂lnϕ∂P,wrapper,p,T)
        tpd_∂delta_d∂T_vapour!(∂lnϕ∂T,wrapper,p,T)
        return lnϕ, ∂lnϕ∂n, ∂lnϕ∂P, ∂lnϕ∂T, V
    end
end

function ∂lnϕ∂n∂P(wrapper::PTFlashWrapper, p, T, z=SA[1.],cache = ∂lnϕ_cache(wrapper,p,T,z,Val{false}());
            phase=:unknown,
            vol0=nothing,
            threaded = true,
            vol = nothing)


    if is_liquid(phase)
        result,aux,logγ,A1,x1,x2,∂lnγ∂P,hconfig = cache
        g_E,lnγ,∂lnγ∂ni = ∂lnγ∂n(__γ_unwrap(wrapper), p, T, z,cache)
        ∂lnγ∂P .= 0
        V = zero(typeof(g_E))
        return lnγ,∂lnγ∂ni,∂lnγ∂P,V
    else
        if vol === nothing
            _vol = volume(gas_model(wrapper),p,T,z;phase,vol0,threaded)
        else
            _vol = vol
        end
        lnϕ, ∂lnϕ∂n, ∂lnϕ∂P, V = ∂lnϕ∂n∂P(gas_model(wrapper), p, T, z,cache;vol = _vol)
        tpd_delta_d_vapour!(lnϕ,wrapper,p,T)
        tpd_∂delta_d∂P_vapour!(∂lnϕ∂P,wrapper,p,T)
        return lnϕ, ∂lnϕ∂n, ∂lnϕ∂P, V
    end
end

function ∂lnϕ∂P(wrapper::PTFlashWrapper, p, T, z=SA[1.], cache = ∂lnϕ_cache(wrapper,p,T,z,Val{false}());
            phase=:unknown,
            vol0=nothing,
            threaded = true,
            vol = volume(wrapper,p,T,z;phase,vol0,threaded))

    if is_liquid(phase)
        result,aux,logγ,A1,x1,x2,∂lnγ∂Pi,hconfig = cache
        ∂lnγ∂Pi .= 0
        V = zero(eltype(∂lnγ∂Pi))
        return ∂lnγ∂Pi,V
    else
        if vol === nothing
            _vol = volume(gas_model(wrapper),p,T,z;phase,vol0,threaded)
        else
            _vol = vol
        end
        ∂lnϕ∂Pi, V = ∂lnϕ∂P(gas_model(wrapper), p, T, z,cache;vol = _vol)
        tpd_∂delta_d∂P_vapour!(∂lnϕ∂Pi,wrapper,p,T)
        return ∂lnϕ∂Pi, V
    end
end

function ∂lnϕ∂T(wrapper::PTFlashWrapper, p, T, z=SA[1.], cache = ∂lnϕ_cache(wrapper,p,T,z,Val{true}());
            phase=:unknown,
            vol0=nothing,
            threaded = true,
            vol = volume(wrapper,p,T,z;phase,vol0,threaded))

    if is_liquid(phase)
        ∂lnϕ∂Ti = ∂lnγ∂T(__γ_unwrap(wrapper),p,T,z,cache)
        V = zero(eltype(∂lnϕ∂Ti))
        return ∂lnϕ∂Ti,V
    else
        if vol === nothing
            _vol = volume(gas_model(wrapper),p,T,z;phase,vol0,threaded)
        else
            _vol = vol
        end
        ∂lnϕ∂Ti, V = ∂lnϕ∂T(gas_model(wrapper), p, T, z, cache;vol = _vol)
        tpd_∂delta_d∂T_vapour!(∂lnϕ∂Ti,wrapper,p,T)
        return ∂lnϕ∂Ti, V
    end
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

#=
10.04327,1616.76,219.54,335.17,394.54
=#