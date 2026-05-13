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
