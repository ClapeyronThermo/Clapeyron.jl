function volume_impl(model::PTFlashWrapper, p, T, z, phase, threaded, vol0)
    if is_unknown(phase) || phase == :stable
        new_phase = identify_phase(model, p, T, z)
        return volume_impl(model.model, p, T, z, new_phase, threaded, vol0)
    else
        volume_impl(model.model, p, T, z, phase, threaded, vol0)
    end
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

function eos_g(wrapper::PTFlashWrapper,p,T,z)
    RT = Rgas(wrapper)*T
    g_E = excess_gibbs_free_energy(__γ_unwrap(wrapper),p,T,z) #excess gibbs
    dg = -tpd_delta_g_vapour(wrapper,p,T,z)*RT #difference between gas fugacity and liquid activity, summed over z
    g_ideal = eos_g(idealmodel(wrapper),p,T,z) #ideal gibbs energy
    g0 = reference_state_eval(wrapper,p,T,z) #reference gibbs energy (equal to reference helmholtz energy)
    return g0 + g_ideal + g_E + dg
end

lb_volume(model::PTFlashWrapper,T,z) = lb_volume(fluid_model(model),T,z)

function PT_property(model::PTFlashWrapper,p,T,z,phase,threaded,vol0,f::F,vol::V) where {F,V}
    if length(model) == 1
        fluidmodel = fluid_model(model)
        return PT_property(fluidmodel,p,T,z,phase,threaded,vol0,f,vol) + Δref(model,fluidmodel,p,T,z,f)
    end

    if phase == :stable || is_unknown(phase)
        new_phase = identify_phase(model,p,T,z,vol0 = vol0)
        return PT_property(model,p,T,z,new_phase,threaded,vol0,f,vol)
    end

    #shortcut for one-component models:

    #=
    Vapour properties are calculated with the fluid model
    Liquid properties are calculated via eos_g(PTFlashWrapper,p,T,z)
    =#
    if is_vapour(phase)
        gasmodel = gas_model(model)
        return PT_property(gasmodel,p,T,z,phase,threaded,vol0,f,vol) + Δref(model,gasmodel,p,T,z,f)

    elseif is_liquid(phase)
        if __γ_unwrap(model) isa IdealLiquidSolution
            #liquid phase + no activity: just delegate to the liquid model, whatever that model may be
            #even for saturated liquid volumes, you can get some props
            return PT_property(liquid_model(model),p,T,z,phase,threaded,vol0,f,vol)
        end

        return PT_property_gibbs(model,p,T,z,f)
        #return PT_property(ActivityModelAresWrapper(model),p,T,z,phase,threaded,vol0,f)
    else
        throw(error("invalid phase specifier: $phase"))
    end
end
