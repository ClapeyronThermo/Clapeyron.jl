"""
    GammaPhi{γ,Φ} <: RestrictedEquilibriaModel

Wrapper struct to signal that a `CompositeModel` uses an activity model in conjunction with a fluid.
"""
struct GammaPhi{γ,Φ} <: RestrictedEquilibriaModel
    components::Vector{String}
    activity::γ
    fluid::EoSVectorParam{Φ}
end

reference_state(model::GammaPhi) = reference_state(model.fluid)

function Base.show(io::IO,mime::MIME"text/plain",model::GammaPhi)
    print(io,"γ-ϕ Model")
    length(model) == 1 && print(io, " with 1 component:")
    length(model) > 1 && print(io, " with ", length(model), " components:")
    println(io)
    show_pairs(io,model.components)
    act = model.activity
    if hasfield(typeof(act),:puremodel)
        print(io,'\n',"Activity Model: ", parameterless_type(act))
    else
        print(io,'\n',"Activity Model: ",typeof(act))
    end
    print(io,'\n',"Fluid Model: ",typeof(model.fluid.model))
    show_reference_state(io,model;space = true)
end

fluid_model(model::GammaPhi) = model.fluid.model
__γ_unwrap(model::GammaPhi) = __γ_unwrap(model.activity)
gas_model(model::GammaPhi) = gas_model(model.fluid.model)
Base.eltype(model::GammaPhi) = Base.promote_eltype(__γ_unwrap(model),gas_model(model))

function excess_gibbs_free_energy(model::GammaPhi,p,T,z)
    return excess_gibbs_free_energy(model.activity,p,T,z)
end

reference_chemical_potential_type(model::GammaPhi) = reference_chemical_potential_type(model.activity)

function volume_impl(model::GammaPhi,p,T,z,phase,threaded,vol0)
    return volume_impl(model.fluid.model,p,T,z,phase,threaded,vol0)
end

molecular_weight(model::GammaPhi,z) = molecular_weight(model.fluid.model,z)
saturation_model(model::GammaPhi) = saturation_model(model.fluid)
idealmodel(model::GammaPhi) = idealmodel(model.fluid.model)

function init_preferred_method(method::typeof(tp_flash),model::GammaPhi,kwargs)
    MichelsenTPFlash(;kwargs...)
end

# Error handling for Activity models that don't provide saturation properties, in the context of VLE.
function ActivitySaturationError(model,method)
    throw(ArgumentError("$method requires $model to be used along with another EoS model that supports saturation properties. If you are using an Activity Model as a raw input, use `CompositeModel(components, liquid = activity_model, fluid = fluid_model)` instead."))
end

function gibbs_solvation(model::GammaPhi,T)
    z = [1.0,1e-30]
    p,v_l,v_v = saturation_pressure(model.fluid[1],T)
    p2,v_l2,v_v2 = saturation_pressure(model.fluid[2],T)
    γ = activity_coefficient(model,p,T,z)
    K = v_v/v_l*γ[2]*p2/p
    return -R̄*T*log(K)
end

#=
some equations require more delicate handling. the general case only works for extensive properties.
=#
function PT_property_gammaphi(model::GammaPhi,p,T,z,f::F,USEP) where F
    ∑z = sum(z)
    z1 = SA[∑z]
    res = zero(Base.promote_eltype(model,p,T,z))
    Vl = zero(res)
    use_p = USEP === Val{true}()
    for i in 1:length(model)
        zi = z[i]
        modeli = model.fluid.pure[i]
        #we suppose this is liquid phase
        Vi = volume(modeli, p, T, z1; phase = :liquid,threaded = false)
        Vl += Vi*zi/∑z
        if use_p
            res += f(modeli,Vi,T,z1,p)*zi/∑z
        else
            res += f(modeli,Vi,T,z1)*zi/∑z
        end
    end
    if !iszero(p)
        if use_p
            res += f(model.activity,0.0,T,z,p)
        else
            res += f(model.activity,0.0,T,z)
        end
        return res
    end
end

function gammaphi_f_hess(model::GammaPhi,p,T,z)
    ∑z = sum(z)
    z1 = SA[∑z]
    _0 = zero(Base.promote_eltype(model,p,T,z))
    _1 = oneunit(_0)
    Vl = zero(_0)
    ∂²A∂T² = ∂²f∂T²(model.activity,_0,T,z)
    ∂²A = SMatrix{2,2}((_0,_0,_0,_0 + ∂²A∂T²))
    for i in 1:length(model)
        zi = z[i]
        modeli = model.fluid.pure[i]
        #eos(model) = sum(xi*eos(model,Vi,T,1.0))

        Vi = volume(modeli, p, T, z1; phase = :liquid, threaded = false)
        Vl += Vi*zi/∑z
        Vii = Vi*zi
        ∂²Aᵢ = f_hess(modeli,Vi,T,z1) .* zi ./ ∑z
        ∂²A = ∂²A .+ ∂²Aᵢ
    end
    return ∂²A,Vl
end

function PT_property_gammaphi(model,p,T,z,::typeof(VT_isobaric_heat_capacity),USEP)
    d²A,V = gammaphi_f_hess(model,p,T,z)
    ∂²A∂V∂T = d²A[1,2]
    ∂²A∂V² = d²A[1,1]
    ∂²A∂T² = d²A[2,2]
    return -T*(∂²A∂T² - ∂²A∂V∂T^2/∂²A∂V²)
end

function PT_property_gammaphi(model,p,T,z,::typeof(VT_adiabatic_index),USEP)
    d²A,V = gammaphi_f_hess(model,p,T,z)
    ∂²A∂V∂T = d²A[1,2]
    ∂²A∂V² = d²A[1,1]
    ∂²A∂T² = d²A[2,2]
    return 1 - ∂²A∂V∂T*∂²A∂V∂T/(∂²A∂V²*∂²A∂T²)
end

function PT_property_gammaphi(model,p,T,z,::typeof(VT_isothermal_compressibility),USEP)
    d²A,V = gammaphi_f_hess(model,p,T,z)
    return -1/V/d²A[1,1]
end

function PT_property_gammaphi(model,p,T,z,::typeof(VT_isentropic_compressibility),USEP)
    d²A,V = gammaphi_f_hess(model,p,T,z)
    ∂²A∂V∂T = d²A[1,2]
    ∂²A∂V² = d²A[1,1]
    ∂²A∂T² = d²A[2,2]
    return 1/V/(∂²A∂V²-∂²A∂V∂T^2/∂²A∂T²)
end

function PT_property_gammaphi(model,p,T,z,::typeof(VT_speed_of_sound),USEP)
    Mr = molecular_weight(model,z)
    d²A,V = gammaphi_f_hess(model,p,T,z)
    ∂²A∂V∂T = d²A[1,2]
    ∂²A∂V² = d²A[1,1]
    ∂²A∂T² = d²A[2,2]
    return V*sqrt((∂²A∂V²-∂²A∂V∂T^2/∂²A∂T²)/Mr)
end

function PT_property_gammaphi(model,p,T,z,::typeof(VT_isobaric_expansivity),USEP)
    d²A,V = gammaphi_f_hess(model,p,T,z)
    ∂²A∂V∂T = d²A[1,2]
    ∂²A∂V² = d²A[1,1]
    return -∂²A∂V∂T/(V*∂²A∂V²)
end

function PT_property_gammaphi(model,p,T,z,::typeof(VT_joule_thomson_coefficient),USEP)
    d²A,V = gammaphi_f_hess(model,p,T,z)
    ∂²A∂V∂T = d²A[1,2]
    ∂²A∂V² = d²A[1,1]
    ∂²A∂T² = d²A[2,2]
    return -(∂²A∂V∂T - ∂²A∂V²*((T*∂²A∂T² + V*∂²A∂V∂T) / (T*∂²A∂V∂T + V*∂²A∂V²)))^-1
end

function PT_property(model::GammaPhi,p,T,z,phase,threaded,vol0,f::F,USEP::Val{UseP}) where {F,UseP}
    if phase == :stable
        #we dont support this in GammaPhi.
        throw(error("automatic phase detection not implemented for $(typeof(model))"))
    end

    #shortcut for one-component models:
    if length(model) == 1
        res_v = PT_property(model.fluid,p,T,z,phase,threaded,vol0,f,USEP)
    end
    #=
    Calculate PT properties as the following:
    prop = ∑xi*prop(pure[i],p,T) + excess_prop(activity,T,x) + prop(reference_state,T)

    Vapour properties are calculated with the fluid model
    =#
    if is_vapour(phase)
        res = PT_property(model.fluid,p,T,z,phase,threaded,vol0,f,USEP)
        return res
    elseif is_liquid(phase) || is_unknown(phase)
       return PT_property_gammaphi(model,p,T,z,f,USEP)
    else
        throw(error("invalid phase specifier: $phase"))
    end
end

function __calculate_reference_state_consts(model::GammaPhi,v,T,p,z,H0,S0,phase)
    ∑z = sum(z)
    S00 = entropy(model,p,T,z,phase = phase)
    a1 = (S00 - S0)#/∑z
    H00 = enthalpy(model,p,T,z,phase = phase)
    a0 = (-H00 + H0)#/∑z
    return a0,a1
end

function PTFlashWrapper(model::GammaPhi,p,T,z,equilibrium)
    fluidmodel = model.fluid
    #check that we can actually solve the equilibria
    if fluidmodel isa IdealModel && !is_lle(equilibrium)
        ActivitySaturationError(model.activity,tp_flash)
    end
    TT = Base.promote_eltype(model,p,T,z)
    wrapper = PTFlashWrapper{TT}(model,equilibrium,fluidmodel.pure)
    if is_vle(equilibrium) || is_unknown(equilibrium)
        update_temperature!(wrapper,T)
    end
    return wrapper
end

function __tpflash_cache_model(model::GammaPhi,p,T,z,equilibrium)
    PTFlashWrapper(model,p,T,z,equilibrium)
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

function __tpflash_gibbs_reduced(wrapper::PTFlashWrapper{<:GammaPhi},p,T,x,y,β,eq,vols)
    model = wrapper.model
    gibbs = zero(Base.promote_eltype(model,p,T,x,β))
    if !isone(β)
        gx,_ = modified_gibbs(wrapper,p,T,x,:l)
        gibbs += gx*(1-β)
    end

    if is_vle(eq) && !iszero(β)
        vv = vols[2]
        gy,_ = modified_gibbs(wrapper,p,T,y,:v,vols[2])
        gibbs += gy*β
    elseif !iszero(β) #lle
        gy,_ = modified_gibbs(wrapper,p,T,y,:l)
        gibbs += gy*β
    end
    return gibbs
end

function K0_lle_init(wrapper::PTFlashWrapper,p,T,z)
    return K0_lle_init(__γ_unwrap(wrapper),p,T,z)
end

function __eval_G_DETPFlash(wrapper::PTFlashWrapper,p,T,x,equilibrium)
    model = wrapper.model
    phase = is_lle(equilibrium) ? :liquid : :unknown
    return modified_gibbs(wrapper,p,T,x,phase)
end

function modified_gibbs(wrapper::PTFlashWrapper,p,T,w,phase = :unknown,vol = NaN)
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

function identify_phase(wrapper::PTFlashWrapper, p::Number, T, w=SA[1.]; phase=:unknown, threaded=true, vol0=nothing, vol = NaN)
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

function tpd_delta_d_vapour!(d,wrapper,p,T)
    ϕsat,sat = wrapper.fug,wrapper.sat
    is_ideal = gas_model(wrapper) isa IdealModel
    RT = Rgas(gas_model(wrapper))*T
    for i in eachindex(d)
        ps,vl,vv = sat[i]
        Δd = log(ps/p)
        is_ideal || (Δd += vl*(p - ps)/RT + log(ϕsat[i]))
        d[i] = d[i] - Δd
    end
    return d
end

function tpd_∂delta_d∂P_vapour!(d,wrapper,p,T)
    ϕsat,sat = wrapper.fug,wrapper.sat
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

function tpd_∂delta_d∂T_vapouri(model,p,T)
    function f(_T)
        RT = Rgas(model)*_T
        ps,vl,vv = saturation_pressure(model,_T)
        Δd = log(ps/p) + VT_lnϕ_pure(gas_model(model),vv,_T,ps)
        if gas_model(model) isa IdealModel
            Δd += vl*(p - ps)/RT
        end
        return Δd
    end
    return Solvers.derivative(f,T)
end


function tpd_∂delta_d∂T_vapour!(d,wrapper,p,T)
    ϕsat,sat = wrapper.fug,wrapper.sat
    pure = wrapper.pures
    for i in eachindex(d)
        dΔddT = tpd_∂delta_d∂T_vapouri(pure[i],p,T)
        d[i] = d[i] - dΔddT
    end
    return d
end

function tpd_delta_g_vapour(wrapper,p,T,w)
    ϕsat,sat = wrapper.fug,wrapper.sat
    is_ideal = gas_model(wrapper) isa IdealModel
    RT = Rgas(gas_model(wrapper))*T
    res = zero(Base.promote_eltype(gas_model(wrapper),p,T,w))
    for i in eachindex(w)
        ps,vl,vv = sat[i]
        Δd = log(ϕsat[i]) + log(ps/p)
        is_ideal || (Δd += vl*(p - ps)/RT)
        res -= w[i]*Δd
    end
    return res
end



function tpd_input_composition(wrapper::PTFlashWrapper{<:GammaPhi},p,T,z,lle,cache = tpd_cache(wrapper,p,T,z,di))

    TT = Base.promote_eltype(wrapper.model,p,T,z)

    pures = wrapper.model.fluid.pure
    model = wrapper.model
    fluidmodel = model.fluid.model
    RT = R̄*T

    d_l,d_v,_,_,_,Hϕ = cache

    TT = Base.promote_eltype(model,p,T,z)

    n = sum(z)
    logsumz = log(n)
    d,vl = tpd_lnϕ_and_v!(last(cache),wrapper,p,T,z,nothing,false,:liquid)
    d_l .= d
    d_l .+= log.(z) .- logsumz

    lle && return copy(d_l),:liquid,vl

    d,vv = tpd_lnϕ_and_v!(last(cache),wrapper,p,T,z,nothing,false,:vapour)
    d_v .= d
    d_v .+= log.(z) .- logsumz
    gr_l = dot(z,d_l)
    gr_v = dot(z,d_v)
    if gr_l < gr_v
        return copy(d_l),:liquid,vl
    else
        return copy(d_v),:vapour,vv
    end
end

function tpd_lnϕ_and_v!(cache,wrapper::PTFlashWrapper,p,T,w,vol0,liquid_overpressure = false,phase = :l,_vol = nothing)
    model = wrapper.model
    RT = R̄*T
    if is_liquid(phase)
        γmodel = __γ_unwrap(model)
        #=
        If the model is not an activity model, then PTFlashWrapper is wrapping
        a normal helmholtz model, we just return lnϕ.
        =#
        if γmodel isa ActivityModel
            logγx = lnγ(γmodel,p,T,w,cache)
            v = zero(eltype(logγx))
            return logγx,v,true
        elseif is_vle(wrapper.equilibrium)
            logγx,v = __lnγ_sat(wrapper,p,T,w,cache)
            return logγx,v,true
        end
    end
    fxy,v,overpressure = tpd_lnϕ_and_v!(cache,gas_model(model),p,T,w,vol0,liquid_overpressure,phase,_vol)
    is_vapour(phase) && !is_lle(wrapper.equilibrium) && tpd_delta_d_vapour!(fxy,wrapper,p,T)

    return fxy,v,overpressure
end

function __lnγ_sat(wrapper::PTFlashWrapper,p,T,w,cache = nothing,vol0 = nothing,vol = volume(wrapper.model,p,T,w,vol0 = vol0,phase = :l))
    model = wrapper.model
    μmix_temp = VT_chemical_potential_res!(cache,model,vol,T,w)
    result,aux,logγ,A1,μmix,x2,x3,hconfig = cache
    μmix .= μmix_temp
    sat = wrapper.sat
    fug = wrapper.fug
    RT = Rgas(model)*T
    for i in 1:length(logγ)
        ϕᵢ = fug[i]
        pᵢ,vpureᵢ,_ = sat[i]

        μᵢ_over_RT = log(ϕᵢ) + log(pᵢ*vpureᵢ/RT)
        logγ[i] = log(vpureᵢ/vol) + μmix[i]/RT - μᵢ_over_RT -  vpureᵢ*(p - pᵢ)/RT
    end
    return logγ,vol
end

function modified_∂lnϕ∂n(wrapper::PTFlashWrapper{<:GammaPhi}, p, T, z, cache; phase = :unknown, vol0 = nothing)
    model = wrapper.model
    if is_vapour(phase)
        lnϕ,∂lnϕ∂n,vol =  modified_∂lnϕ∂n(gas_model(model),p,T,z,cache;phase,vol0)
        tpd_delta_d_vapour!(lnϕ,wrapper,p,T)
        return lnϕ,∂lnϕ∂n,vol
    elseif is_liquid(phase)
        g_E,lnγ,∂lnγ∂ni = ∂lnγ∂n(__γ_unwrap(model),p,T,z,cache)
        return lnγ,∂lnγ∂ni,zero(g_E)
    else
        throw(error("invalid specification for phase: $phase"))
    end
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
    fug = wrapper.fug
    for i in 1:nc
        psi,vli,_ = sat[i]
        zi = z[i]
        ∑zlogϕsat += zi*log(fug[i])
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

export GammaPhi
