"""
    GammaPhi{γ,Φ} <: RestrictedEquilibriaModel

wrapper struct to signal that a `CompositeModel` uses an activity model in conjunction with a fluid.
"""
struct GammaPhi{γ,Φ} <: RestrictedEquilibriaModel
    components::Vector{String}
    activity::γ
    fluid::EoSVectorParam{Φ}
end

__gas_model(model::GammaPhi) = model.fluid

function init_preferred_method(method::typeof(saturation_pressure),model::GammaPhi,kwargs)
    return init_preferred_method(saturation_pressure,model.fluid,kwargs)
end

function init_preferred_method(method::typeof(saturation_temperature),model::GammaPhi,kwargs)
    return init_preferred_method(saturation_temperature,model.fluid,kwargs)
end

function activity_coefficient(model::GammaPhi,p,T,z=SA[1.]; phase = :unknown, threaded=true)
    return activity_coefficient(model.activity,p,T,z)
end

function saturation_pressure(model::GammaPhi,T,method::SaturationMethod)
    return saturation_pressure(model.fluid,T,method)
end

crit_pure(model::GammaPhi) = crit_pure(model.fluid)
x0_sat_pure(model::GammaPhi,T) = x0_sat_pure(model.fluid)
x0_psat(model::GammaPhi,T) = x0_psat(model.fluid,T)

function saturation_temperature(model::GammaPhi,p,method::SaturationMethod)
    return saturation_temperature(model.fluid,p,method)
end

function init_preferred_method(method::typeof(tp_flash),model::GammaPhi,kwargs)
    RRTPFlash(;kwargs...)
end

# Error handling for Activity models that don't provide saturation properties, in the context of VLE.
function ActivitySaturationError(model,method)
    throw(ArgumentError("$method requires $model to be used in conjuction with another EoS model that supports saturation properties. If you are using an Activity Model as a raw input, use `CompositeModel(components, liquid = activity_model, fluid = fluid_model)` instead."))
end

function gibbs_solvation(model::GammaPhi,T)
    z = [1.0,1e-30]
    p,v_l,v_v = saturation_pressure(model.fluid[1],T)
    p2,v_l2,v_v2 = saturation_pressure(model.fluid[2],T)
    γ = activity_coefficient(model,p,T,z)
    K = v_v/v_l*γ[2]*p2/p
    return -R̄*T*log(K)
end

function PTFlashWrapper(model::GammaPhi,p,T::Number,equilibrium::Symbol)
    fluidmodel = model.fluid
    #check that we can actually solve the equilibria
    if fluidmodel isa IdealModel && !is_lle(equilibrium)
        ActivitySaturationError(model.activity,tp_flash)
    end
    pures = fluidmodel.pure
    RT = R̄*T
    if fluidmodel.model isa IdealModel
        vv = RT/p
        nan = zero(vv)/zero(vv)
        sats = fill((nan,nan,vv),length(model))
        ϕpure = fill(one(vv),length(model))
        g_pure = [VT_gibbs_free_energy(__gas_model(pures[i]),vv,T) for i in 1:length(model)]
        return PTFlashWrapper(model.components,model,sats,ϕpure,g_pure)
    else
        sats = saturation_pressure.(pures,T)
        vv_pure = last.(sats)
        p_pure = first.(sats)
        μpure = only.(VT_chemical_potential_res.(__gas_model.(pures),vv_pure,T))
        ϕpure = exp.(μpure ./ RT .- log.(p_pure .* vv_pure ./ RT))
        g_pure = [VT_gibbs_free_energy(__gas_model(pures[i]),vv_pure[i],T) for i in 1:length(model)]
        return PTFlashWrapper(model.components,model,sats,ϕpure,g_pure)
    end

end

__tpflash_cache_model(model::GammaPhi,p,T,z,equilibrium) = PTFlashWrapper(model,p,T,equilibrium)

function update_K!(lnK,wrapper::PTFlashWrapper{<:GammaPhi},p,T,x,y,volx,voly,phasex,phasey,β = nothing,inx = FillArrays.Fill(true,length(x)),iny = inx)
    model = wrapper.model
    fluidmodel = model.fluid.model
    sats = wrapper.sat
    g_pures = wrapper.μ
    n = length(model)
    #crits = wrapper.crit
    fug = wrapper.fug
    RT = R̄*T
    γx = activity_coefficient(model, p, T, x)
    volx = volume(fluidmodel, p, T, x, phase = phasex, vol0 = volx)
    _0 = zero(eltype(lnK))

    if β === nothing
        _0 = zero(eltype(lnK))
        gibbs = _0/_0
    else
        gibbs = _0
        for i in eachindex(x)
            if inx[i]
                g_E_x = x[i]*RT*log(γx[i])
                g_ideal_x = x[i]*RT*log(x[i])
                g_pure_x = x[i]*g_pures[i]
                gibbs += (g_E_x + g_ideal_x + g_pure_x)*(1-β)/RT
            end
        end
    end

    if is_vapour(phasey)
        lnϕy, voly = lnϕ(__gas_model(fluidmodel), p, T, y; phase=phasey, vol0=voly)
        for i in eachindex(lnK)
            if iny[i]
                ϕli = fug[i]
                p_i = sats[i][1]
                lnK[i] = log(γx[i]*p_i*ϕli/p) - lnϕy[i] + volx*(p - p_i)/RT
                gibbs += β*y[i]*(log(y[i]) + lnϕy[i])
            end
        end
    else
        γy = activity_coefficient(model, p, T, y)
        lnK .= log.(γx./γy)
        voly = volume(fluidmodel, p, T, y, phase = phasey, vol0 = voly)
        if β !== nothing
            for i in eachindex(y)
                if iny[i]
                    g_E_y = y[i]*RT*log(γy[i])
                    g_ideal_y = y[i]*RT*(log(y[i]))
                    g_pure_y = y[i]*g_pures[i]
                    gibbs += (g_E_y + g_ideal_y + g_pure_y)*β/RT
                end
            end
        end
    end

    return lnK,volx,voly,gibbs
end

function __tpflash_gibbs_reduced(wrapper::PTFlashWrapper{<:GammaPhi},p,T,x,y,β,eq)
    pures = wrapper.model.fluid.pure
    model = wrapper.model
    fluidmodel = model.fluid.model
    g_pures = wrapper.μ

    γx = activity_coefficient(model.activity, p, T, x)
    RT = R̄*T
    n = length(model)
    g_E_x = sum(x[i]*RT*log(γx[i]) for i ∈ 1:n)
    g_ideal_x = sum(x[i]*RT*(log(x[i])) for i ∈ 1:n)
    g_pure_x = dot(x,g_pures)
    gibbs = (g_E_x + g_ideal_x + g_pure_x)*(1-β)/RT
    if is_vle(eq)
        gibbs += gibbs_free_energy(__gas_model(fluidmodel),p,T,y,phase =:v)*β/R̄/T
    else #lle
        γy = activity_coefficient(model.activity, p, T, y)
        g_E_y = sum(y[i]*RT*log(γy[i]) for i ∈ 1:n)
        g_ideal_y = sum(y[i]*R̄*T*(log(y[i])) for i ∈ 1:n)
        g_pure_y = dot(y,g_pures)
        gibbs += (g_E_y + g_ideal_y + g_pure_y)*β/RT
    end
    return gibbs
    #(gibbs_free_energy(model,p,T,x)*(1-β)+gibbs_free_energy(model,p,T,y)*β)/R̄/T
end

#TODO: derive expressions for this

function dgibbs_obj!(model::PTFlashWrapper{<:GammaPhi}, p, T, z, phasex, phasey,
    nx, ny, vcache, ny_var = nothing, in_equilibria = FillArrays.Fill(true,length(z)), non_inx = in_equilibria, non_iny = in_equilibria;
    F=nothing, G=nothing, H=nothing)
    throw(error("γ-ϕ Composite Model don't support gibbs energy optimization in MichelsenTPFlash."))
end

function K0_lle_init(wrapper::PTFlashWrapper{<:GammaPhi},p,T,z)
    return K0_lle_init(wrapper.model.activity,p,T,z)
end

function __eval_G_DETPFlash(wrapper::PTFlashWrapper{<:GammaPhi},p,T,x,equilibrium)
    model = wrapper.model
    phase = is_lle(equilibrium) ? :liquid : :unknown
    n = length(model)
    g_pures = wrapper.μ
    R = Rgas()
    RT = R*T
    γx = activity_coefficient(model.activity, p, T, x)
    g_E_x = sum(x[i]*RT*log(γx[i]) for i ∈ 1:n)
    g_ideal_x = sum(x[i]*RT*(log(x[i])) for i ∈ 1:n)
    g_pure_x = dot(x,g_pures)
    gl = (g_E_x + g_ideal_x + g_pure_x)
    vl = volume(model.fluid.model,p,T,x,phase = :l)
    if phase == :liquid
        return gl,vl
    else
        throw(error("γ-ϕ Composite Model does not support VLE calculation with `DETPFlash`. if you want to calculate LLE equilibria, try using `DETPFlash(equilibrium = :lle)`"))
        #=  
        vv = volume(model.fluid.model,p,T,x,phase = :v)
        gv = VT_gibbs_free_energy(model.fluid.model, vv, T, x)
        if gv > gl
            return gl,vl
        else
            return gv,vv
        end
        =#
    end
end

export GammaPhi
