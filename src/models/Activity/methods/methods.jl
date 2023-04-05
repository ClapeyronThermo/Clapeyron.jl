"""
    __gas_model(model::EoSModel)

internal function.
provides the model used to calculate gas properties.
normally, this is the identity, but `CompositeModel` has a gas model by itself.
"""
__gas_model(model::EoSModel) = model

include("bubble_point.jl")
include("dew_point.jl")
include("LLE_point.jl")


struct ActivityPTFlashWrapper{T,S,C,F} <: EoSModel
    components::Vector{String}
    model::T
    sat::Vector{S}
    crit::Vector{C}
    fug::Vector{F}
end

Base.length(model::ActivityPTFlashWrapper) = length(model.model)


function ActivityPTFlashWrapper(model::ActivityModel,T::Number) 
    pures = model.puremodel.pure
    sats = saturation_pressure.(pures,T)
    crits = crit_pure.(pures)
    vv_pure = last.(sats)
    RT = R̄*T
    p_pure = first.(sats)
    μpure = only.(VT_chemical_potential_res.(__gas_model.(pures),vv_pure,T))
    ϕpure = exp.(μpure ./ RT .- log.(p_pure .* vv_pure ./ RT))
    return ActivityPTFlashWrapper(model.components,model,sats,crits,ϕpure)
end

__tpflash_cache_model(model,p,T,z) = ActivityPTFlashWrapper(model,T)

function update_K!(lnK,wrapper::ActivityPTFlashWrapper,p,T,x,y,volx,voly,phasex,phasey,β = nothing)
    model = wrapper.model
    pures = wrapper.model.puremodel.pure
    sats = wrapper.sat
    n = length(model)
    #crits = wrapper.crit
    fug = wrapper.fug
    RT = R̄*T
    γx = activity_coefficient(model, p, T, x)
    volx = volume(model.puremodel.model, p, T, x, phase = phasex, vol0 = volx)

    if β === nothing
        _0 = zero(eltype(lnK))
        gibbs = _0/_0
    else
        g_E_x = sum(x[i]*RT*log(γx[i]) for i ∈ 1:n)
        g_ideal_x = sum(x[i]*RT*(log(x[i])) for i ∈ 1:n)
        g_pure_x = sum(x[i]*VT_gibbs_free_energy(pures[i],wrapper.sat[i][2],T) for i ∈ 1:n)    
        gibbs = (g_E_x + g_ideal_x + g_pure_x)*(1-β)/RT
    end
    
    if is_vapour(phasey)
        lnϕy, voly = lnϕ(model, p, T, y; phase=phasey, vol0=voly)
        for i in eachindex(lnK)
            ϕli = fug[i]
            p_i = sats[i][1]
            lnK[i] = log(γx*ϕli*exp(volx*(p - p_i)/RT)/exp(lnϕy[i])/p)
            gibbs += β*y[i]*log(y[i] + lnϕy[i])
        end
    else
        γy = activity_coefficient(model, p, T, y)
        lnK .= log.(γx./γy)
        for i in eachindex(lnK)
            lnK[i] = log(γx[i]/γy[i])
        end
        voly = volume(model.puremodel.model, p, T, y, phase = phasey, vol0 = voly)
        if β !== nothing
            g_E_y = sum(y[i]*RT*log(γy[i]) for i ∈ 1:n)
            g_ideal_y = sum(y[i]*R̄*T*(log(y[i])) for i ∈ 1:n)
            g_pure_y = sum(y[i]*VT_gibbs_free_energy(pures[i],wrapper.sat[i][2],T) for i ∈ 1:n)    
            gibbs += (g_E_y + g_ideal_y + g_pure_y)*β/RT
        end
    end
    
    return lnK,volx,voly,gibbs
end

function __tpflash_gibbs_reduced(wrapper::ActivityPTFlashWrapper,p,T,x,y,β,eq)
    pures = wrapper.model.puremodel.pure
    model = wrapper.model
    γx = activity_coefficient(model, p, T, x)
    RT = R̄*T
    n = length(model)
    volx = volume(model.puremodel.model, p, T, x, phase = :liquid)
    g_E_x = sum(x[i]*RT*log(γx[i]) for i ∈ 1:n)
    g_ideal_x = sum(x[i]*RT*(log(x[i])) for i ∈ 1:n)
    g_pure_x = sum(x[i]*VT_gibbs_free_energy(pures[i],wrapper.sat[i][2],T) for i ∈ 1:n)    
    gibbs = (g_E_x + g_ideal_x + g_pure_x)*(1-β)/RT
    if is_vle(eq)
        gibbs += gibbs_free_energy(model.puremodel.model,p,T,y)*β/R̄/T
    else #lle
        voly = volume(model.puremodel.model, p, T, y, phase = :liquid)
        γy = activity_coefficient(model, p, T, y)
        g_E_y = sum(y[i]*RT*log(γy[i]) for i ∈ 1:n)
        g_ideal_y = sum(y[i]*R̄*T*(log(y[i])) for i ∈ 1:n)
        g_pure_y = sum(y[i]*VT_gibbs_free_energy(pures[i],wrapper.sat[i][2],T) for i ∈ 1:n)    
        gibbs += (g_E_y + g_ideal_y + g_pure_y)*β/RT
    end
    return gibbs
    #(gibbs_free_energy(model,p,T,x)*(1-β)+gibbs_free_energy(model,p,T,y)*β)/R̄/T
end