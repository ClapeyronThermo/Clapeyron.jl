struct GEDepartureParam <: EoSParam
    vref::SingleParam{Float64}
end

struct GEDeparture{𝔸} <: MultiFluidDepartureModel
    components::Vector{String}
    params::GEDepartureParam
    activity::𝔸
    references::Vector{String}
end

function GEDeparture_constructor(f)
    function departure(components;userlocations = String[],verbose = false)
        return GEDeparture(components;activity = f,userlocations = userlocations,verbose = verbose)
    end
    return departure
end

GEDeparture(f::Function) = GEDeparture_constructor(f)

GEDeparture(f::Type{T}) where T <:ActivityModel = GEDeparture_constructor(T)

#activity is already initialized, construct directly
function GEDeparture(act::ActivityModel)
    init_activity = act
    if has_groups(init_activity)
        components = init_activity.groups.components
    else
        components = init_activity.components
    end
    vref = SingleParam("reference volume",components)
    pkgparams = GEDepartureParam(vref)
    references = ["10.1021/acs.iecr.1c01186","10.1016/j.fluid.2018.04.015"]
    return GEDeparture(components,pkgparams,init_activity,references)
end
"""
GEDeparture <: MultiFluidDepartureModel
    GEDeparture(components;
    activity = UNIFAC,
    userlocations = String[],
    verbose = false)

## Input parameters
none
- `k1`: Pair Parameter (`Float64`) - binary, T-dependent interaction parameter `[K⁻¹]`
## Model parameters
- `vref`: Single Parameter (`Float64`, calculated) - Reference pure molar volume `[m³·mol⁻¹]`
## Input models
- `activity`: activity coefficient model

## Description

Departure that uses the residual excess gibbs energy from an activity model:

```
aᵣ = ∑xᵢaᵣᵢ(δ,τ) + Δa
Δa = gᴱᵣ/RT - log(1+bρ)/log(1+bρref) * ∑xᵢ(aᵣᵢ(δref,τ) - aᵣᵢ(δrefᵢ,τᵢ))
τᵢ = Tcᵢ/T
δref = ρref/ρr
δrefᵢ = ρrefᵢ/ρcᵢ
b = 1/1.17ρref
```

## References
1. Jäger, A., Breitkopf, C., & Richter, M. (2021). The representation of cross second virial coefficients by multifluid mixture models and other equations of state. Industrial & Engineering Chemistry Research, 60(25), 9286–9295. [doi:10.1021/acs.iecr.1c01186](https://doi.org/10.1021/acs.iecr.1c01186)
"""
function GEDeparture(components::AbstractVector; activity = UNIFAC, userlocations = String[], verbose::Bool=false)
    init_activity = init_model(activity,components,userlocations,verbose)
    comps = init_activity.components
    vref = SingleParam("reference volume",comps)
    pkgparams = GEDepartureParam(vref)
    references = ["10.1021/acs.iecr.1c01186","10.1016/j.fluid.2018.04.015"]
    return GEDeparture(comps,pkgparams,init_activity,references)
end

function multiparameter_a_res(model,V,T,z,departure::GEDeparture,δ,τ,∑z = sum(z))
    lnδ = log(δ)
    lnτ = log(τ)
    ∑z⁻¹ = 1/∑z
    aᵣ = multiparameter_a_res0(model,V,T,z,δ,τ,lnδ,lnτ,∑z)
    Vᵣ = v_scale(model,z,model.mixing,∑z)
    _0 = zero(aᵣ)
    isone(length(z)) && return aᵣ
    vref = departure.params.vref.values
    v̄ref = dot(z,vref)*∑z⁻¹
    ρref = 1/v̄ref
    b = 0.8547008547008548*v̄ref #(1/1.17)
    δref = Vᵣ*ρref
    lnδref = log(δref)
    Δa = zero(aᵣ)
    m = model.pures
    Tinv = 1/T
    Tc = model.params.Tr.values
    Vc = model.params.Vr.values
    for i in @comps
        mᵢ = m[i]
        τᵢ = Tc[i]*Tinv
        δrefᵢ = Vc[i]/vref[i]
        Δa += z[i]*(reduced_a_res(mᵢ,δref,τ) - reduced_a_res(mᵢ,δrefᵢ,τᵢ))
    end
    Δa *= ∑z⁻¹
    R = Rgas(model)
    ρ = ∑z/V
    gᴱ = excess_g_res(departure.activity,V,T,z)
    lnb = log1p(b*ρ)/log1p(b*ρref)
    return aᵣ + lnb*(gᴱ*∑z⁻¹/(R*T) - Δa)
end

function lb_volume(model::MultiFluid{A,M,GEDeparture}, T, z) where {A,M}
    vref = model.departure.vref
    v̄ref = dot(z,vref)/sum(z)
    return 0.8547008547008548*v̄ref
end

function recombine_departure!(model::MultiFluid,dep::GEDeparture)
    m = model.pures
    for i in 1:length(m)
        mi = m[i]
        p_atm = 101325.0
        pref = p_atm
        p_triple = mi.properties.ptp
        if p_triple > p_atm
            pref = p_triple + 100.0
        end
        sat = saturation_temperature(mi,pref)
        dep.params.vref[i] = sat[2]
    end
end

export GEDeparture