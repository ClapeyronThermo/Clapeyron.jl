struct GEDepartureParam <: EoSParam
    vref::SingleParam{Float64}
end

struct GEDeparture{𝔸} <: MultiFluidDepartureModel
    components::Vector{String}
    params::GEDepartureParam
    activity::𝔸
    references::Vector{String}
end

function GEDeparture(f::F) where F <: Union{Function,ActivityModel}
    function departure(components;userlocations = String[],verbose = false)
        return GEDeparture(components;activity = f,userlocations = userlocations,verbose = verbose)
    end
    return departure
end

function GEDeparture(components; activity = Wilson, userlocations=String[], verbose::Bool=false)
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
    Vᵣ = δ*V*∑z⁻¹
    _0 = zero(aᵣ)
    isone(length(z)) && return aᵣ
    gᴱ = excess_gibbs_free_energy(departure.activity,V,T,z)
    vref = departure.params.vref.values
    v̄ref = dot(z,vref)*∑z⁻¹
    ρref = 1/v̄ref
    b = 0.8547008547008548*v̄ref #(1/1.17)
    δref = Vᵣ*ρref
    lnδref = log(δref)
    Δa = zero(aᵣ)
    m = model.pures
    Tinv = 1/T
    Tc = model.params.Tc.values
    Vc = model.params.Vc.values
    for i in @comps
        mᵢ = m[i]
        τᵢ = Tc[i]*Tinv
        δrefᵢ = Vc[i]*v̄ref
        Δa += z[i]*(reduced_a_res(mᵢ,δref,τ,lnδref,lnτ) - reduced_a_res(mᵢ,δrefᵢ,τᵢ))
    end
    Δa *= ∑z⁻¹
    R = Rgas(model)
    ρ = ∑z/V
    lnb = log1p(b*ρ)/log1p(b*ρref)
    return aᵣ + lnb*(gᴱ/(R*T) - Δa)
end

function lb_volume(model::EmpiricMultiFluid{A,M,GEDeparture},z = SA[1.0]) where {A,M}
    vref = model.departure.vref
    v̄ref = dot(z,vref)/sum(z)
    return 0.8547008547008548*v̄ref
end
