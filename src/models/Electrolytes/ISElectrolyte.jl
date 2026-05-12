"""
    MeanIonicApproach(model::ESElectrolyteModel;salts = nothing)


Given an explicit solvent model, returns an implicit solvent model, where the charged components are paired to form binary salts.
If no `salts` argument is specified, the salt pairings will be created via [`Clapeyron.auto_binary_salts`](@ref).
## Example

´´´julia-repl
julia> system = ePCSAFT(["water","acetonitrile"],["sodium","chloride"])
Explicit Electrolyte Model with 4 components:
 "water"
 "acetonitrile"
 "sodium" (+1)
 "chloride" (-1)
Neutral Model: pharmaPCSAFT{BasicIdeal, Float64}
Ion Model: hsdDH{ConstRSP}
RSP Model: ConstRSP

julia> salt_system = MeanIonicApproach(system)
MeanIonicApproach{ePCSAFT{BasicIdeal, pharmaPCSAFT{BasicIdeal, Float64}, hsdDH{ConstRSP}}} with 3 components:
 "water"
 "acetonitrile"
 "sodium.chloride"
´´´

"""
struct MeanIonicApproach{M} <: ISElectrolyteModel
    components::Vector{String}
    model::M
    salt::SaltParam
end

struct ISElectrolyteIdealWrapper{M} <: IdealModel
    components::Vector{String}
    model::M
    salt::SaltParam
end

salt_compositions(m::Union{ISElectrolyteIdealWrapper,MeanIonicApproach},z) = salt_compositions(m.salt,z) 
ion_compositions(m::Union{ISElectrolyteIdealWrapper,MeanIonicApproach},z) = ion_compositions(m.salt,z) 

function MeanIonicApproach(model::ESElectrolyteModel;salts = nothing)
    salt = SaltParam(model,salts)
    components = salt.implicit_components
    return MeanIonicApproach(components,model,salt)
end

function PT_property(model::MeanIonicApproach,p,T,z,phase,threaded,vol0,f::F,vol::V) where {F,V}
    w = ion_compositions(model,z)
    PT_property(model.model,p,T,w,phase,threaded,vol0,f,vol)
end

function VT_pressure(model::MeanIonicApproach,V,T,z)
    w = ion_compositions(model,z)
    return VT_pressure(model.model,V,T,z)
end

#=
function a_res(model::MeanIonicApproach, V, T, z)
    w = ion_compositions(model,z)
    return a_res(model.model,V,T,w)
end
=#

Rgas(model::MeanIonicApproach) = Rgas(model.model)
Rgas(model::ISElectrolyteIdealWrapper) = Rgas(model.model)

reference_state(model::MeanIonicApproach) = reference_state(model.model)

function tp_flash_K0!(K,model::ISElectrolyteModel,p,T,z)
    neutral = ones(Bool,length(model))
    isalts = model.salt.isalts
    neutral[isalts] .= false
    pures = split_model(model,neutral)
    psat = first.(extended_saturation_pressure.(pures,T))
    K .= 0
    Kview = @view K[neutral]
    Kview .= psat ./ p
    return K
end

function each_split_model(model::MeanIonicApproach,I_salt)
    salt_i,I_ion = IS_each_split_model(model.salt,I_salt)
    return MeanIonicApproach(model.components[I_salt],each_split_model(model.model,I_ion),salt_i)
end

function volume_impl(model::MeanIonicApproach, p, T, z, phase, threaded, vol0)
    w = ion_compositions(model,z)
    return volume(model.model, p, T, w, phase=phase, threaded=threaded, vol0=vol0)
end

function mw(model::MeanIonicApproach)
    return mw(model.neutralmodel)
end

function ∂lnϕ_cache(model::MeanIonicApproach, p, T, z, BB::Val{B}) where B
    w = ion_compositions(model,z)
    return ∂lnϕ_cache(model.model, p, T, w, BB)
end

function tpd_lnϕ_and_v!(cache,wrapper::MeanIonicApproach,p,T,w,vol0,liquid_overpressure = false,phase = :liquid,_vol = nothing)
    ww = ion_compositions(wrapper,w)
    ∑w = sum(ww)
    lnϕw,v,overpressure = tpd_lnϕ_and_v!(cache,wrapper.model,p,T,ww,vol0,liquid_overpressure,phase,_vol)
    
    salt = wrapper.salt
    if iszero(length(salt.isalts))
        return lnϕw,v,overpressure
    end
    lnϕw .+= log.(ww)
    lnϕw .-= log(∑w/sum(w)) #does normalize by sum(w) matter?
    
    Z = wrapper.model.charge
    lnϕz = similar(lnϕw,length(w))
    E = eachcol(salt.E)
    for i in 1:length(w)
        Ei = E[i]
        if i in salt.isalts
            iref = findfirst(Ei)
            lnϕref = lnϕw[iref]
            lnϕww = @view lnϕw[Ei]
            ZZ = @view Z[Ei]
            lnϕ1,lnϕ2 = lnϕww
            Z1,Z2 = ZZ
            if Z1 > 0
                lnϕz[i] = lnϕ1 + Z1/(Z2 - Z1) * (lnϕ1 - lnϕ2)
            else
                lnϕz[i] = lnϕ2 + Z2/(Z1 - Z2) * (lnϕ2 - lnϕ1)
            end
        else
            lnϕz[i] = only(@view(lnϕw[Ei]))
        end
    end
    lnϕz -= log.(w)

    return lnϕz,v,overpressure
end

function modified_lnϕ(wrapper::MeanIonicApproach, p, T, z, cache; phase = :unknown, vol0 = nothing)
    lnϕz,v,_ = tpd_lnϕ_and_v!(cache,wrapper,p,T,z,vol0,false,phase,nothing)
    return lnϕz,v
end

function tpd_input_composition(wrapper::MeanIonicApproach,p,T,z,lle,cache = tpd_cache(wrapper,p,T,z,z))
    d_l,d_v,_,_,_,Hϕ = cache
    TT = Base.promote_eltype(wrapper.model,p,T,z)
    w = ion_compositions(wrapper,z)
    n = sum(w)
    logsumz = log(n)
    d,vl = tpd_lnϕ_and_v!(Hϕ,wrapper,p,T,z,nothing,false,:liquid)
    d_l .= d
    d_l .+= log.(z) .- logsumz
    
    lle && return copy(d_l),:liquid,vl
    d,vv = tpd_lnϕ_and_v!(Hϕ,wrapper,p,T,z,nothing,false,:vapour)
    d_v .= d
    d_v .+= log.(z) .- logsumz
    gr_l = dot(z,d_l)
    gr_v = dot(z,d_v)
    isnan(gr_v) && (gr_v = Inf*one(gr_v))
    if gr_l < gr_v
        return copy(d_l),:liquid,vl
    else
        return copy(d_v),:vapour,vv
    end
end

function ∂lnϕ∂P(wrapper::MeanIonicApproach, p, T, z=SA[1.], cache = ∂lnϕ_cache(wrapper,p,T,z,Val{false}());
            phase=:unknown,
            vol0=nothing,
            threaded = true,
            vol = volume(wrapper,p,T,z;phase,vol0,threaded))

    ∂lnϕ∂Pi = ∂lnϕ∂P(wrapper.model,p,T,w,cache;phase,vol0,threaded,vol)
    E = eachrow(wrapper.salt.E)
    for i in 1:length(∂lnϕ∂Pi)
        Ei = E[i]
        ∂lnϕ∂Pi[i] = dot(Ei,∂lnϕ∂Pi)/sum(Ei)
    end
    return ∂lnϕ∂Pi, V
end


function ∂lnϕ∂T(wrapper::MeanIonicApproach, p, T, z=SA[1.], cache = ∂lnϕ_cache(wrapper,p,T,z,Val{true}());
            phase=:unknown,
            vol0=nothing,
            threaded = true,
            vol = volume(wrapper,p,T,z;phase,vol0,threaded))

    ∂lnϕ∂Ti = ∂lnϕ∂Ti(wrapper.model,p,T,w,cache;phase,vol0,threaded,vol)
    E = eachrow(wrapper.salt.E)
    for i in 1:length(∂lnϕ∂Ti)
        Ei = E[i]
        ∂lnϕ∂Ti[i] = dot(Ei,∂lnϕ∂Ti)/sum(Ei)
    end
    return ∂lnϕ∂Ti, V
end

export MeanIonicApproach
