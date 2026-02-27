"""
    ISElectrolyteWrapper(model::ESElectrolyteModel;salts = nothing)


Given en explicit solvent model, returns an implicit solvent model, where are the charged components are paired to form binary salts.
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

julia> salt_system = ISElectrolyteWrapper(system)
ISElectrolyteWrapper{ePCSAFT{BasicIdeal, pharmaPCSAFT{BasicIdeal, Float64}, hsdDH{ConstRSP}}} with 3 components:
 "water"
 "acetonitrile"
 "sodium.chloride"
´´´

"""
struct ISElectrolyteWrapper{M} <: ISElectrolyteModel
    components::Vector{String}
    model::M
    salt::SaltParam
end

struct ISElectrolyteIdealWrapper{M} <: IdealModel
    components::Vector{String}
    model::M
    salt::SaltParam
end

salt_compositions(m::Union{ISElectrolyteIdealWrapper,ISElectrolyteWrapper},z) = salt_compositions(m.salt,z) 
ion_compositions(m::Union{ISElectrolyteIdealWrapper,ISElectrolyteWrapper},z) = ion_compositions(m.salt,z) 

function ISElectrolyteWrapper(model::ESElectrolyteModel;salts = nothing)
    salt = SaltParam(model,salts)
    components = salt.implicit_components
    return ISElectrolyteWrapper(components,model,salt)
end

function PT_property(model::ISElectrolyteWrapper,p,T,z,phase,threaded,vol0,f::F,USEP::Val{UseP}) where {F,UseP}
    w = ion_compositions(model,z)
    PT_property(model.model,p,T,w,phase,threaded,vol0,f,USEP)
end

function VT_pressure(model::ISElectrolyteWrapper,V,T,z)
    w = ion_compositions(model,z)
    return VT_pressure(model.model,V,T,z)
end



#=
function a_res(model::ISElectrolyteWrapper, V, T, z)
    w = ion_compositions(model,z)
    return a_res(model.model,V,T,w)
end
=#

Rgas(model::ISElectrolyteWrapper) = Rgas(model.model)
Rgas(model::ISElectrolyteIdealWrapper) = Rgas(model.model)

reference_state(model::ISElectrolyteWrapper) = reference_state(model.model)

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

function each_split_model(model::ISElectrolyteWrapper,I_salt)
    salt_i,I_ion = IS_each_split_model(model.salt,I_salt)
    return ISElectrolyteWrapper(model.components[I_salt],each_split_model(model.model,I_ion),salt_i)
end

function volume_impl(model::ISElectrolyteWrapper, p, T, z, phase, threaded, vol0)
    w = ion_compositions(model,z)
    return volume(model.model, p, T, w, phase=phase, threaded=threaded, vol0=vol0)
end

function mw(model::ISElectrolyteWrapper)
    return mw(model.neutralmodel)
end

function ∂lnϕ_cache(model::ISElectrolyteWrapper, p, T, z, BB::Val{B}) where B
    w = ion_compositions(model,z)
    return ∂lnϕ_cache(model.model, p, T, w, BB)
end

function tpd_lnϕ_and_v!(cache,wrapper::ISElectrolyteWrapper,p,T,w,vol0,liquid_overpressure = false,phase = :liquid,_vol = nothing)
    lnϕw,v,overpressure = tpd_lnϕ_and_v!(cache,wrapper.model,p,T,w,vol0,false,phase,_vol)
    if iszero(count(!iszero,wrapper.model.charge))
        return lnϕw,v,overpressure
    end
    lnϕz = similar(lnϕw,length(lnϕw) - 1)
    idx = zeros(Bool,length(lnϕz))
    salt = wrapper.salt
    nions = length(wrapper.model)
    E = eachrow(salt.E)
    for i in 1:length(lnϕz)
        Ei = E[i]
        lnϕz[i] = dot(Ei,lnϕw)/sum(Ei)
    end
    return lnϕz,v,overpressure
end


export ISElectrolyteWrapper
