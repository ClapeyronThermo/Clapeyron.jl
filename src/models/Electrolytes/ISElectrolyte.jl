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

function ISElectrolyteWrapper(model::ESElectrolyteModel;salts = nothing)
    salt = SaltParam(model,salts)
    components = salt.implicit_components
    return ISElectrolyteWrapper(components,model,salt)
end

function a_res(model::ISElectrolyteWrapper, V, T, z)
    w = to_ion(model.salt,z)
    return a_res(model.model,V,T,w)
end

#=
function a_res(model::ESElectrolyteWrapper, V, T, z)
    w = to_salt(model.salt,z)
    return a_res(model.model,V,T,w)
end =#

function idealmodel(model::ISElectrolyteWrapper)
    return ISElectrolyteIdealWrapper(model.components,idealmodel(model.model),model.salt)
end

function a_ideal(model::ISElectrolyteIdealWrapper,V,T,z)
    w = to_ion(model.salt,z)
    return a_ideal(model.model,V,T,w)
end

Rgas(model::ISElectrolyteWrapper) = Rgas(model.model)
Rgas(model::ISElectrolyteIdealWrapper) = Rgas(model.model)

#=
function eos_impl(model::ISElectrolyteWrapper,V,T,z)
    w = to_ion(model.salt,z)
    return eos_impl(model.model,V,T,w)
end=#

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
    w = to_ion(model.salt,z)
    return volume(model.model, p, T, w, phase=phase, threaded=threaded, vol0=vol0)
end

function lb_volume(model::ISElectrolyteWrapper,T,z)
    w = to_ion(model.salt,z)
    return lb_volume(model.model,T,w)
end

function mw(model::ISElectrolyteWrapper)
    return mw(model.neutralmodel)
end

function p_scale(model::ISElectrolyteWrapper,z)
    w = to_ion(model.salt,z)
    return p_scale(model.model,w)
end

function T_scale(model::ISElectrolyteWrapper,z)
    w = to_ion(model.salt,z)
    return T_scale(model.model,w)
end

export ISElectrolyteWrapper
