


struct ISElectrolite{T<:IdealModel,c<:EoSModel,i<:IonModel} <: ISElectrolyteModel
    components::Array{String,1}
    idealmodel::T
    neutralmodel::c
    ionmodel::i
    references::Array{String,1}
end


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

function ISElectrolyteWrapper(model::ESElectrolyteModel)
    salt = SaltParam(model)
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

function eos_impl(model::ISElectrolyteWrapper,V,T,z)
    w = to_ion(model.salt,z)
    ∑z = sum(z)
    return ∑z*Rgas(model)*T*eos_impl(model,V,T,w)
end

function tp_flash_K0!(K,model::ISElectrolyteModel,p,T,z)
    neutral = zeros(Bool,length(model))
    r = eachrow(model.salt.mat)
    for i in 1:length(model)
        if count(!iszero,r[i]) == 1
            neutral[i] = true
        end
    end
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

export ISElectrolyteWrapper
