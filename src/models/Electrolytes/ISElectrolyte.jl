struct ISElectrolyteWrapper{M} <: ISElectrolyteModel
    components::Vector{String}
    model::M
    salt::SaltParam
end

struct ESElectrolyteWrapper{M} <: ESElectrolyteModel
    components::Vector{String}
    charge::Vector{Int64}
    model::M
    salt::SaltParam
end

struct ISElectrolite{T<:IdealModel,c<:EoSModel,i<:IonModel} <: ISElectrolyteModel
    components::Array{String,1}
    idealmodel::T
    neutralmodel::c
    ionmodel::i
    references::Array{String,1}
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

function a_res(model::ESElectrolyteWrapper, V, T, z)
    w = to_salt(model.salt,z)
    return a_res(model.model,V,T,w)
end

idealmodel(model::ISElectrolyteWrapper) = idealmodel(model.model)
idealmodel(model::ESElectrolyteWrapper) = idealmodel(model.model)

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

export ISElectrolyteWrapper