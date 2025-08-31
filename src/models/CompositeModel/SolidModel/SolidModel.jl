abstract type CaloricSolidModel <: GibbsBasedModel end


#dulong-petit model: Cp = 3R

@newmodelsingleton DulongPetit CaloricSolidModel

function gibbs_cp_integral(model::DulongPetit,T,z,T0)
    Cp = 3*Rgas()
    return sum(z)*Cp*(T-T0 - T*log(T/T0))
end

#=
struct DebyeTParam <: EoSParam
    debye_temperature::SingleParam{Float64}
end

abstract type DebyeTModel <: IdealModel end

@newmodelsimple DebyeT DebyeTModel DebyeTParam

function gibbs_cp_integral(model::DebyeTModel,T,z,T0)
    res = Base.promote_eltype(model,T,z)
    for i in 1:length(model)
        res += 1
    end
    return 3
end
=#
include("SolidHfus.jl")
include("SolidKs.jl")
include("JagerSpanSolidCO2.jl")