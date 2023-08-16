abstract type SolidSolubilityModel <: EoSModel end

struct SolidSolubilityParam <: EoSParam
    Hfus::SingleParam{Float64}
    Tm::SingleParam{Float64}
end

@newmodelsimple SolidSolubility SolidSolubilityModel SolidSolubilityParam

SolidSolubility
default_locations(::Type{SolidSolubility}) = String[]
default_references(::Type{SolidSolubility}) = String[]

function chemical_potential(model::SolidSolubilityModel, p, T, z)
    Hfus = model.params.Hfus.values
    Tm = model.params.Tm.values
    μ = zeros(length(z))
    for i in @comps
        if !isnothing(Tm[i])
            μ[i] = -Hfus[i]/(Rgas()*Tm[i])*(1-T/Tm[i])
        end
    end
    return μ
end

export SolidSolubility