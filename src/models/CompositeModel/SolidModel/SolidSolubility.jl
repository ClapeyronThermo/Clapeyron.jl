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
    return @. Hfus*T*(1/Tm-1/T)
end

export SolidSolubility