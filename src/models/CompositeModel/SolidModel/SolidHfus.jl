abstract type SolidHfusModel <: EoSModel end

struct SolidHfusParam <: EoSParam
    Hfus::SingleParam{Float64}
    Tm::SingleParam{Float64}
end

@newmodelsimple SolidHfus SolidHfusModel SolidHfusParam

SolidHfus
default_locations(::Type{SolidHfus}) = ["solids/fusion.csv"]
default_references(::Type{SolidHfus}) = String[]

function chemical_potential(model::SolidHfusModel, p, T, z)
    Hfus = model.params.Hfus.values
    Tm = model.params.Tm.values
    return @. Hfus*T*(1/Tm-1/T)
end

export SolidHfus