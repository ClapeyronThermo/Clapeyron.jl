abstract type SolidKsModel <: GibbsBasedModel end

struct SolidKsParam <: EoSParam
    Gform::SingleParam{Float64}
    Hform::SingleParam{Float64}
    Tref::SingleParam{Float64}
end

@newmodelsimple SolidKs SolidKsModel SolidKsParam

"""
    SolidKsModel <: EoSModel

    SolidKs(components;
    userlocations = String[],
    verbose::Bool=false)

## Parameters

  - `Gform`: Single Parameter (`Float64`) - Formation Gibbs energy in water at infinite dilution 1 bar and the reference temperature`[J·mol⁻¹]`
  - `Hform`: Single Parameter (`Float64`) -Formation enthalpy in water at infinite dilution 1 bar and the reference temperature`[J·mol⁻¹]`
  - `Tref`: Single Parameter (`Float64`) - Reference temperature `[K]`

## Description

Approximation of the excess chemical potential in the solid phase, using enthalpies and gibbs energies of formation where the liquid reference is at infinite dilution in water:

```
ln(xᵢγᵢ) = -Gformᵢ*T/Trefᵢ - Hformᵢ*(1 - T/Trefᵢ)
```
"""
SolidKs
default_locations(::Type{SolidKs}) = ["solids/formation.csv"]
default_references(::Type{SolidKs}) = String[]
default_ignore_missing_singleparams(::Type{SolidKs}) = ["CpSL"]

function volume_impl(model::SolidKsModel, p, T, z, phase, threaded, vol0)
    _0 = zero(T + first(z))
    return _0/_0
end

function chemical_potential_impl(model::SolidKsModel, p, T, z, phase, threaded, vol0)
    Gform = model.params.Gform.values
    Hform = model.params.Hform.values
    Tref = model.params.Tref.values
    return @. -Gform*T/Tref - Hform*(1 - T/Tref)
end

function eos_g(model::SolidKsModel, p, T, z)
    Gform = model.params.Gform.values
    Hform = model.params.Hform.values
    Tref = model.params.Tref.values
    g = zero(Base.promote_eltype(model, T, z))
    for i in 1:length(model)
        Trefi = Tref[i]
        μi = -Gform[i]*T/Trefi - Hform[i]*(1 - T/Trefi)
        g += z[i]*μi
    end
    return g
end

sle_Δμ(solid::SolidKsModel, model, p, T, j) = sle_pure_solid_μ(solid, p, T, j)/(Rgas(model.fluid)*T)

#formation properties refer to the ions at infinite dilution in water
function sle_μ_ref(solid::SolidKsModel, fluid::ESElectrolyteModel, p, T)
    zref = fill(1e-30, length(fluid))
    zref[find_water_indx(fluid)] = 1.0
    zref ./= sum(zref)
    return chemical_potential(fluid, p, T, zref; phase=:l) .- Rgas(fluid)*T .* log.(zref)
end

export SolidKs
