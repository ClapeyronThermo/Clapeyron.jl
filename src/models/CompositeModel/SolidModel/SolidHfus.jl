abstract type SolidHfusModel <: GibbsBasedModel end

struct SolidHfusParam <: EoSParam
    Hfus::SingleParam{Float64}
    Tm::SingleParam{Float64}
    CpSL::SingleParam{Float64}
end

@newmodelsimple SolidHfus SolidHfusModel SolidHfusParam

"""
    SolidHfusModel <: EoSModel

    SolidHfus(components;
    userlocations = String[],
    verbose::Bool=false)

## Parameters

- `Hfus`: Single Parameter (`Float64`) - Enthalpy of Fusion at 1 bar `[J·mol⁻¹]`
- `Tm`: Single Parameter (`Float64`) - Melting Temperature `[K]`
- `CpSL`: Single Parameter (`Float64`) (optional) - Heat Capacity of the Solid-Liquid Phase Transition `[J·mol⁻¹·K⁻¹]`

## Description

Approximation of the excess chemical potential in the solid phase (`CpSL` is not necessary by default):
```
ln(xᵢγᵢ) = Hfusᵢ*T*(1/Tmᵢ-1/T)-CpSLᵢ/R̄*(Tmᵢ/T-1-log(Tmᵢ/T))
```
"""
SolidHfus
default_locations(::Type{SolidHfus}) = ["solids/fusion.csv"]
default_references(::Type{SolidHfus}) = String[]
default_ignore_missing_singleparams(::Type{SolidHfus}) = ["CpSL"]

sle_T_ref(model::SolidHfusModel) = model.params.Tm.values

function chemical_potential_impl(model::SolidHfusModel,p,T,z,phase,threaded,vol0)
    Hfus = model.params.Hfus.values
    Tm = model.params.Tm.values
    CpSL = model.params.CpSL.values
    return @. Hfus*T*(1/Tm-1/T)-CpSL/Rgas()*(Tm/T-1-log(Tm/T))
end

p_scale(model::SolidHfusModel,z) = 101325.0
T_scale(model::SolidHfusModel,z) = dot(model.params.Tm.values,z)/sum(z) 

function eos_g(model::SolidHfusModel,p,T,z)
    Hfus = model.params.Hfus.values
    Tm = model.params.Tm.values
    CpSL = model.params.CpSL.values
    g = zero(Base.promote_eltype(model,T,z))
    for i in 1:length(model)
        Tmi = Tm[i]
        μi = Hfus[i]*T*(1/Tmi-1/T)-CpSL[i]/Rgas(model)*(Tmi/T-1-log(Tmi/T))
        g += z[i]*μi
    end
    return g
end

function init_preferred_method(method::typeof(melting_pressure),model::CompositeModel{<:GibbsBasedModel,<:SolidHfusModel},kwargs...)
    return MeltingCorrelation()
end

function init_preferred_method(method::typeof(melting_pressure),model::SolidHfusModel,kwargs...)
    return MeltingCorrelation()
end

function melting_pressure_impl(model::CompositeModel{<:EoSModel,<:SolidHfusModel},T,method::MeltingCorrelation)
    P,vs,_ = melting_pressure_impl(model.solid,T,method)
    vl = volume(model.fluid,P,T,phase = :l)
    return P,vs,vl
end

function melting_pressure_impl(model::SolidHfusModel,T,method::MeltingCorrelation)
    Hfus = model.params.Hfus.values[1]
    Tm = model.params.Tm.values[1]
    Pm = 1e5
    logP = log(Pm) - Hfus*(1/T - 1/Tm)/Rgas()
    P = exp(logP)
    #=
    dPdT = P*dlogPdT = P*Hfus/T^2/R
    =#
    nan = zero(P)/zero(P)
    return P, nan, nan
end

function init_preferred_method(method::typeof(melting_temperature),model::CompositeModel{<:GibbsBasedModel,<:SolidHfusModel},kwargs...)
    return MeltingCorrelation()
end

function init_preferred_method(method::typeof(melting_temperature),model::SolidHfusModel,kwargs...)
    return MeltingCorrelation()
end


function melting_temperature(model::CompositeModel{<:Any,<:SolidHfusModel},P,method::ThermodynamicMethod)
    return melting_temperature_impl(model,P,method)
end

function melting_temperature_impl(model::CompositeModel{<:EoSModel,<:SolidHfusModel},P,method::MeltingCorrelation)
    T,vs,_ = melting_temperature_impl(model.solid,P,method)
    vl = volume(fluid_model(model),P,T,SA[1.0],phase = :l)
    return T,vs,vl
end

function melting_temperature_impl(model::SolidHfusModel,P,method::MeltingCorrelation)
    Hfus = model.params.Hfus.values[1]
    Tm = model.params.Tm.values[1]
    Pm = 1e5
    T = 1/(1/Tm - log(P/Pm)*Rgas(model)/Hfus)
    nan = zero(T)/zero(T)
    return T, nan, nan
end

function x0_eutectic_point(model::CompositeModel{<:EoSModel,<:SolidHfusModel},p)
    return x0_eutectic_point(model.solid,p)
end

function x0_eutectic_point(model::SolidHfusModel,p)
    Hfus = model.params.Hfus.values
    Tm = model.params.Tm.values
    #Hfus correlation
    Tmax = 2/minimum(Tm)
    R = Rgas()
    f(tinv) = 1 - exp(-Hfus[1]/R*(tinv -1/Tm[1])) - exp(-Hfus[2]/R*(tinv -1/Tm[2]))
    prob = Roots.ZeroProblem(f,(zero(Tmax),Tmax))
    T0inv = Roots.solve(prob)
    T0 = 1/T0inv
    x0 = exp(-Hfus[1]/Rgas()*(1/T0-1/Tm[1]))
    return [T0/200.,x0]
end

export SolidHfus
