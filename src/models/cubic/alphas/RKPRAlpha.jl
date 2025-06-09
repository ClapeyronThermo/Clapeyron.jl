abstract type RKPRAlphaModel <: AlphaModel end

struct RKPRAlphaParam <: EoSParam
    k1::SingleParam{Float64}
    k2::SingleParam{Float64}
    acentricfactor::SingleParam{Float64}
end

@newmodelsimple RKPRAlpha RKPRAlphaModel RKPRAlphaParam
export RKPRAlpha

"""
    RKPRAlpha <: RKPRAlphaModel
    
    RKPRAlpha(components;
    userlocations = String[],
    verbose::Bool=false)

## Input Parameters

- `k1`: Single Parameter (`Float64`)
- `k1`: Single Parameter (`Float64`)
- `acentricfactor`: Single Parameter (`Float64`) (optional), acentric factor. used to estimate the value of `k1` and `k2` when not provided.

## Description

Cubic alpha `(α(T))` model. Default for [`RKPR`](@ref) EoS.
```
αᵢ = ((k₂ᵢ + 1)/(k₂ᵢ + Trᵢ))^k₁ᵢ
```

in the case of RKPR EoS. if no values of `k1` and `k2` are provided, the following correlation is used:
```
k₁ = (12.504Zc -2.7238) + (7.4513*Zc + 1.9681)ωᵢ + (-2.4407*Zc + 0.0017)ωᵢ^2
Trᵢ = T/Tcᵢ
k₂ = 2
```
for Other cubic EoS, the k is calculated via the Leivobici alpha:
```
k₂ = 2
k₁ = log(αᵢ(leivovici,Tr = 0.7))/log(3/2.7)
```
## Model Construction Examples
```
# Using the default database
alpha = RKPRAlpha("water") #single input
alpha = RKPRAlpha(["water","ethanol"]) #multiple components

# Using user-provided parameters

# Passing files or folders
alpha = RKPRAlpha(["neon","hydrogen"]; userlocations = ["path/to/my/db","critical/acentric.csv"])

# Passing parameters directly
alpha = RKPRAlpha(["neon","hydrogen"];userlocations = (;acentricfactor = [-0.03,-0.21]))
```
"""
RKPRAlpha
default_locations(::Type{RKPRAlpha}) = critical_data()

function transform_params(::Type{RKPRAlphaParam},params,components)
    n = length(components)
    k1 = get!(params,"k1") do
        SingleParam("k1",components,zeros(n),fill(true,n))
    end

    k2 = get!(params,"k2") do
        SingleParam("k2",components,zeros(n),fill(true,n))
    end
    return params
end

fast_build_alpha(::Type{RKPRAlpha}) = true


function α_function(model::CubicModel,V,T,z,alpha_model::RKPRAlphaModel)
    k1 = alpha_model.params.k1.values
    k2 = alpha_model.params.k2.values
    Tc = model.params.Tc.values
    α = zeros(typeof(1.0*T),length(Tc))
    for i in @comps
        Tr = T/Tc[i]
        k2i = k2[i]
        k1i = k1[i]
        α[i] = ((k2i + 1)/(k2i + Tr))^k1i
    end
    return α
end

function α_function(model::CubicModel,V,T,z::SingleComp,alpha_model::RKPRAlphaModel)
    k1 = alpha_model.params.k1.values[1]
    k2 = alpha_model.params.k2.values[1]
    Tc = model.params.Tc.values[1]
    Tr = T/Tc
    α = ((k2 + 1)/(k2 + Tr))^k1
    return SA[α]
end

function recombine_alpha!(model::DeltaCubicModel,alpha::RKPRAlphaModel)
    ω = alpha.params.acentricfactor.values
    ω = alpha.params.acentricfactor.values
    k1 = alpha.params.k1
    k2 = alpha.params.k2
    for i in 1:length(model)
        if k2.ismissingvalues[i]
            k2[i] = 2.0
            k2i = 2.0
        else
            k2i = k2[i]
        end

        if k1.ismissingvalues[i]
            poly = α_m_leibovici(model,i)
            mi = evalpoly(ω[i],poly)
            αi = (1+mi*(1-sqrt(0.7)))^2
            k1[i] = log(αi)/log((k2i + 1)/(k2i + 0.7))
        end
    end
end

function recombine_alpha!(model::RKPRModel,alpha::RKPRAlphaModel)
    Pc = model.params.Pc.values
    Tc = model.params.Tc.values
    ω = alpha.params.acentricfactor.values
    k1 = alpha.params.k1
    k2 = alpha.params.k2
    for i in 1:length(model)
        if k2.ismissingvalues[i]
            k2[i] = 2.0
            k2i = 2.0
        else
            k2i = k2[i]
        end

        if k1.ismissingvalues[i] && k2i == 2.0
            z = FillArrays.OneElement(i, length(model))
            Δ1,Δ2 = cubic_Δ(model,z)
            b = lb_volume(model,Tc[i],z)
            Δ1,Δ2 = cubic_Δ(model,z)
            B = b*Pc[i]/(Rgas(model)*Tc[i])
            Zc = (1 + (Δ1 + Δ2 + 1)*B)/3 #Pc
            k1[i] = evalpoly(ω[i],(12.504*Zc -2.7238,7.4513*Zc + 1.9681,-2.4407*Zc + 0.0017))
        elseif k1.ismissingvalues[i] && k2i != 2.0
            poly = α_m_leibovici(model,i)
            mi = evalpoly(ω[i],poly)
            αi = (1+mi*(1-sqrt(0.7)))^2
            k1[i] = log(αi)/log((k2i + 1)/(k2i + 0.7))
        end
    end
end