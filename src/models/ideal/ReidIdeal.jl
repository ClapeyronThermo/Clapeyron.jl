abstract type PolynomialIdealModel <: IdealModel end

struct ReidIdealParam <: EoSParam
    coeffs::SingleParam{NTuple{5,Float64}}
    reference_state::ReferenceState
end

function reid_coeffs(a::SingleParameter,b::SingleParameter,c::SingleParameter,d::SingleParameter,e::SingleParameter,comps::Vector{String})
    comps = a.components
    a = a.values
    b = b.values
    c = c.values
    d = d.values
    _e = e == nothing ? FillArrays.Zeros(length(comps)) : e.values

    return reid_coeffs(a,b,c,d,_e,comps)
end

function reid_coeffs(a::AbstractVector,b::AbstractVector,c::AbstractVector,d::AbstractVector,comps::Vector{String})
    n = length(a)
    coeffs = [(a[i],b[i],c[i],d[i],zero(a[i])) for i in 1:n]
    SingleParam("Reid Coefficients",comps,coeffs)
end

function reid_coeffs(a::AbstractVector,b::AbstractVector,c::AbstractVector,d::AbstractVector,e::AbstractVector,comps::Vector{String})
    n = length(a)
    coeffs = [(a[i],b[i],c[i],d[i],e[i]) for i in 1:n]
    SingleParam("Reid Coefficients",comps,coeffs)
end

abstract type ReidIdealModel <: PolynomialIdealModel end
@newmodelsimple ReidIdeal ReidIdealModel ReidIdealParam

"""
    ReidIdeal <: IdealModel

    ReidIdeal(components; 
    userlocations = String[],
    reference_state = nothing,
    verbose = false)

## Input parameters

- `a`: Single Parameter (`Float64`) - polynomial coefficient
- `b`: Single Parameter (`Float64`) - polynomial coefficient
- `c`: Single Parameter (`Float64`) - polynomial coefficient
- `d`: Single Parameter (`Float64`) - polynomial coefficient
- `e`: Single Parameter (optional) (`Float64`)  - polynomial coefficient

## Model parameters

- `coeffs`: Single Parameter (`NTuple{5,Float64}`)

## Description

Reid Ideal Model. Helmholtz energy obtained via integration of specific heat capacity:

```
Cpᵢ(T) = aᵢ  + bᵢT + cᵢT^2 + dᵢT^3 + eᵢT^4
Cp(T) = ∑Cpᵢxᵢ
```

## Model Construction Examples
```
# Using the default database
idealmodel = ReidIdeal("water") #single input
idealmodel = ReidIdeal(["water","ethanol"]) #multiple components

# Using user-provided parameters

# Passing files or folders
idealmodel = ReidIdeal(["neon","hydrogen"]; userlocations = ["path/to/my/db","reid.csv"])

# Passing parameters directly
idealmodel = ReidIdeal(["water","butane"];
            userlocations = (a = [32.24, 9.487], 
                        b = [0.00192, 0.3313], 
                        c = [1.06e-5, -0.0001108],
                        d = [-3.6e-9, -2.822e-9])
                        ) #e is not used
```

"""
ReidIdeal

export ReidIdeal
default_locations(::Type{ReidIdeal}) = ["ideal/ReidIdeal.csv"]
default_ignore_missing_singleparams(::Type{ReidIdeal}) = ["e"]
function transform_params(::Type{ReidIdeal},params,components)
    a,b,c,d = params["a"],params["b"],params["c"],params["d"]
    e = get(params,"e") do
        SingleParam("e",components)
    end
    params["coeffs"] = reid_coeffs(a,b,c,d,e,components)
    return params
end

recombine_impl!(model::ReidIdealModel) = model

evalcoeff(::ReidIdealModel,coeffs,T,lnT = log(T)) = evalpoly(T,coeffs)

function eval∫coeff(::ReidIdealModel,coeffs,T,lnT = log(T))
    n = length(coeffs)
    div1 = NTuple{n,Int}(1:n)
    ∫poly = coeffs ./ div1
    return evalpoly(T,∫poly)*T
end

function eval∫coeffT(::ReidIdealModel,coeffs,T,lnT = log(T))
    n = length(coeffs)
    div1 = NTuple{n-1,Int}(1:(n-1))
    A = first(coeffs)
    coeffs1 = coeffs[2:end]
    ∫polyT = coeffs1 ./ div1
    return evalpoly(T,∫polyT)*T + A*lnT
end

function a_ideal(model::PolynomialIdealModel, V, T, z)
    #x = z/sum(z)
    polycoeff = model.params.coeffs.values
    #return sum(x[i]*(log(z[i]/V) + 1/(R̄*T)*(sum(polycoeff[k][i]/k*(T^k-298^k) for k in 1:4)) -
    #    1/R̄*((polycoeff[k][1]-R̄)*log(T/298)+sum(polycoeff[k][i]/(k-1)*(T^(k-1)-298^(k-1)) for k in 2:4))) for i in @comps)
    V⁻¹ = 1/V
    res = zero(V+T+first(z))
    Σz = sum(z)
    RT = R̄*T
    R̄⁻¹ = 1/R̄
    RT⁻¹ = 1/RT
    T0 = 298.
    lnT0 = log(T0)
    lnT = log(T)
    @inbounds for i in @comps
        coeffs = polycoeff[i] 
        H = (eval∫coeff(model,coeffs,T,lnT) - eval∫coeff(model,coeffs,T0,lnT0))*RT⁻¹
        TS = (eval∫coeffT(model,coeffs,T,lnT) - eval∫coeffT(model,coeffs,T0,lnT0))*R̄⁻¹
        α₀ᵢ = H - TS + lnT - lnT0
        res += z[i]*α₀ᵢ
        res += xlogx(z[i],V⁻¹)
    end
    return res/Σz
end

function ∂²f∂T²(model::PolynomialIdealModel,V,T,z)
    coeff = model.params.coeffs.values
    Cp = zero(T+first(z))
    Σz = sum(z)
    for i in @comps
        pol = coeff[i]
        Cp +=z[i]*evalcoeff(model,pol,T)
    end
    Cv = Cp - Σz*Rgas(model)
    return -Cv/T
end