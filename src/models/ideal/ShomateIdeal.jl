

abstract type ShomateIdealModel <: ReidIdealModel end
@newmodelsimple ShomateIdeal ShomateIdealModel ReidIdealParam

"""
    ShomateIdeal <: ShomateIdealModel
    ShomateIdeal(components; 
    userlocations::Array{String,1}=String[], 
    verbose = false)

## Input parameters

- `a`: Single Parameter (`Float64`) - polynomial coefficient
- `b`: Single Parameter (`Float64`) - polynomial coefficient
- `c`: Single Parameter (`Float64`) - polynomial coefficient
- `d`: Single Parameter (`Float64`) - polynomial coefficient
- `e`: Single Parameter (optional) (`Float64`)  - polynomial coefficient for 1/T^2

## Model parameters

- `coeffs`: Single Parameter (`NTuple{5,Float64}`)

## Description

Shomate Ideal Model. Helmholtz energy obtained via integration of specific heat capacity:

```
Cpᵢ(T) = aᵢ  + bᵢT + cᵢT^2 + dᵢT^3 + eᵢT^-2
Cp(T) = ∑Cpᵢxᵢ
```

## Model Construction Examples
```
# Using the default database
idealmodel = ShomateIdeal("water") #single input
idealmodel = ShomateIdeal(["water","ethanol"]) #multiple components

# Using user-provided parameters

# Passing files or folders
idealmodel = ShomateIdeal(["neon","hydrogen"]; userlocations = ["path/to/my/db","shomate.csv"])

# Passing parameters directly
idealmodel = ShomateIdeal(["water","butane"];
            userlocations = (a = [32.24, 9.487], 
                        b = [0.00192, 0.3313], 
                        c = [1.06e-5, -0.0001108],
                        d = [-3.6e-9, -2.822e-9])
                        ) #e is not used
```

"""
ShomateIdeal
export ShomateIdeal

default_locations(::Type{ShomateIdeal}) = ["ideal/ReidIdeal.csv","ideal/ShomateIdeal.csv"]
default_ignore_missing_singleparams(::Type{ShomateIdeal}) = ["e"]
function transform_params(::Type{ShomateIdeal},params,components)
    a,b,c,d = params["a"],params["b"],params["c"],params["d"]
    e = get(params,"e") do
        SingleParam("e",components)
    end
    params["coeffs"] = reid_coeffs(a,b,c,d,e,components)
    return params
end

function evalcoeff(::ShomateIdealModel,coeffs,T,lnT = log(T))
    A,B,C,D,E = coeffs
    poly = (A,B,C,D)
    return evalpoly(T,poly) + E/(T*T)
end

function eval∫coeff(::ShomateIdealModel,coeffs,T,lnT = log(T))
    A,B,C,D,E = coeffs
    ∫poly = (A,B/2,C/3,D/4)
    return evalpoly(T,∫poly)*T - E/T
end

function eval∫coeffT(::ShomateIdealModel,coeffs,T,lnT = log(T))
    A,B,C,D,E = coeffs
    ∫polyT = (B,C/2,D/3)
    return A*lnT + evalpoly(T,∫polyT)*T - E/(2*T*T)
end

