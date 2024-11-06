abstract type PolynomialCpLiquidModel <: LiquidCpModel end

@newmodelsimple PolynomialCpLiquid PolynomialCpLiquidModel ReidIdealParam

"""
    PolynomialCpLiquid <: LiquidCpModel

    PolynomialCpLiquid(components; 
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

Incompressible fluid caloric properties. Helmholtz energy obtained via integration of heat capacity:

```
Cpᵢ(T) = aᵢ + bᵢT + cᵢT^2 + dᵢT^3 + eᵢT^4
Cp(T) = ∑Cpᵢxᵢ
```

## Model Construction Examples
```
# Using the default database
model = PolynomialCpLiquid("water") #single input
model = PolynomialCpLiquid(["water","ethanol"]) #multiple components

# Using user-provided parameters

# Passing files or folders
idealmodel = PolynomialCpLiquid(["neon","hydrogen"]; userlocations = ["path/to/my/db","cpl.csv"])

# Passing parameters directly
idealmodel = PolynomialCpLiquid(["water","butane"];
            userlocations = (a = [32.24, 9.487], 
                        b = [0.00192, 0.3313], 
                        c = [1.06e-5, -0.0001108],
                        d = [-3.6e-9, -2.822e-9])
                        ) #e is not used
```

"""
PolynomialCpLiquid

export PolynomialCpLiquid
#default_locations(::Type{PolynomialCpLiquid}) = ["ideal/ReidIdeal.csv"]
default_ignore_missing_singleparams(::Type{PolynomialCpLiquid}) = ["e"]
function transform_params(::Type{PolynomialCpLiquid},params,components)
    a,b,c,d = params["a"],params["b"],params["c"],params["d"]
    e = get(params,"e") do
        SingleParam("e",components)
    end
    params["coeffs"] = reid_coeffs(a,b,c,d,e,components)
    return params
end

recombine_impl!(model::PolynomialCpLiquid) = model

function a_res(model::PolynomialCpLiquidModel, V, T, z)
    #x = z/sum(z)
    polycoeff = caloric_coefficients(model)
    res = zero(Base.promote_eltype(model,T,z))
    Σz = sum(z)
    RT = R̄*T
    R̄⁻¹ = 1/R̄
    RT⁻¹ = 1/RT
    T0 = 298.
    lnT = log(T)
    @inbounds for i in @comps
        coeffs = polycoeff[i] 
        H = eval∫coeff(model,coeffs,T,lnT)*RT⁻¹
        TS = eval∫coeffT(model,coeffs,T,lnT)*R̄⁻¹
        α₀ᵢ = H - TS - lnT*R̄⁻¹
        res += z[i]*α₀ᵢ
    end
    return res/Σz
end
