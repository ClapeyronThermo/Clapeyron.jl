struct CPLNGEstIdealParam <: EoSParam
    coeffs::SingleParam{NTuple{5,Float64}}
    Mw::SingleParam{Float64}
end

@newmodelsimple CPLNGEstIdeal ReidIdealModel CPLNGEstIdealParam

"""
    CPLNGEstIdeal <: ReidIdealModel

    CPLNGEstIdeal(components; 
    userlocations = String[],
    reference_state = nothing,
    verbose = false)

## Input parameters

- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`

## Model parameters

- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `coeffs`: Single Parameter (`NTuple{5,Float64}`) - polynomial coefficients

## Description

Estimation of Reid polynomial, using the molecular weight as input:

```
Cpᵢ(T) = aᵢ  + bᵢT + cᵢT^2 + dᵢT^3
Cp(T) = ∑Cpᵢxᵢ
a = -10.9602   * γ₀ + 25.9033
b = 2.1517e-1  * γ₀ - 6.8687e-2 
c = -1.3337e-4 * γ₀ + 8.6387e-5
d = 3.1474e-8  * γ₀ -2.8396e-8
γ₀ = Mw/Mw(air)
```

## Model Construction Examples
```
# Using the default database
idealmodel = CPLNGEstIdeal("water") #single input
idealmodel = CPLNGEstIdeal(["water","ethanol"]) #multiple components

# Using user-provided parameters

# Passing files or folders
idealmodel = CPLNGEstIdeal(["neon","hydrogen"]; userlocations = ["path/to/my/db","mw.csv"])

# Passing parameters directly
idealmodel = CPLNGEstIdeal(["neon","hydrogen"];userlocations = (;Mw = [20.17, 2.]))
```

## References
1. Kareem, L. A., Iwalewa, T. M., & Omeke, J. E. (2014). Isobaric specific heat capacity of natural gas as a function of specific gravity, pressure and temperature. Journal of Natural Gas Science and Engineering, 19, 74–83. [doi:10.1016/j.jngse.2014.04.011]("http://doi.org/10.1016/j.jngse.2014.04.011")
"""
CPLNGEstIdeal

default_locations(::Type{CPLNGEstIdeal}) = mw_data()
default_references(::Type{CPLNGEstIdeal}) = ["10.1016/j.jngse.2014.04.011"]
function transform_params(::Type{CPLNGEstIdeal},params,components)
    Mw = params["Mw"]   
    γ₀ = Mw ./ 28.9647
    a = -10.9602   .* γ₀ .+ 25.9033
    b = 2.1517e-1  .* γ₀ .- 6.8687e-2 
    c = -1.3337e-4 .* γ₀ .+ 8.6387e-5
    d = 3.1474e-8  .* γ₀ .- 2.8396e-8
    params["coeffs"] = reid_coeffs(a,b,c,d,components)
    return params
end

export CPLNGEstIdeal