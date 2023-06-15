struct CPLNGEstIdealParam <: EoSParam
    coeffs::SingleParam{NTuple{4,Float64}}
    Mw::SingleParam{Float64}
end

@newmodelsimple CPLNGEstIdeal ReidIdealModel CPLNGEstIdealParam


"""
    CPLNGEstIdeal <: ReidIdealModel
    CPLNGEstIdeal(components; 
    userlocations::Array{String,1}=String[], 
    verbose=false)

## Input parameters

- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`

## Model parameters

- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `coeffs`: Single Parameter (`NTuple{4,Float64}`) - polynomial coefficients

## Description

Estimation of Reid polynomial, using `Tc`, `Pc` and `MW` as inputs

```
Cpᵢ(T) = aᵢ  + bᵢT + cᵢT^2 + dᵢT^3
Cp(T) = ∑Cpᵢxᵢ
a = -10.9602   * γ₀ + 25.9033
b = 2.1517e-1  * γ₀ - 6.8687e-2 
c = -1.3337e-4 * γ₀ + 8.6387e-5
d = 3.1474e-8  * γ₀ -2.8396e-8
```
"""
CPLNGEstIdeal

export CPLNGEstIdeal
function CPLNGEstIdeal(components::Array{String,1}; userlocations::Array{String,1}=String[], verbose=false)
    params = getparams(components, ["properties/molarmass.csv"]; userlocations=userlocations, verbose=verbose)
    Mw = params["Mw"]
    
    γ₀ = Mw ./ 28.9647
#=
a1 = -10.9602
a2 = 25.9033
b1 = 2.1517e-1
b2 = -6.8687e-2
c1 = -1.3337e-4
c2 = 8.6387e-5
d1 = 3.1474e-8
d2 = -2.8396e-8 =#

a = -10.9602   .* γ₀ .+ 25.9033
b = 2.1517e-1  .* γ₀ .- 6.8687e-2 
c = -1.3337e-4 .* γ₀ .+ 8.6387e-5
d = 3.1474e-8  .* γ₀ .- 2.8396e-8
ReidIdealParam(a,b,c,d,components)
    packagedparams = ReidIdealParam(a, b, c, d)
    references = String["http://dx.doi.org/10.1016/j.jngse.2014.04.011"] #  Fill this up.
    
    return ReidIdeal(packagedparams; references=references)
end

#TODO,add a dependency of a,b,c,d parameters
recombine_impl!(model::ReidIdealModel) = model