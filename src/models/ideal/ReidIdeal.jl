struct ReidIdealParam <: EoSParam
    coeffs::SingleParam{NTuple{4,Float64}}
end

function ReidIdealParam(a::SingleParam,b::SingleParam,c::SingleParam,d::SingleParam)
    comps = a.components
    a = a.values
    b = b.values
    c = c.values
    d = d.values
    return ReidIdealParam(a,b,c,d,comps)
end

function ReidIdealParam(a,b,c,d,comps)
    n = length(a)
    coeffs = [(a[i],b[i],c[i],d[i]) for i in 1:n]
    ReidIdealParam(SingleParam("Reid Coefficients",comps,coeffs))
end

abstract type ReidIdealModel <: IdealModel end
@newmodelsimple ReidIdeal ReidIdealModel ReidIdealParam

"""
    ReidIdeal <: IdealModel
    ReidIdeal(components; 
    userlocations::Array{String,1}=String[], 
    verbose=false)

## Input parameters

- `a`: Single Parameter (`Float64`)
- `b`: Single Parameter (`Float64`)
- `c`: Single Parameter (`Float64`)
- `d`: Single Parameter (`Float64`)

## Model parameters

- `coeffs`: Single Parameter (`NTuple{4,Float64}`)

## Description

Reid Ideal Model. Helmholtz energy obtained via integration of specific heat capacity:

```
Cpᵢ(T) = aᵢ  + bᵢT + cᵢT^2 + dᵢT^3
Cp(T) = ∑Cpᵢxᵢ
```
"""
ReidIdeal

export ReidIdeal
function ReidIdeal(components::Array{String,1}; userlocations::Array{String,1}=String[], verbose=false)
    params = getparams(components, ["ideal/ReidIdeal.csv"]; userlocations=userlocations, verbose=verbose)
    a = params["a"]
    b = params["b"]
    c = params["c"]
    d = params["d"]
    packagedparams = ReidIdealParam(a, b, c, d)
    references = String[] #  Fill this up.
    return ReidIdeal(packagedparams; references=references)
end

function a_ideal(model::ReidIdealModel, V, T, z)
    #x = z/sum(z)
    polycoeff = model.params.coeffs.values
    #return sum(x[i]*(log(z[i]/V) + 1/(R̄*T)*(sum(polycoeff[k][i]/k*(T^k-298^k) for k in 1:4)) -
    #    1/R̄*((polycoeff[k][1]-R̄)*log(T/298)+sum(polycoeff[k][i]/(k-1)*(T^(k-1)-298^(k-1)) for k in 2:4))) for i in @comps)

    res = zero(V+T+first(z))
    Σz = sum(z)
    @inbounds for i in @comps
        c0,c1,c2,c3 = polycoeff[i]
        pol1 = (c0,c1*0.5,c2/3,c3*0.25)
        pol2 = (c1,c2*0.5,c3/3)
        lnV = log(z[i]/V)
        lnT = (R̄ - c0)*(log(T/298))/R̄
        k1 = evalpoly(T,pol1)/R̄ - 298*evalpoly(298,pol1)/(R̄*T)
        k2 = T*evalpoly(T,pol2) - 298*evalpoly(298,pol2)
        res += z[i]*(lnV+lnT+k1-k2/R̄)
        #res +=x[i]*(log(z[i]/V) + 1/(R̄*T)*(sum(polycoeff[k]/k*(T^k-298^k) for k in 1:4)) -
        #1/R̄*((polycoeff[1]-R̄)*log(T/298)+sum(polycoeff[k]/(k-1)*(T^(k-1)-298^(k-1)) for k in 2:4)))
    end
    return res/Σz
end

function VT_isobaric_heat_capacity(model::ReidIdealModel,V,T,z=SA[1.])
    coeff = model.params.coeffs.values
    res = zero(T+first(z))
    Σz = sum(z)
    for i in @comps
        pol = coeff[i]
        res +=z[i]*evalpoly(T,pol)
    end
    return res/Σz
end