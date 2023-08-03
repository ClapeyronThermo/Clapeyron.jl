struct ReidIdealParam <: EoSParam
    coeffs::SingleParam{NTuple{4,Float64}}
end



function reid_coeffs(a::SingleParam,b::SingleParam,c::SingleParam,d::SingleParam)
    comps = a.components
    a = a.values
    b = b.values
    c = c.values
    d = d.values
    return reid_coeffs(a,b,c,d,comps)
end

function reid_coeffs(a,b,c,d,comps)
    n = length(a)
    coeffs = [(a[i],b[i],c[i],d[i]) for i in 1:n]
    SingleParam("Reid Coefficients",comps,coeffs)
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
default_locations(::Type{ReidIdeal}) = ["ideal/ReidIdeal.csv"]
function transform_params(::Type{ReidIdeal},params)
    a,b,c,d = params["a"],params["b"],params["c"],params["d"]
    params["coeffs"] = reid_coeffs(a,b,c,d)
    return params
end
recombine_impl!(model::ReidIdealModel) = model


function a_ideal(model::ReidIdealModel, V, T, z)
    #x = z/sum(z)
    polycoeff = model.params.coeffs.values
    #return sum(x[i]*(log(z[i]/V) + 1/(R̄*T)*(sum(polycoeff[k][i]/k*(T^k-298^k) for k in 1:4)) -
    #    1/R̄*((polycoeff[k][1]-R̄)*log(T/298)+sum(polycoeff[k][i]/(k-1)*(T^(k-1)-298^(k-1)) for k in 2:4))) for i in @comps)
    V⁻¹ = 1/V
    res = zero(V+T+first(z))
    Σz = sum(z)
    @inbounds for i in @comps
        ci = polycoeff[i]
        n = length(ci)
        c0 = first(ci)
        cii = last(ci,n-1)
        div1 = NTuple{n,Int}(1:n)
        div2 = NTuple{n-1,Int}(1:n)
        R̄⁻¹= 1/R̄
        pol1 = ci ./ div1
        pol2 = cii ./ div2
        lnT = (1 - c0*R̄⁻¹)*(log(T/298))
        H = (evalpoly(T,pol1) - 298*evalpoly(298,pol1)/T)*R̄⁻¹
        TS = (T*evalpoly(T,pol2) - 298*evalpoly(298,pol2))*R̄⁻¹
        res += z[i]*(lnT+H-TS)
        res += xlogx(z[i],V⁻¹)
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
