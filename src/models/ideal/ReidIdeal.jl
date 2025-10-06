abstract type PolynomialIdealModel <: IdealModel end

struct ReidIdealParam <: EoSParam
    a::SingleParam{Float64}
    b::SingleParam{Float64}
    c::SingleParam{Float64}
    d::SingleParam{Float64}
    e::SingleParam{Float64}
    coeffs::SingleParam{NTuple{5,Float64}}
    reference_state::ReferenceState
    Mw::SingleParam{Float64}
end

function reid_coeffs(a,b,c,d,e,comps)
    _coeffs = fill((0.0,0.0,0.0,0.0,0.0),length(comps))
    coeffs = SingleParam("Reid Coefficients",comps,_coeffs)
    return reid_coeffs!(coeffs,a,b,c,d,e)
end

reid_coeffs(a,b,c,d,comps) = reid_coeffs(a,b,c,d,FillArrays.Zeros(length(a)),comps)

function reid_coeffs!(coeffs,a,b,c,d,e)
    for i in 1:length(coeffs)
        coeffs[i] = (a[i],b[i],c[i],d[i],e[i])
    end
    return coeffs
end

reid_coeffs!(coeffs,a,b,c,d,e::Nothing) = reid_coeffs(coeffs,a,b,c,d,FillArrays.Zeros(length(coeffs)))


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
- `Mw`: Single Parameter (`Float64`) (Optional) - Molecular Weight `[g·mol⁻¹]`

## Model parameters

- `a`: Single Parameter (`Float64`) - polynomial coefficient
- `b`: Single Parameter (`Float64`) - polynomial coefficient
- `c`: Single Parameter (`Float64`) - polynomial coefficient
- `d`: Single Parameter (`Float64`) - polynomial coefficient
- `e`: Single Parameter (optional) (`Float64`)  - polynomial coefficient for 1/T^2
- `coeffs`: Single Parameter (`NTuple{5,Float64}`)
- `Mw`: Single Parameter (`Float64`) (Optional) - Molecular Weight `[g·mol⁻¹]`

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
                        d = [-3.6e-9, -2.822e-9],
                        Mw = [18.01, 58.12])
                        ) #e is not used
```

"""
ReidIdeal

export ReidIdeal
default_locations(::Type{ReidIdeal}) = ["ideal/ReidIdeal.csv","properties/molarmass.csv"]
default_ignore_missing_singleparams(::Type{ReidIdeal}) = ["e","Mw"]
function transform_params(::Type{ReidIdeal},params,components)
    a,b,c,d = params["a"],params["b"],params["c"],params["d"]
    e = get(params,"e") do
        SingleParam("e",components)
    end
    params["coeffs"] = reid_coeffs(a,b,c,d,e,components)
    return params
end

function recombine_impl!(model::ReidIdealModel)
    p = model.params
    reid_coeffs!(p.coeffs,p.a,p.b,p.c,p.d,p.e)
end

evalcoeff(::ReidIdealModel,coeffs,T,lnT = log(T)) = evalpoly(T,coeffs)

function eval∫coeff(::ReidIdealModel,coeffs,T,lnT = log(T))
    return Solvers.evalpolyint(T,coeffs)
end

function eval∫coeffT(::ReidIdealModel,coeffs,T,lnT = log(T))
    A = first(coeffs)
    coeffs1 = coeffs[2:end]
    return Solvers.evalpolyint(T,coeffs1) + A*lnT
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

#used for gibbs based models
function gibbs_cp_integral(model::PolynomialIdealModel,T,z,T0)
        #x = z/sum(z)
    polycoeff = model.params.coeffs.values
    #return sum(x[i]*(log(z[i]/V) + 1/(R̄*T)*(sum(polycoeff[k][i]/k*(T^k-298^k) for k in 1:4)) -
    #    1/R̄*((polycoeff[k][1]-R̄)*log(T/298)+sum(polycoeff[k][i]/(k-1)*(T^(k-1)-298^(k-1)) for k in 2:4))) for i in @comps)
    res = zero(T+first(z))
    Σz = sum(z)
    RT = R̄*T
    R̄⁻¹ = 1/R̄
    RT⁻¹ = 1/RT
    lnT0 = log(T0)
    lnT = log(T)
    @inbounds for i in @comps
        coeffs = polycoeff[i] 
        H = (eval∫coeff(model,coeffs,T,lnT) - eval∫coeff(model,coeffs,T0,lnT0))*RT⁻¹
        TS = (eval∫coeffT(model,coeffs,T,lnT) - eval∫coeffT(model,coeffs,T0,lnT0))*R̄⁻¹
        α₀ᵢ = H - TS + lnT - lnT0
        res += z[i]*α₀ᵢ
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