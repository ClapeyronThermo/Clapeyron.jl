struct NasaIdealParam <: EoSParam
    coeffs::SingleParam{NTuple{19,Float64}}
    reference_state::ReferenceState
    Mw::SingleParam{Float64}
end

abstract type NasaIdealModel <: IdealModel end
@newmodelsimple NasaIdeal NasaIdealModel NasaIdealParam

"""
    NasaIdeal <: IdealModel

    NasaIdeal(components;
    userlocations = String[],
    reference_state = nothing,
    verbose = false)

## CSV columns (per species)

- `Tmid`: Temperature splitting the low/high fits (typically 1000 K)
- `a1l,...,a7l`: NASA `a1..a7` for the **low-T** range
- `b1l,b2l`    : NASA `a8,a9` (integration constants) for the **low-T** range
- `a1h,...,a7h`: NASA `a1..a7` for the **high-T** range
- `b1h,b2h`    : NASA `a8,a9` (integration constants) for the **high-T** range

All `a*`/`b*` are **dimensionless** (NASA form). Cp is expressed as:

```
Cp/R = a1*T^-2 + a2*T^-1 + a3 + a4*T + a5*T^2 + a6*T^3 + a7*T^4
```

and the integrated forms are:

```
H/(R*T) = -a1/T^2 + a2*log(T)/T + a3 + a4*T/2 + a5*T^2/3 + a6*T^3/4 + a7*T^4/5 + b1/T
S/R     = -a1/(2*T^2) - a2/T + a3*log(T) + a4*T + a5*T^2/2 + a6*T^3/3 + a7*T^4/4 + b2
```

These are applied piecewise using the low/high sets with the split at `Tmid` for each component.
"""
NasaIdeal

export NasaIdeal


default_locations(::Type{NasaIdeal}) = ["ideal/NasaIdeal.csv","properties/molarmass.csv"]
default_ignore_missing_singleparams(::Type{NasaIdeal}) = ["Mw"]

# pack low/high + Tmid into a single 19-tuple per component:
# (a1l..a7l, b1l, b2l, a1h..a7h, b1h, b2h, Tmid)
function nasa_coeffs(a_low,b_low,a_high,b_high,Tmid,comps)
    _coeffs = fill(ntuple(_->0.0,19), length(comps))
    coeffs = SingleParam("NASA Coefficients", comps, _coeffs)
    return nasa_coeffs!(coeffs,a_low,b_low,a_high,b_high,Tmid)
end

function nasa_coeffs!(coeffs,a_low,b_low,a_high,b_high,Tmid)
    for i in 1:length(coeffs)
        c = (a_low[1][i],a_low[2][i],a_low[3][i],a_low[4][i],a_low[5][i],a_low[6][i],a_low[7][i],
             b_low[1][i],b_low[2][i],
             a_high[1][i],a_high[2][i],a_high[3][i],a_high[4][i],a_high[5][i],a_high[6][i],a_high[7][i],
             b_high[1][i],b_high[2][i],
             Tmid[i])
        coeffs[i] = c
    end
    return coeffs
end

function transform_params(::Type{NasaIdeal}, params, components)
    # Collect column SingleParams for a's and b's
    a_low  = ntuple(k->params["a$(k)l"], 7)
    b_low  = (params["b1l"], params["b2l"])
    a_high = ntuple(k->params["a$(k)h"], 7)
    b_high = (params["b1h"], params["b2h"])
    Tmid   = params["Tmid"]
    params["coeffs"] = nasa9_coeffs(a_low,b_low,a_high,b_high,Tmid,components)
    return params
end

@inline function _split19(c::NTuple{19,Float64})
    a1l,a2l,a3l,a4l,a5l,a6l,a7l, b1l,b2l,
    a1h,a2h,a3h,a4h,a5h,a6h,a7h, b1h,b2h, Tmid = c
    return (a1l,a2l,a3l,a4l,a5l,a6l,a7l,b1l,b2l), (a1h,a2h,a3h,a4h,a5h,a6h,a7h,b1h,b2h), Tmid
end

@inline function _pick(c::NTuple{19,Float64}, T)
    low, high, Tmid = _split19(c)
    return T <= Tmid ? low : high
end

# Cp(T) in J/mol/K
function evalcoeff(::NasaIdealModel, c::NTuple{19,Float64}, T, lnT = log(T))
    a1,a2,a3,a4,a5,a6,a7, b1,b2 = _pick(c,T)
    Tinv = inv(T)
    Tinv2 = Tinv*Tinv
    return R̄*(a1*Tinv2 + a2*Tinv + a3 + a4*T + a5*T*T + a6*T^3 + a7*T^4)
end

# H = ∫ Cp dT  (J/mol)   [indefinite up to a constant; differences cancel constants]
function eval∫coeff(::NasaIdealModel, c::NTuple{19,Float64}, T, lnT = log(T))
    a1,a2,a3,a4,a5,a6,a7, b1,b2 = _pick(c,T)
    T2 = T*T
    T3 = T2*T
    T4 = T3*T
    T5 = T4*T
    return R̄*(-a1/T + a2*lnT + a3*T + a4*T2/2 + a5*T3/3 + a6*T4/4 + a7*T5/5 + b1)
end

# S = ∫ Cp/T dT  (J/mol/K)
function eval∫coeffT(::NasaIdealModel, c::NTuple{19,Float64}, T, lnT = log(T))
    a1,a2,a3,a4,a5,a6,a7, b1,b2 = _pick(c,T)
    Tinv = inv(T)
    Tinv2 = Tinv*Tinv
    T2 = T*T
    T3 = T2*T
    T4 = T3*T
    return R̄*(-a1*(0.5)*Tinv2 - a2*Tinv + a3*lnT + a4*T + a5*T2/2 + a6*T3/3 + a7*T4/4 + b2)
end
