abstract type LKPmodModel <: LKPModel end

@newmodel LKPmod LKPmodModel LKPParam false

"""
    LKPmod <: LKPModel
    LKPmod(components;
        idealmodel=BasicIdeal,
        verbose=false)

## Input parameters
- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Pressure `[Pa]`
- `Vc`: Single Parameter (`Float64`) (optional) - Critical Volume `[m^3]`
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `acentricfactor`: Single Parameter (`Float64`) - Acentric Factor (no units)
- `k`: Pair Parameter (`Float64`) (optional) - binary interaction parameter (no units)

## Input models
- `idealmodel`: Ideal Model

## Description
Lee-Kesler-Plöker equation of state. corresponding states using interpolation between a simple, spherical fluid (methane, `∅`)  and a reference fluid (n-octane, `ref`):

Coefficients of the original EoS were modified in [2] guarantee a correct behaviour of the EoS in the liquid phase.
```
αᵣ = (1 - ωᵣ)*αᵣ(δr,τ,params(∅)) + ωᵣ*αᵣ(δr,τ,params(ref))
τ = Tr/T
δr = Vr/V/Zr
Zr = Pr*Vr/(R*Tr)
Pr = (0.2905 - 0.085*ω̄)*R*Tr/Vr
ωᵣ = (ω̄ - ω(∅))/(ω(ref) - ω(∅))
ω̄ = ∑xᵢωᵢ
Tr = ∑xᵢ*xⱼ*Tcᵢⱼ*Vcᵢⱼ^η * (1-kᵢⱼ)
Vr = ∑xᵢ*xⱼ*Tcᵢⱼ*Vcᵢⱼ
Tcᵢⱼ = √Tcᵢ*Tcⱼ
Vcᵢⱼ = 0.125*(∛Vcᵢ + ∛Vcⱼ)^3
η = 0.25
```

## Model Construction Examples
```julia
# Using the default database
model = LKPmod("water") #single input
model = LKPmod(["water","ethanol"]) #multiple components
model = LKPmod(["water","ethanol"], idealmodel = ReidIdeal) #modifying ideal model

# Passing a prebuilt model

my_idealmodel = MonomerIdeal(["neon","hydrogen"];userlocations = (;Mw = [20.17, 2.]))
model =  LKPmod(["neon","hydrogen"],idealmodel = my_idealmodel)

# User-provided parameters, passing files or folders
model = LKPmod(["neon","hydrogen"]; userlocations = ["path/to/my/db","lkp/my_k_values.csv"])

# User-provided parameters, passing parameters directly

model = LKPmod(["neon","hydrogen"];
        userlocations = (;Tc = [44.492,33.19],
                        Pc = [2679000, 1296400],
                        Vc = [4.25e-5, 6.43e-5],
                        Mw = [20.17, 2.],
                        acentricfactor = [-0.03,-0.21]
                        k = [0. 0.18; 0.18 0.]) #k,l can be ommited in single-component models.
                    )
```

## References
1. Plöcker, U., Knapp, H., & Prausnitz, J. (1978). Calculation of high-pressure vapor-liquid equilibria from a corresponding-states correlation with emphasis on asymmetric mixtures. Industrial & Engineering Chemistry Process Design and Development, 17(3), 324–332. [doi:10.1021/i260067a020](https://doi.org/10.1021/i260067a020)
2. Sabozin, F., Jäger, A., & Thol, M. (2024). Enhancement of the Lee–Kesler–Plöcker equation of state for calculating thermodynamic properties of long-chain alkanes. International Journal of Thermophysics, 45(5). [doi:10.1007/s10765-024-03360-0](https://doi.org/10.1007/s10765-024-03360-0)
"""
LKPmod

export LKPmod

default_references(::Type{LKPmod}) = ["10.1021/i260067a020","10.1007/s10765-024-03360-0"]
default_locations(::Type{LKPmod}) = ["properties/critical.csv","properties/molarmass.csv"]
default_ignore_missing_singleparams(::Type{LKPmod}) = ["Vc"]

function transform_params(::Type{LKPmod},params,components)
    k = get(params,"k",nothing)
    if k === nothing
        nc = length(components)
        params["k"] = PairParam("k",components)
    end
    _Vc = get(params,"Vc",nothing)
    if _Vc === nothing
        Vc = SingleParam("Vc",components)
        params["Vc"] = Vc
    else
        Vc = _Vc
    end
    Tc,Pc,ω = params["Tc"],params["Pc"],params["acentricfactor"]
    for i in 1:length(Vc.values)
        if Vc.ismissingvalues[i]
            Vc.values[i] = (0.2905 - 0.085*ω[i])*Rgas()*Tc[i]/Pc[i]
        end
    end
    return params
end

lkp_params_simple(model::LKPmodModel) = (0.1331199, 0.3392959, 0.0786113, 0.0498273, 0.0218093, 0.010958, 0.0050041, 0.0309082, 1.9876201e-5, 3.4930069e-5, 0.585946, 0.0677684, 0.01142)
lkp_params_reference(model::LKPmodModel) = (0.0243243, 0.0640205, 0.0899694, 0.2313499, 0.0647721, 0.0928313, 0.0154748, 0.04441, 2.0525725e-5, 3.5470136e-5, 1.4003447, 0.0286862, 0.3978)
