abstract type LKPSJTModel <: LKPModel end

struct LKPSJT{I} <: LKPSJTModel 
    components::Vector{String}
    params::LKPParam
    methane::SingleFluidResidualParam
    octane::SingleFluidResidualParam
    idealmodel::I
    references::Vector{String}
end

default_references(::Type{LKPSJT}) = ["10.1021/i260067a020","10.1007/s10765-024-03360-0"]
default_locations(::Type{LKPSJT}) = ["properties/critical.csv","properties/molarmass.csv"]

function transform_params(::Type{LKPSJT},params,components)
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
    for i in 1:length(Vc)
        if Vc.ismissingvalues[i]
            Vc[i] = (0.2905 - 0.085*ω[i])*Rgas()*Tc[i]/Pc[i]
        end
    end
    return params
end

"""
    LKPSJT <: LKPModel
    LKPSJT(components;
        version = 2,
        idealmodel=BasicIdeal,
        userlocations = String[],
        ideal_userlocations = String[],
        verbose=false,
        reference_state = nothing)
    
    enhancedLKP(components;
    idealmodel=BasicIdeal,
    verbose=false)

## Input parameters
- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Pressure `[Pa]`
- `Vc`: Single Parameter (`Float64`) (optional) - Critical Volume `[m³]`
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `acentricfactor`: Single Parameter (`Float64`) - Acentric Factor (no units)
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Parameter (no units)

## Input models
- `idealmodel`: Ideal Model

## Description
Lee-Kesler-Plöker-equation of state, Sabozin-Jäger-Thol enhancement. Corresponding states using interpolation between a simple, spherical fluid (methane, `∅`)  and a reference fluid (n-octane, `ref`):

Instead of using the original BWR formulation, multiparameter equations of state are used.
```
αᵣ = (1 - ωᵣ)*αᵣ(simple,δ,τ) + ωᵣ*αᵣ(reference,δ,τ)
τ = Tr/T
δ = Vr/V
ωᵣ = (ω̄ - ω(∅))/(ω(ref) - ω(∅))
ω̄ = ∑xᵢωᵢ
Tr = ∑xᵢ*xⱼ*Tcᵢⱼ*Vcᵢⱼ^η * (1-kᵢⱼ)
Vr = ∑xᵢ*xⱼ*Tcᵢⱼ*Vcᵢⱼ
Tcᵢⱼ = √Tcᵢ*Tcⱼ
Vcᵢⱼ = 0.125*(∛Vcᵢ + ∛Vcⱼ)^3
η = 0.25
```

If `Vcᵢ` is not provided, it can be estimated from the critical temperature and pressure:

```
Zcᵢ = Zᴵ + Zᴵᴵ ω
```

There are two versions of this EoS, that can be selected via the `version` argument.

- First version (`version = 1`): the simple fluid is the methane multiparameter EoS, and the reference fluid is octane. 
  For estimating the critical volumes, `Zc = 0.2905 - 0.085ω` is used. 
- Second version (`version = 1`): the simple fluid is replaced by the propane multiparameter EoS. 
  For estimating the critical volumes, Zc is an interpolation of the critical compressibilities from the simple and reference fluids: `Zc(ω0) = Z0, Zc(ωref) = Zref`
- `version = 2` : interpolation between the simple and the reference fluid.

## Model Construction Examples
```julia
# Using the default database
model = LKPSJT("water") #single input
model = LKPSJT(["water","ethanol"],version = 1) #multiple components, using original version
model = LKPSJT(["water","ethanol"], idealmodel = ReidIdeal) #modifying ideal model
model = enhancedLKP(["water","ethanol"], idealmodel = ReidIdeal) #enhancedLKP is just an alias to LKPSJT

# Passing a prebuilt model

my_idealmodel = MonomerIdeal(["neon","hydrogen"];userlocations = (;Mw = [20.17, 2.]))
model =  LKPSJT(["neon","hydrogen"],idealmodel = my_idealmodel)

# User-provided parameters, passing files or folders
model = LKPSJT(["neon","hydrogen"]; userlocations = ["path/to/my/db","lkp/my_k_values.csv"])

# User-provided parameters, passing parameters directly

model = LKPSJT(["neon","hydrogen"];
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
LKPSJT
const enhancedLKP = LKPSJT
@doc (@doc LKPSJT) enhancedLKP

function LKPSJT(components;
    idealmodel = BasicIdeal,
    userlocations = String[], 
    ideal_userlocations = String[],
    reference_state = nothing,
    verbose = false,
    version = 2,)

    reference_model = lkpsjt(version,verbose = verbose)
    return LKP(components; idealmodel, reference_model, userlocations, ideal_userlocations, reference_state, verbose)
end

export LKPSJT, enhancedLKP

struct LKPMultiparameter{S}  <: LKPReferenceModel
    simple::S
    reference::S
    ZI::Float64
    ZII::Float64
end

model_showname(::LKPMultiparameter) = "LKPMultiparameter"

function show_info_lkp(io,model::LKPMultiparameter)
    print(io,"LKP-SJT (")
    print(io,"simple fluid = ")
    print(io,only(component_list(model.simple)))
    print(io,", reference fluid = ")
    print(io,only(component_list(model.reference)))
    print(io,")")
end

function lkpsjt(version;verbose = false)
    if version == 1
        simple = SingleFluid("methane",verbose = verbose)
    else
        simple = PropaneRef()
    end
    
    reference = SingleFluid("octane",verbose = verbose)
    version == 1 && return LKPMultiparameter(simple,reference,0.2905,- 0.085)
    
    p0,pref = simple.properties,reference.properties
    ω0   = simple.properties.acentricfactor
    ωref = reference.properties.acentricfactor
    Zc0   = p0.Pc / (p0.Rgas * p0.Tc * p0.rhoc)
    Zcref = pref.Pc / (pref.Rgas * pref.Tc * pref.rhoc)
    ZI  = Zc0 - ((Zc0 - Zcref) / (ω0 - ωref)) * ω0
    ZII = (Zc0 - Zcref) / (ω0 - ωref)
    return LKPMultiparameter(simple,reference,ZI,ZII)
end

function lkp_Zc(model::LKPMultiparameter,ω)
    return model.ZI + model.ZII*ω
end

function a_res_lkp(model::LKPMultiparameter,V,T,data)
    T̄,V̄,ω̄ = data
    δ = V̄/V
    τ = T̄/T
    αr_0    = reduced_a_res(model.simple,δ,τ)
    αr_ref  = reduced_a_res(model.reference,δ,τ)
    ω0 = model.simple.properties.acentricfactor
    ωref = model.reference.properties.acentricfactor
    ωᵣ = (ω̄ - ω0)/(ωref - ω0)
    return (1 - ωᵣ)*αr_0 + ωᵣ*αr_ref
end