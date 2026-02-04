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

Instead of using the original BWR formulation, the reference equations of state of methane and octane are used.
```
αᵣ = (1 - ωᵣ)*αᵣ(δ,τ,params(∅)) + ωᵣ*αᵣ(δ,τ,params(ref))
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

## Model Construction Examples
```julia
# Using the default database
model = LKPSJT("water") #single input
model = LKPSJT(["water","ethanol"]) #multiple components
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
    verbose = false)

    formatted_components = format_components(components)
    params = getparams(formatted_components, default_locations(LKPSJT); userlocations = userlocations, verbose = verbose,ignore_missing_singleparams = ["Vc"])
    params = transform_params(LKPSJT,params,formatted_components)
    packagedparams = build_eosparam(LKPParam,params)
    references = default_references(LKPSJT)
    methane = SingleFluid("methane",verbose = verbose).residual
    octane = SingleFluid("octane",verbose = verbose).residual
    init_idealmodel = init_model(idealmodel,formatted_components,references,verbose)
    model = LKPSJT(formatted_components,packagedparams,methane,octane,init_idealmodel,references)
    set_reference_state!(model,reference_state;verbose)
    return model
end

export LKPSJT, enhancedLKP

lkp_params_simple(model::LKPSJTModel) = (false, 0.01142)
lkp_params_reference(model::LKPSJTModel) = (true, 0.3978)

function reduced_a_res_lkp(model::LKPSJTModel,δ,τ,δr,params)
    is_ref = first(params)
    if is_ref
        #use the reference equations of state.
        return reduced_a_res(model.octane,δ,τ)
    else
        return reduced_a_res(model.methane,δ,τ)
    end
end
