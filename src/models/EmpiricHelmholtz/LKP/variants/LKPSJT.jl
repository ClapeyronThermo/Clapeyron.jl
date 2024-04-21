abstract type LKPSJTModel <: LKPModel end

struct SKPSJT{I} <: LKPSJTModel 
    components::Vector{String}
    params::LKPParam
    methane::SingleFluidResidualParam
    octane::SingleFluidResidualParam
    idealmodel::I
    references::Vector{String}
end

default_references(::Type{SKPSJT}) = ["10.1021/i260067a020","10.1007/s10765-024-03360-0"]
default_locations(::Type{SKPSJT}) = ["properties/critical.csv","properties/molarmass.csv"]

function transform_params(::Type{SKPSJT},params,components)
    k = get(params,"k",nothing)
    if k === nothing
        nc = length(components)
        params["k"] = PairParam("k",components)
    end
    Vc = get(params,"Vc",nothing)
    if Vc === nothing
        params["Vc"] = SingleParam("Vc",components)
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
    SKPSJT <: LKPModel
    SKPSJT(components;
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
Lee-Kesler-Plöker-equation of state, Sabozin-Jäger-Thol improvement. Corresponding states using interpolation between a simple, spherical fluid (methane, `∅`)  and a reference fluid (n-octane, `ref`):

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
model = SKPSJT("water") #single input
model = SKPSJT(["water","ethanol"]) #multiple components
model = SKPSJT(["water","ethanol"], idealmodel = ReidIdeal) #modifying ideal model

# Passing a prebuilt model

my_idealmodel = MonomerIdeal(["neon","hydrogen"];userlocations = (;Mw = [20.17, 2.]))
model =  SKPSJT(["neon","hydrogen"],idealmodel = my_idealmodel)

# User-provided parameters, passing files or folders
model = SKPSJT(["neon","hydrogen"]; userlocations = ["path/to/my/db","lkp/my_k_values.csv"])

# User-provided parameters, passing parameters directly

model = SKPSJT(["neon","hydrogen"];
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
SKPSJT

function SKPSJT(components;
    idealmodel = BasicIdeal,
    userlocations = String[], 
    pure_userlocations = String[],
    reference_state = nothing,
    verbose = false)

    formatted_components = format_components(components)
    params = getparams(formatted_components, default_locations(SKPSJT); userlocations = userlocations, verbose = verbose)
    params = transform_params(SKPSJT,params,formatted_components)
    packagedparams = build_eosparam(LKPParam,params)
    references = default_references(SKPSJT)
    methane = SingleFluid("methane",verbose = verbose).residual
    octane = SingleFluid("octane",verbose = verbose).residual
    init_idealmodel = init_model(idealmodel,formatted_components,references,verbose,reference_state)
    model = SKPSJT(formatted_components,packagedparams,methane,octane,init_idealmodel,references)
    return model
end

export SKPSJT

lkp_params_simple(model::LKPSJTModel) = (false, 0.01142)
lkp_params_reference(model::LKPSJTModel) = (true, 0.3978)

function reduced_a_res_lkp(model::LKPSJTModel,δ,τ,δr,params)
    is_ref = first(params_reference)
    if is_ref
        #use the reference equations of state.
        return reduced_a_res(model.octane,δ,τ)
    else
        return reduced_a_res(model.methane,δ,τ)
    end
end