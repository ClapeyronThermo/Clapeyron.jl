struct tcTranslationParam <: EoSParam
    acentricfactor::SingleParam{Float64}
    ZRA::SingleParam{Float64}
    v_shift::SingleParam{Float64}
end

@newmodelsimple tcTranslation ConstantTranslationModel tcTranslationParam

"""

    tcTranslation <: tcTranslationModel

    tcTranslation(components;
    userlocations = String[],
    verbose::Bool=false)

## Input Parameters

- `ZRA`: Single Parameter (`Float64`) - Rackett compressibility factor.

## Model Parameters

- `ZRA`: Single Parameter (`Float64`) - Rackett compressibility factor.
- `v_shift`: Single Parameter (`Float64`) - Volume shift `[m³·mol⁻¹]`.
- `acentricfactor`: Single Parameter (`Float64`) - Acentric factor.  

## Description

Translation model used in the tc-PR and tc-RK equation of state:

```
V = V₀ + mixing_rule(cᵢ)
```

For PR, if no Translation Volume parameters are provided can be estimated from acentric factor (ω) or Rackett compressibility factor (ZRA) (From the 2025 paper (default), both the 2018 and 2016 versions are also available):

From acentric factor (ω):

2025 version :

```
cᵢ(ωᵢ) = (R·T_c)/(P_c) * (-0.014471 + 0.067498*ωᵢ - 0.084852*ωᵢ² + 0.067298*ωᵢ³ - 0.017366*ωᵢ⁴)

```

2018 version :

```
cᵢ(ωᵢ) = (R·T_c)/(P_c) * (0.0096 + 0.0048*ωᵢ)

```

2016 version :

```
cᵢ(ωᵢ) = (R·T_c)/(P_c) * (-0.0065 + 0.0198*ωᵢ)

```

From Rackett compressibility  factor (ZRA):

2018 version :

```
cᵢ(ZRA) = (R·T_c)/(P_c) * (0.1975 - 0.7325*ZRA)

```

2016 version :

```
cᵢ(ZRA) = (R·T_c)/(P_c) * (0.1398 - 0.5294*ZRA)

```

For RK, if no Translation Volume parameters are provided can be estimated from acentricfactor (ω) or Rackett compressibility factor (ZRA) (From the 2018 paper (default), 2016 version is also available):

From acentric factor (ω):

2018 version :

```
cᵢ(ωᵢ) = (R·T_c)/(P_c) * (0.0227 + 0.0093*ωᵢ)

```

2016 version :

```
cᵢ(ωᵢ) = (R·T_c)/(P_c) * (0.0096 + 0.0172*ωᵢ)

```

From Rackett compressibility  factor (ZRA):

2018 version :

```
cᵢ(ZRA) = (R·T_c)/(P_c) * (0.2150 - 0.7314*ZRA)

```

2016 version :

```
cᵢ(ZRA) = (R·T_c)/(P_c) * (0.1487 - 0.5052*ZRA)

```

## Model Construction Examples
```
# Using the default database
translation = tcTranslation("water") #single input
translation = tcTranslation(["water","ethanol"]) #multiple components

# Using user-provided parameters

# Passing files or folders
translation = tcTranslation(["neon","hydrogen"]; userlocations = ["path/to/my/db","critical/Vc.csv"])

# Passing parameters directly
translation = tcTranslation(["neon","hydrogen"];userlocations = (;Vc = [4.25e-5, 6.43e-5]))
```

## References

1. Le Guennec, Y., Privat, R., & Jaubert, J.-N. (2016). Development of the translated-consistent tc-PR and tc-RK cubic equations of state for a safe and accurate prediction of volumetric, energetic and saturation properties of pure compounds in the sub- and super-critical domains. Fluid Phase Equilibria, 429, 301–312. [doi:10.1016/j.fluid.2016.09.003](http://dx.doi.org/10.1016/j.fluid.2016.09.003)
2. Pina-Martinez, A., Le Guennec, Y., Privat, R., Jaubert, J.-N., & Mathias, P. M. (2018). Analysis of the combinations of property data that are suitable for a safe estimation of consistent twu α-function parameters: Updated parameter values for the translated-consistent tc-PR and tc-RK cubic equations of state. Journal of Chemical and Engineering Data, 63(10), 3980–3988. [doi:10.1021/acs.jced.8b00640](http://dx.doi.org/10.1021/acs.jced.8b00640)
3. Paes, F., de Souza Batalha, G., Citrangolo Destro, F., Fournet, R., Privat, R., Jaubert, J-N., Sirjean, B. (2025). Integrating Solvent Effects into the Prediction of Kinetic Constants Using a COSMO-Based Equation of State. Journal of Chemical Theory and Computation, 21, 3625-3648. [doi:10.1021/acs.jctc.5c00133](https://doi.org/10.1021/acs.jctc.5c00133)
4. Magoulas, K., & Tassios, D. (1990). Thermophysical properties of n-Alkanes from C1 to C20 and their prediction for higher ones. Fluid Phase Equilibria, 56, 119–140. [doi:10.1016/0378-3812(90)85098-u](https://doi.org/10.1016/0378-3812(90)85098-u)

"""
tcTranslation

default_locations(::Type{tcTranslation}) = critical_data()
default_references(::Type{tcTranslation}) = ["10.1016/0378-3812(82)80002-2"]

function transform_params(::Type{tcTranslation},params,components)
    nc = length(components)
    #just initialization
    ZRA = get!(params,"ZRA") do
        SingleParam("ZRA",components,zeros(nc),fill(true,nc))
    end

    v_shift = get!(params,"v_shift") do
        SingleParam("v_shift",components,zeros(nc),fill(true,nc))
    end

    acentricfactor = get!(params,"acentricfactor") do
        SingleParam("acentricfactor",components,zeros(nc),fill(true,nc))
    end
    return params
end

function recombine_translation!(model::CubicModel,translation_model::tcTranslation)
    c = translation_model.params.v_shift
    translation!(c,model,translation_model)
    return translation_model
end

function translation!(c,model::PRModel,translation_model::tcTranslation)
    Tc = model.params.Tc.values
    Pc = model.params.Pc.values
    ZRA = translation_model.params.ZRA
    w = translation_model.params.acentricfactor
    for i in 1:length(model)
        if c.ismissingvalues[i]       
            R = Rgas(model)
            RTp = (R*Tc[i]/Pc[i])
            if !ZRA.ismissingvalues[i]
                #c[i] = RTp*(0.1398 - 0.5294*ZRA[i]) #2016 Version
                c[i] = RTp*(0.1975 - 0.7325*ZRA[i]) #2018 Version
            elseif !w.ismissingvalues[i]
                #c[i] = RTp*(-0.0065 + 0.0198*w[i]) #2016 Version
                #c[i] = RTp*(0.0096 + 0.0048 *w) #2018 Version
                c[i] = RTp*evalpoly(w[i],(-0.014471,0.067498,-0.084852,0.067298,-0.017366)) #2025 Version
            else
                throw(error("tcTranslation: cannot estimate v_shift: missing acentricfactor or ZRA parameter"))
            end
        end
    end
    return c
end

function translation!(c,model::RKModel,translation_model::tcTranslation)
    Tc = model.params.Tc.values
    Pc = model.params.Pc.values
    ZRA = translation_model.params.ZRA
    w = translation_model.params.acentricfactor
    for i in 1:length(model)
        if c.ismissingvalues[i]       
            R = Rgas(model)
            RTp = (R*Tc[i]/Pc[i])
            if !ZRA.ismissingvalues[i]
                #c[i] = RTp*(0.1487 - 0.5052*ZRA[i]) #2016 Version
                c[i] = RTp*(0.2150 - 0.7314*ZRA[i]) #2018 Version
            elseif !w.ismissingvalues[i]
                #c[i] = RTp*(0.0096 + 0.0172*w[i]) #2016 Version
                c[i] = RTp*(0.0227 + 0.0093*w[i]) #2018 Version
            else
                throw(error("tcTranslation: cannot estimate v_shift: missing acentricfactor or ZRA parameter"))
            end
        end
    end
    return c
end

export tcTranslation
