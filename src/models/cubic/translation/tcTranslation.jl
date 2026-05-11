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

- `ZRA`: Single Parameter (`Float64`) - Rackett compressibility factor 

## Model Parameters

- `ZRA`: Single Parameter (`Float64`) - Rackett compressibility factor
- `v_shift`: Single Parameter (`Float64`) - Volume shift `[m³·mol⁻¹]`


## Description

Translation model used in the tc-PR and tc-RK equation of state:
```
V = V₀ + mixing_rule(cᵢ)
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

1. Péneloux A, Rauzy E, Fréze R. (1982) A consistent correction for Redlich‐Kwong‐Soave volumes. Fluid Phase Equilibria 1, 8(1), 7–23. [doi:10.1016/0378-3812(82)80002-2](https://doi.org/10.1016/0378-3812(82)80002-2)
2. Ahlers, J., & Gmehling, J. (2001). Development of an universal group contribution equation of state. Fluid Phase Equilibria, 191(1–2), 177–188. [doi:10.1016/s0378-3812(01)00626-4](https://doi.org/10.1016/s0378-3812(01)00626-4) 

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
                c[i] = RTp*(0.1398 - 0.5294*ZRA[i])
            elseif !w.ismissingvalues[i]
                c[i] = RTp*evalpoly(w[i],(-0.014471,0.067498,-0.084852,0.067298,-0.017366))
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
                c[i] = RTp*(0.1487 - 0.5052*ZRA[i])
            
            elseif !w.ismissingvalues[i]
                c[i] = RTp*(0.0227 + 0.0093*w[i])
            else
                throw(error("tcTranslation: cannot estimate v_shift: missing acentricfactor or ZRA parameter"))
            end
        end
    end
    return c
end

export tcTranslation
