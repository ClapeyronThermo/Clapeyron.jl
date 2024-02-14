function XiangDeitersConsts(ω,θ)
    d = [1, 1, 1, 2, 3, 7, 1, 1, 2, 5, 1, 1, 4, 2]
    t = [0.25, 1.25, 1.5, 1.375, 0.25, 0.875, 0, 2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5]
    l = [1, 1, 1, 1, 2, 2, 2, 3]
    a0 = [8.5740489E-01,  -3.2863233E+00, 1.6480939E+00,  -5.4524817E-02, 6.1623592E-02, 2.7389266E-04,  -6.0655087E-02,
    -3.1811852E-02, -1.1550422E-01, -1.8610466E-02, -1.8348671E-01, 5.5071325E-03, -1.2268039E-02, -5.0433436E-03]
    a1 = [5.6200117E-01, 3.2439544E+00, -4.9628768E+00, -2.2132851E-01, 9.3388356E-02, 2.4940171E-05,  -1.7519204E-01,
    8.9325660E-01, 2.9886613E+00, 1.0881387E-01,  -6.7166746E-01, 1.4477326E-01, -2.8716809E-01, -1.1478402E-01]
    a2 = [-8.1680511E+01, 4.6384732E+02,  -2.7970850E+02, 2.9317364E+01, -2.2324825E+01, -5.0932691E-02, -7.2836590E+00,
    -2.2063100E+02, -3.0435126E+02, 5.8514719E+00,  1.7995451E+02, -1.0178400E+02, 4.0848053E+01,  1.2411984E+01]
    n = a0 .+ a1 .* ω + a2 .* θ
    return SingleFluidResidualParam(n,t,d,l)
end

"""
    XiangDeiters::SingleFluid
    XiangDeiters(component;
        idealmodel = BasicIdeal,
        userlocations = String[],
        ideal_userlocations = String[],
        Rgas = nothing,
       verbose = false)

## Input parameters
- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Pressure `[Pa]`
- `Vc`: Single Parameter (`Float64`) - Critical Volume `[m3/mol]`
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `acentricfactor`: Single Parameter (`Float64`)
  
## Input models
- `idealmodel`: Ideal Model

## Description

Xiang-Deiters model. estimates a single component Empiric EoS Model from critical values and the acentric factor.
```
Zc = PcVc/RTc
θ = (Zc - 0.29)^2
aᵣ = a₀(δ,τ) + ω*a₁(δ,τ) + θ*a₂(δ,τ)
```

`Rgas` can be used to set the value of the gas constant that is used during property calculations.

## Model Construction Examples
```julia
# Using the default database
model = XiangDeiters("water") #single input
model = XiangDeiters(["water"]) #single input, as a vector
model = XiangDeiters(["water"], idealmodel = ReidIdeal) #modifying ideal model

# Passing a prebuilt model

my_idealmodel = MonomerIdeal(["ethane"])
model = XiangDeiters(["ethane"],idealmodel = my_idealmodel)

# User-provided parameters, passing files or folders
model = XiangDeiters(["hydrogen"]; userlocations = ["path/to/my/db","critical.csv"])

# User-provided parameters, passing parameters directly

model = XiangDeiters(["hydrogen"];
        userlocations = (;Tc = [44.492],
                        Pc = [2679000],
                        Vc = [4.25e-5],
                        Mw = [2.0],
                        acentricfactor = [-0.21])
                    )
```

## references
1. Xiang, H. W., & Deiters, U. K. (2008). A new generalized corresponding-states equation of state for the extension of the Lee–Kesler equation to fluids consisting of polar and larger nonpolar molecules. Chemical Engineering Science, 63(6), 1490–1496. [doi:10.1016/j.ces.2007.11.029](https://doi.org/10.1016/j.ces.2007.11.029)
"""
function XiangDeiters(component;
    idealmodel = BasicIdeal,
    userlocations = String[],
    ideal_userlocations = String[],
    Rgas = nothing,
    verbose = false)

    components = [get_only_comp(component)]
    params = getparams(components, ["properties/critical.csv", "properties/molarmass.csv"]; userlocations = userlocations, verbose = verbose)
    Pc = params["Pc"][1]
    Vc = params["Vc"][1]
    Mw = params["Mw"][1]
    Tc = params["Tc"][1]
    acentricfactor = params["acentricfactor"][1]

    Zc = Pc*Vc/(R̄*Tc)
    θ = (Zc - 0.29)^2
    rhoc = 1/Vc
    lb_volume = 1/(3.25*rhoc)
    residual = XiangDeitersConsts(acentricfactor,θ)
    if Rgas === nothing
        R = R̄
    else
        R = Rgas
    end
    properties = SingleFluidProperties(Mw,Tc,rhoc,lb_volume,Tc,Pc,rhoc,NaN,NaN,NaN,NaN,acentricfactor,R)
    ancillaries = propane_ancillary_cs(components,Tc,Pc,Vc)
    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    ideal_dict_data = Clapeyron.idealmodel_to_json_data(init_idealmodel; Tr = Tc, Vr = Vc)
    ideal = _parse_ideal(ideal_dict_data,verbose)
    SingleFluid(components,properties,ancillaries,ideal,residual,String["10.1016/j.ces.2007.11.029"])
end

export XiangDeiters