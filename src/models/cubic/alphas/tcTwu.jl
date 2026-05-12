struct tcTwuAlphaParam <: EoSParam
    M::SingleParam{Float64}
    N::SingleParam{Float64}
    L::SingleParam{Float64}
    acentricfactor::SingleParam{Float64}
end

@newmodelsimple tcTwuAlpha TwuAlphaModel tcTwuAlphaParam
default_references(::Type{tcTwuAlpha}) = ["10.1016/0378-3812(80)80003-3","10.1016/j.fluid.2016.09.003","10.1021/acs.jced.8b00640","10.1002/aic.17518","10.1021/acs.iecr.1c03003"]
default_ignore_missing_singleparams(::Type{tcTwuAlpha}) = ["N","N","L","acentricfactor","Tc","Pc","Vc","ZRA"]

"""
    tcTwuAlpha::TwuAlpha

    tcTwuAlpha(components::Vector{String};
    userlocations = String[],
    verbose::Bool=false)

## Input Parameters
- `M`: Single Parameter (`Float64`)
- `N`: Single Parameter (`Float64`)
- `L`: Single Parameter (`Float64`)

## Model Parameters
- `M`: Single Parameter (`Float64`)
- `N`: Single Parameter (`Float64`)
- `L`: Single Parameter (`Float64`)

## Description
Twu alpha `(α(T))` model, used in the tc-PR and tc-RK models.

```
αᵢ = Trᵢ^(N*(M-1))*exp(L*(1-Trᵢ^(N*M))
Trᵢ = T/Tcᵢ
N = 2
```

For PR, if no Twu parameters are provided (From the 2022 paper (default), both the 2018 and 2016 versions are also available):

2022 version :

```
L(ω) = 0.0297*ω² + 0.7536*ω + 0.0544
M(ω) = 0.1401*ω² - 0.1785*ω + 0.8678

```

2018 version :

```
L(ω) = 0.0925*ω² + 0.6693*ω + 0.0728
M(ω) = 0.1695*ω² - 0.2258*ω + 0.8788

```

2016 version :

```
L(ω) = 0.1290*ω² + 0.6039*ω + 0.877
M(ω) = 0.1760*ω² - 0.2600*ω + 0.8884

```

For RK, if no Twu parameters are provided (From the 2018 paper (default), the 2016 version is also available):

2018 version :

```
L(ω) = 0.0611*ω² + 0.7535*ω + 0.1359
M(ω) = 0.1709*ω² - 0.2063*ω + 0.8787

```

2016 version :

```
L(ω) = 0.0947*ω² + 0.6871*ω + 0.1508
M(ω) = 0.1615*ω² - 0.2349*ω + 0.8876

```


## Model Construction Examples
```
# Using the default database
alpha = tcTwuAlpha("water") #single input
alpha = tcTwuAlpha(["water","ethanol"]) #multiple components

# Using user-provided parameters

# Passing files or folders
alpha = tcTwuAlpha(["neon","hydrogen"]; userlocations = ["path/to/my/db","twu88.csv"])

# Passing parameters directly
alpha = tcTwuAlpha(["neon","hydrogen"];
    userlocations = (;L = [0.40453, 156.21],
                    M = [0.95861, -0.0062072],
                    N = [0.8396, 5.047]) #if we don't pass N, then is assumed N = 2
                )
```

## References
1. Twu, C. H., Lee, L. L., & Starling, K. E. (1980). Improved analytical representation of argon thermodynamic behavior. Fluid Phase Equilibria, 4(1–2), 35–44. [doi:10.1016/0378-3812(80)80003-3](https://doi.org/10.1016/0378-3812(80)80003-3)
2. Le Guennec, Y., Privat, R., & Jaubert, J.-N. (2016). Development of the translated-consistent tc-PR and tc-RK cubic equations of state for a safe and accurate prediction of volumetric, energetic and saturation properties of pure compounds in the sub- and super-critical domains. Fluid Phase Equilibria, 429, 301–312. [doi:10.1016/j.fluid.2016.09.003](http://dx.doi.org/10.1016/j.fluid.2016.09.003)
3. Pina-Martinez, A., Le Guennec, Y., Privat, R., Jaubert, J.-N., & Mathias, P. M. (2018). Analysis of the combinations of property data that are suitable for a safe estimation of consistent twu α-function parameters: Updated parameter values for the translated-consistent tc-PR and tc-RK cubic equations of state. Journal of Chemical and Engineering Data, 63(10), 3980–3988. [doi:10.1021/acs.jced.8b00640](http://dx.doi.org/10.1021/acs.jced.8b00640)
4. Piña-Martinez, A., Privat, R., & Jaubert, J.-N. (2022). Use of 300,000 pseudo‐experimental data over 1800 pure fluids to assess the performance of four cubic equations of state: SRK , PR , tc ‐RK , and tc ‐PR. AIChE Journal. American Institute of Chemical Engineers, 68(2). [doi:10.1002/aic.17518](https://doi.org/10.1002/aic.17518)

"""
tcTwuAlpha

function transform_params(::Type{tcTwuAlpha},params,components)
    nc = length(components)
    #just initialization
    N = get!(params,"N") do
        SingleParam("N",components,zeros(nc),fill(true,nc))
    end

    M = get!(params,"M") do
        SingleParam("N",components,zeros(nc),fill(true,nc))
    end

    L = get!(params,"L") do
        SingleParam("L",components,zeros(nc),fill(true,nc))
    end

    acentricfactor = get!(params,"acentricfactor") do
        SingleParam("acentricfactor",components,zeros(nc),fill(true,nc))
    end
    return params
end

function recombine_alpha!(model::PRModel,alpha::tcTwuAlpha)
    N = alpha.params.N
    M = alpha.params.M
    L = alpha.params.L
    w = alpha.params.acentricfactor
    nc = length(model)
    for i in 1:nc
        L.ismissingvalues[i] && w.ismissingvalues[i] && throw(error("tcTwuAlpha: cannot estimate L: missing acentricfactor parameter"))
        M.ismissingvalues[i] && w.ismissingvalues[i] && throw(error("tcTwuAlpha: cannot estimate M: missing acentricfactor parameter"))
        N.ismissingvalues[i] && w.ismissingvalues[i] && throw(error("tcTwuAlpha: cannot estimate N: missing acentricfactor parameter"))
        ω = w[i]

        if N.ismissingvalues[i] && !(w.ismissingvalues[i])
            N[i] = 2
        end

        if L.ismissingvalues[i] && !(w.ismissingvalues[i])
            L[i] = evalpoly(ω,(0.0544,0.7536,0.0297)) #2021 Version
            #L[i] = evalpoly(ω,(0.0728,0.6693,0.0925)) #2018 Version
            #L[i] = evalpoly(ω,(0.877,0.6039,0.1290)) #2016 Version
        end

        if M.ismissingvalues[i] && !(w.ismissingvalues[i])
            M[i] = evalpoly(ω,(0.8678,-0.1785,0.1401)) #2021 Version
            #M[i] = evalpoly(ω,(0.8788,0.2258,0.1695)) #2018 Version
            #M[i] = evalpoly(ω,(0.8884,-0.2600,0.1760)) #2016 Version
        end       
    end
end

function recombine_alpha!(model::RKModel,alpha::tcTwuAlpha)
    N = alpha.params.N
    M = alpha.params.M
    L = alpha.params.L
    w = alpha.params.acentricfactor
    nc = length(model)
    for i in 1:nc
        L.ismissingvalues[i] && w.ismissingvalues[i] && throw(error("tcTwuAlpha: cannot estimate L: missing acentricfactor parameter"))
        M.ismissingvalues[i] && w.ismissingvalues[i] && throw(error("tcTwuAlpha: cannot estimate M: missing acentricfactor parameter"))
        N.ismissingvalues[i] && w.ismissingvalues[i] && throw(error("tcTwuAlpha: cannot estimate N: missing acentricfactor parameter"))
        ω = w[i]

        if N.ismissingvalues[i] && !(w.ismissingvalues[i])
            N[i] = 2
        end

        if L.ismissingvalues[i] && !(w.ismissingvalues[i])
            L[i] = evalpoly(ω,(0.1359,0.7535,0.0611)) #2018 Version
            #L[i] = evalpoly(ω,(0.1508,0.6871,0.0947)) #2016 Version
        end

        if M.ismissingvalues[i] && !(w.ismissingvalues[i])
            M[i] = evalpoly(ω,(0.8787,-0.2063,0.1709)) #2018 Version
            #M[i] = evalpoly(ω,(0.8876,-0.2349,0.1615)) #2018 Version
        end       
    end
end

export tcTwuAlpha