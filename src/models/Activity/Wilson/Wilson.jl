struct WilsonParam <: EoSParam
    g::PairParam{Float64}
    Tc::SingleParam{Float64}
    Pc::SingleParam{Float64}
    ZRA::SingleParam{Float64}
    Mw::SingleParam{Float64}
end

abstract type WilsonModel <: ActivityModel end

struct Wilson{c<:EoSModel} <: WilsonModel
    components::Array{String,1}
    params::WilsonParam
    puremodel::EoSVectorParam{c}
    references::Array{String,1}
end

export Wilson

"""
    Wilson <: ActivityModel
    Wilson(components;
    puremodel = PR,
    userlocations = String[],
    pure_userlocations = String[],
    verbose = false,
    reference_state = nothing)

## Input parameters
- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Pressure `[Pa]`
- `acentricfactor`: Single Parameter (`Float64`) - Acentric Factor
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `g`: Pair Parameter (`Float64`, asymetrical, defaults to `0`) - Binary Interaction Parameter

## model parameters
- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Pressure `[Pa]`
- `ZRA`: Single Parameter (`Float64`) - Rackett compresibility factor
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `g`: Pair Parameter (`Float64`, asymetrical, defaults to `0`) - Binary Interaction Parameter

## Input models
- `puremodel`: model to calculate pure pressure-dependent properties

## Description
Wilson activity model, with Rackett correlation for liquid volume:
```
Gᴱ = nRT∑xᵢlog(∑xⱼjΛᵢⱼ)
Λᵢⱼ = exp(-gᵢⱼ/T)*Vⱼ/Vᵢ
ZRAᵢ = 0.29056 - 0.08775ωᵢ
Vᵢ = (RTcᵢ/Pcᵢ)ZRAᵢ^(1 + (1-T/Tcᵢ)^2/7)
```

## Model Construction Examples
```
# Using the default database
model = Wilson(["water","ethanol"]) #Default pure model: PR
model = Wilson(["water","ethanol"],puremodel = BasicIdeal) #Using Ideal Gas for pure model properties
model = Wilson(["water","ethanol"],puremodel = PCSAFT) #Using Real Gas model for pure model properties

# Passing a prebuilt model

my_puremodel = AbbottVirial(["water","ethanol"]; userlocations = ["path/to/my/db","critical.csv"])
mixing = Wilson(["water","ethanol"],puremodel = my_puremodel)

# Using user-provided parameters

# Passing files or folders
model = Wilson(["water","ethanol"];userlocations = ["path/to/my/db","wilson.csv"])

# Passing parameters directly
model = Wilson(["water","ethanol"],
        userlocations = (g = [0.0 3988.52; 1360.117 0.0],
                        Tc = [647.13, 513.92],
                        Pc = [2.19e7, 6.12e6],
                        acentricfactor = [0.343, 0.643],
                        Mw = [18.015, 46.069])
                        )
```

## References
1. Wilson, G. M. (1964). Vapor-liquid equilibrium. XI. A new expression for the excess free energy of mixing. Journal of the American Chemical Society, 86(2), 127–130. [doi:10.1021/ja01056a002](https://doi.org/10.1021/ja01056a002)
"""
Wilson

default_locations(::Type{Wilson}) = ["properties/critical.csv", "properties/molarmass.csv","Activity/Wilson/Wilson_unlike.csv"]

function Wilson(components;
    puremodel = PR,
    userlocations = String[],
    pure_userlocations = String[],
    verbose = false,
    reference_state = nothing)
    
    formatted_components = format_components(components)
    params = getparams(formatted_components, default_locations(Wilson); userlocations = userlocations, asymmetricparams=["g"], ignore_missing_singleparams=["g"], verbose = verbose)
    g  = params["g"]
    Tc        = params["Tc"]
    pc        = params["Pc"]
    Mw        = params["Mw"]
    ZRA       = params["acentricfactor"]
    ZRA.values .*= -0.08775
    ZRA.values .+= 0.29056
    
    _puremodel = init_puremodel(puremodel,formatted_components,pure_userlocations,verbose)
    packagedparams = WilsonParam(g,Tc,pc,ZRA,Mw)
    references = String["10.1021/ja01056a002"]
    model = Wilson(formatted_components,packagedparams,_puremodel,references)
    set_reference_state!(model,reference_state,verbose = verbose)
    return model
end
#=
function activity_coefficient(model::WilsonModel,p,T,z)
    return activity_coefficient_wilson(model,p,T,z)
end

function activity_coefficient_wilson(model::WilsonModel,p,T,z,Vi = wilson_volume(model,T))
    Λ = (Vi' ./ Vi) .*exp.(-model.params.g.values/R̄/T)
    x = z ./ sum(z)
    lnγ = 1 .- log.(sum(x[i]*Λ[:,i] for i ∈ @comps)) .-sum(x[j] .*Λ[j,:] ./(sum(x[i]*Λ[j,i] for i ∈ @comps)) for j ∈ @comps)
    return exp.(lnγ)
end
=#

function excess_gibbs_free_energy(model::WilsonModel,p,T,z)
    excess_g_wilson(model::WilsonModel,p,T,z)
end

function excess_g_res(model::WilsonModel,p,T,z)
    excess_g_res_wilson(model,p,T,z)
end

function excess_g_res_wilson(model::WilsonModel,p,T,z,V = wilson_volume(model,T))
    g_E = excess_g_wilson(model,p,T,z,V)
    g_comb = zero(g_E)
    zV = dot(z,V)
    zVinv = 1/zV
    for i in 1:length(model)
        g_comb += z[i]*(log(V[i]*zVinv))
    end
    g_comb = g_comb*Rgas(model)*T
    return g_E - g_comb
end

function excess_g_wilson(model::WilsonModel,p,T,z,V = wilson_volume(model,T))
    g = model.params.g.values
    _0 = zero(T+first(z))
    n = sum(z)
    invn = 1/n
    invRT = 1/(R̄*T)
    res = _0
    for i ∈ @comps
        ∑xΛ = _0
        xi = z[i]*invn
        for j ∈ @comps
            Λij = exp(-g[i,j]*invRT)*V[j]/V[i]
            ∑xΛ += Λij*z[j]*invn
        end
        res += xi*log(∑xΛ)
    end
    return -n*res*R̄*T
end

function wilson_volume(model::Wilson,T)
    #a^b^c is too slow to be done on a quadratic loop
    ZRA = model.params.ZRA.values
    Tc  = model.params.Tc.values
    Pc  = model.params.Pc.values
    V = zeros(typeof(1.0*T),length(model))
    for i ∈ @comps
        Tci = Tc[i]
        Tri = T/Tci
        V[i] = (R̄ *Tci/Pc[i])*ZRA[i]^(1 + (1-Tri)^2/7)
    end
    return V
end