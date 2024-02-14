struct NRTLParam <: EoSParam
    a::PairParam{Float64}
    b::PairParam{Float64}
    c::PairParam{Float64}
    Mw::SingleParam{Float64}
end

abstract type NRTLModel <: ActivityModel end

struct NRTL{c<:EoSModel} <: NRTLModel
    components::Array{String,1}
    params::NRTLParam
    puremodel::EoSVectorParam{c}
    references::Array{String,1}
end


export NRTL
"""
    NRTL <: ActivityModel

    function NRTL(components;
    puremodel=PR,
    userlocations = String[],
    pure_userlocations = String[],
    verbose = false)

## Input parameters
- `a`: Pair Parameter (`Float64`, asymetrical, defaults to `0`) - Interaction Parameter
- `b`: Pair Parameter (`Float64`, asymetrical, defaults to `0`) - Interaction Parameter
- `c`: Pair Parameter (`Float64`, asymetrical, defaults to `0`) - Interaction Parameter
- `Mw`: Single Parameter (`Float64`) (Optional) - Molecular Weight `[g/mol]`

## Input models
- `puremodel`: model to calculate pure pressure-dependent properties

## Description
NRTL (Non Random Two Fluid) activity model:
```
Gᴱ = nRT∑[xᵢ(∑τⱼᵢGⱼᵢxⱼ)/(∑Gⱼᵢxⱼ)]
Gᵢⱼ exp(-cᵢⱼτᵢⱼ)
τᵢⱼ = aᵢⱼ + bᵢⱼ/T
```

## Model Construction Examples
```
# Using the default database
model = NRTL(["water","ethanol"]) #Default pure model: PR
model = NRTL(["water","ethanol"],puremodel = BasicIdeal) #Using Ideal Gas for pure model properties
model = NRTL(["water","ethanol"],puremodel = PCSAFT) #Using Real Gas model for pure model properties

# Passing a prebuilt model

my_puremodel = AbbottVirial(["water","ethanol"]; userlocations = ["path/to/my/db","critical.csv"])
mixing = NRTL(["water","ethanol"],puremodel = my_puremodel)

# Using user-provided parameters

# Passing files or folders
model = NRTL(["water","ethanol"];userlocations = ["path/to/my/db","nrtl_ge.csv"])

# Passing parameters directly
model = NRTL(["water","ethanol"],
userlocations = (a = [0.0 3.458; -0.801 0.0],
                b = [0.0 -586.1; 246.2 0.0],
                c = [0.0 0.3; 0.3 0.0])
            )
```

## References
1. Renon, H., & Prausnitz, J. M. (1968). Local compositions in thermodynamic excess functions for liquid mixtures. AIChE journal. American Institute of Chemical Engineers, 14(1), 135–144. [doi:10.1002/aic.690140124](https://doi.org/10.1002/aic.690140124)
"""
NRTL

default_locations(::Type{NRTL}) = ["properties/molarmass.csv","Activity/NRTL/NRTL_unlike.csv"]

function NRTL(components; puremodel=PR,
    userlocations = String[], 
    pure_userlocations = String[],
    verbose = false)

    formatted_components = format_components(components)
    params = getparams(formatted_components, default_locations(NRTL); userlocations = userlocations, asymmetricparams=["a","b"], ignore_missing_singleparams=["a","b","Mw"], verbose = verbose)
    a  = params["a"]
    b  = params["b"]
    c  = params["c"]
    Mw  = get(params,"Mw",SingleParam("Mw",formatted_components))
    
    _puremodel = init_puremodel(puremodel,components,pure_userlocations,verbose)
    packagedparams = NRTLParam(a,b,c,Mw)
    references = String["10.1002/aic.690140124"]
    model = NRTL(formatted_components,packagedparams,_puremodel,references)
    return model
end

#=
function activity_coefficient(model::NRTLModel,p,T,z)
    a = model.params.a.values
    b = model.params.b.values
    c = model.params.c.values
    x = z ./ sum(z)
    τ = @. a+b/T
    G = @. exp(-c*τ)
    lnγ = sum(x[j]*τ[j,:].*G[j,:] for j ∈ @comps)./sum(x[k]*G[k,:] for k ∈ @comps)+sum(x[j]*G[:,j]/sum(x[k]*G[k,j] for k ∈ @comps).*(τ[:,j] .-sum(x[m]*τ[m,j]*G[m,j] for m ∈ @comps)/sum(x[k]*G[k,j] for k ∈ @comps)) for j in @comps)
    return exp.(lnγ)
end
=#


function excess_g_res(model::NRTLModel,p,T,z)
    a = model.params.a.values
    b  = model.params.b.values
    c  = model.params.c.values
    _0 = zero(T+first(z))
    n = sum(z)
    invn = 1/n
    invT = 1/(T)
    res = _0 
    for i ∈ @comps
        ∑τGx = _0
        ∑Gx = _0
        xi = z[i]*invn
        for j ∈ @comps
            xj = z[j]*invn
            τji = a[j,i] + b[j,i]*invT
            Gji = exp(-c[j,i]*τji)
            Gx = xj*Gji
            ∑Gx += Gx
            ∑τGx += Gx*τji
        end
        res += xi*∑τGx/∑Gx
    end
    return n*res*R̄*T
end

excess_gibbs_free_energy(model::NRTLModel,p,T,z) = excess_g_res(model,p,T,z)