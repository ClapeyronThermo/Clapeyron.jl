struct UNIQUACParam <: EoSParam
    a::PairParam{Float64}
    r::SingleParam{Float64}
    q::SingleParam{Float64}
    q_p::SingleParam{Float64}
    Mw::SingleParam{Float64}
end

abstract type UNIQUACModel <: ActivityModel end

struct UNIQUAC{c<:EoSModel} <: UNIQUACModel
    components::Array{String,1}
    params::UNIQUACParam
    puremodel::EoSVectorParam{c}
    references::Array{String,1}
end

export UNIQUAC

"""
    UNIQUACModel <: ActivityModel

    UNIQUAC(components;
    puremodel = PR,
    userlocations = String[], 
    pure_userlocations = String[],
    verbose = false)

## Input parameters
- `a`: Pair Parameter (`Float64`, asymetrical, defaults to `0`) - Binary Interaction Energy Parameter
- `r`: Single Parameter (`Float64`)  - Normalized Van der Vals volume
- `q`: Single Parameter (`Float64`) - Normalized Surface Area
- `q_p`: Single Parameter (`Float64`) - Modified Normalized Surface Area 
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`


## Input models
- `puremodel`: model to calculate pure pressure-dependent properties

UNIQUAC (Universal QuasiChemical Activity Coefficients) activity model:

```
Gᴱ = nRT(gᴱ(comb) + gᴱ(res))
gᴱ(comb) = ∑[xᵢlog(Φᵢ/xᵢ) + 5qᵢxᵢlog(θᵢ/Φᵢ)]
gᴱ(res) = -∑xᵢqᵖᵢlog(∑θᵖⱼτⱼᵢ)
θᵢ = qᵢxᵢ/∑qᵢxᵢ
θᵖ = qᵖᵢxᵢ/∑qᵖᵢxᵢ
Φᵢ = rᵢxᵢ/∑rᵢxᵢ
τᵢⱼ = exp(-aᵢⱼ/T)
```

## Model Construction Examples
```
# Using the default database
model = UNIQUAC(["water","ethanol"]) #Default pure model: PR
model = UNIQUAC(["water","ethanol"],puremodel = BasicIdeal) #Using Ideal Gas for pure model properties
model = UNIQUAC(["water","ethanol"],puremodel = PCSAFT) #Using Real Gas model for pure model properties

# Passing a prebuilt model

my_puremodel = AbbottVirial(["water","ethanol"]; userlocations = ["path/to/my/db","critical.csv"])
mixing = UNIQUAC(["water","ethanol"],puremodel = my_puremodel)

# Using user-provided parameters

# Passing files or folders
model = UNIQUAC(["water","ethanol"];userlocations = ["path/to/my/db","uniquac_ge.csv"])

# Passing parameters directly
model = UNIQUAC(["water","ethanol"],
        userlocations = (a = [0.0 378.1; 258.4 0.0], 
                        r = [0.92, 2.11],
                        q = [1.4, 1.97],
                        q_p = [1.0, 0.92],
                        Mw = [18.015, 46.069])
                    )
```

## References

1. Abrams, D. S., & Prausnitz, J. M. (1975). Statistical thermodynamics of liquid mixtures: A new expression for the excess Gibbs energy of partly or completely miscible systems. AIChE journal. American Institute of Chemical Engineers, 21(1), 116–128. [doi:10.1002/aic.690210115](https://doi.org/10.1002/aic.690210115)

"""
UNIQUAC

default_locations(::Type{UNIQUAC}) = ["Activity/UNIQUAC/UNIQUAC_like.csv", "properties/molarmass.csv","Activity/UNIQUAC/UNIQUAC_unlike.csv"]

function UNIQUAC(components;
    puremodel = PR,
    userlocations = String[], 
    pure_userlocations = String[],
    verbose = false)

    formatted_components = format_components(components)
    params = getparams(formatted_components, default_locations(UNIQUAC); userlocations = userlocations, asymmetricparams=["a"], ignore_missing_singleparams=["a"], verbose = verbose)
    a  = params["a"]
    r  = params["r"]
    q  = params["q"]
    q_p = params["q_p"]
    Mw  = params["Mw"]
    
    _puremodel = init_puremodel(puremodel,components,pure_userlocations,verbose)
    packagedparams = UNIQUACParam(a,r,q,q_p,Mw)
    references = String[]
    model = UNIQUAC(formatted_components,packagedparams,_puremodel,references)
    return model
end
#=
function activity_coefficient(model::UNIQUACModel,p,T,z)
    q = model.params.q.values
    q_p = model.params.q_p.values
    r = model.params.r.values

    x = z ./ sum(z)

    Φ = x.*r/sum(x[i]*r[i] for i ∈ @comps)
    θ = x.*q/sum(x[i]*q[i] for i ∈ @comps)
    θ_p = x.*q_p/sum(x[i]*q_p[i] for i ∈ @comps)
    τ = Ψ(model,p,T,z)
    lnγ_comb = @. log(Φ/x)+(1-Φ/x)-5*q*(log(Φ/θ)+(1-Φ/θ))
    lnγ_res  = q_p.*(1 .-log.(sum(θ_p[i]*τ[i,:] for i ∈ @comps)) .-sum(θ_p[i]*τ[:,i]/sum(θ_p[j]*τ[j,i] for j ∈ @comps) for i ∈ @comps))
    return exp.(lnγ_comb+lnγ_res)
end
=#

function excess_g_comb(model::UNIQUACModel,p,T,z=SA[1.0])
    _0 = zero(eltype(z))
    r = model.params.r.values
    q = model.params.q.values
    
    n = sum(z)
    invn = 1/n
    Φm = dot(r,z)*invn
    θm = dot(q,z)*invn
    G_comp = _0
    for i ∈ @comps
        xi = z[i]*invn
        Φi = r[i]/Φm
        θi = q[i]/θm
        G_comp += xi*log(Φi) + 5*q[i]*xi*log(θi/Φi)
    end
    return n*G_comp
end

function excess_g_res(model::UNIQUACModel,p,T,z=SA[1.0])
    _0 = zero(T+first(z))
    q_p = model.params.q_p.values
    a = model.params.a.values
    n = sum(z)
    invn = 1/n
    invT = 1/T
    θpm = dot(q_p,z)*invn
    G_res = _0
    for i ∈ @comps
        q_pi = q_p[i]
        xi = z[i]*invn
        ∑θpτ = _0
        for j ∈ @comps
            θpj = q_p[j]*z[j]/θpm*invn
            τji = exp(-a[j,i]*invT)
            ∑θpτ += θpj*τji
        end
        G_res += q_pi*xi*log(∑θpτ)
    end
    return -n*G_res
end

function excess_gibbs_free_energy(model::UNIQUACModel,p,T,z)
    g_comp = excess_g_comb(model,p,T,z)
    g_res = excess_g_res(model,p,T,z)
    return (g_comp+g_res)*R̄*T 
end