struct FloryHugginsParam <: EoSParam
    a::PairParam{Float64}         # Flory-Huggins interaction parameter
    b::PairParam{Float64}       # Flory-Huggins volume parameter
    Mw::SingleParam{Float64}        # Molecular weight [g·mol⁻¹]
    N::SingleParam{Float64}         # Degree of polymerization
    v::SingleParam{Float64}         # Monomer volume (can be different for each component)
end

abstract type FloryHugginsModel <: ActivityModel end

struct FloryHuggins{c<:EoSModel} <: FloryHugginsModel
    components::Array{String,1}
    params::FloryHugginsParam
    puremodel::EoSVectorParam{c}
    references::Array{String,1}
end

const FH = FloryHuggins
export FH, FloryHuggins

"""
    FloryHuggins <: ActivityModel
    FloryHuggins(components, N;
    puremodel = BasicIdeal,
    userlocations = String[],
    pure_userlocations = String[],
    verbose = false,
    reference_state = nothing)
## Input parameters
- `N`: Single Parameter (`Float64`) - Degree of Polymerization
- `v`: Single Parameter (`Float64`) - Monomer Volume
- `Mw`: Single Parameter (`Float64`) - Molecular Weight
- `a`: Pair Parameter (`Float64`, defaults to `0`) - Interaction Parameter
- `b`: Pair Parameter (`Float64`, defaults to `0`) - Interaction Parameter

## Input models
- `puremodel`: model to calculate pure pressure-dependent properties
## Description
Flory-Huggins activity coefficient model:
```
Gᴱ = nRT·(∑[xᵢlog(rᵢ)]+N∑[ϕᵢϕⱼχᵢⱼ])
```
## References
1. Flory, P. J. (1953). "Principles of Polymer Chemistry". Cornell University Press.
2. Huggins, M. L. (1941). "Solutions of Long-Chain Compounds". Journal of Chemical Physics, 9(5), 440-440.
"""
FloryHuggins
default_locations(::Type{FloryHuggins}) = ["properties/molarmass.csv","Activity/FH/FH_unlike.csv","Activity/FH/FH_like.csv"]

function FloryHuggins(components, N; puremodel=BasicIdeal,
    userlocations = String[], 
    pure_userlocations = String[],
    verbose = false,
    reference_state = nothing)

    formatted_components = format_components(components)
    params = getparams(formatted_components, default_locations(FloryHuggins); userlocations = userlocations, ignore_missing_singleparams=["a","b","v","N"], verbose = verbose)

    a = params["a"]
    b = params["b"]
    Mw  = get(params,"Mw",SingleParam("Mw",formatted_components))
    N   = SingleParam("N", formatted_components, N)

    # Handle missing monomer volumes: set all to one if missing
    if haskey(params, "v") && !__ismissing(params["v"])
        v = params["v"]
        v.values .*= 1 ./Mw.values
    else
        v = SingleParam("v", formatted_components, fill(1.0, length(formatted_components)))
    end

    Mw.values .*= N.values

    _puremodel = init_puremodel(puremodel,components,pure_userlocations,verbose)
    packagedparams = FloryHugginsParam(a,b,Mw,N,v)
    references = String["10.1063/1.1723621"]
    model = FloryHuggins(formatted_components,packagedparams,_puremodel,references)
    set_reference_state!(model,reference_state,verbose = verbose)
    return model
end

# Excess Gibbs free energy for Flory-Huggins in molar units (per mole of mixture)
function excess_g_res(model::FloryHugginsModel, p, T, z)
    a = model.params.a.values
    b = model.params.b.values
    N = model.params.N.values
    v = model.params.v.values
    n = sum(z)
    
    x = z ./ n
    V = sum(x .* v .* N)
    NT = sum(x .* N)    
    v0 = V / NT
    ϕ = x .* v .* N ./ V  # Volume fraction of each component
    res = 0.0
    for i ∈ @comps
        xi = x[i]
        ϕi = ϕ[i]
        ri = v[i]*N[i]/V
        res += xi * log(ri)
        for j ∈ i+1:length(model.components)
            ϕj = ϕ[j]
            χij = a[i,j] + b[i,j]/T
            res += ϕi * ϕj * χij * NT
        end
    end
    return R̄ * T * res * n
end

excess_gibbs_free_energy(model::FloryHugginsModel, p, T, z) = excess_g_res(model, p, T, z)