abstract type sPCSAFTModel <: PCSAFTModel end
@newmodel sPCSAFT sPCSAFTModel PCSAFTParam
default_references(::Type{sPCSAFT}) = ["10.1021/ie020753p"]
default_locations(::Type{sPCSAFT}) = ["SAFT/PCSAFT", "SAFT/PCSAFT/sPCSAFT"]
function transform_params(::Type{sPCSAFT},params)
    sigma = params["sigma"]
    sigma.values .*= 1E-10
    return saft_lorentz_berthelot(params)
end
export sPCSAFT

"""
    sPCSAFT <: PCSAFTModel

    sPCSAFT(components;
    idealmodel = BasicIdeal,
    userlocations = String[],
    ideal_userlocations = String[],
    reference_state = nothing,
    verbose = false,
    assoc_options = AssocOptions())

## Input parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Single Parameter (`Float64`) - Segment Diameter [`A°`]
- `epsilon`: Single Parameter (`Float64`) - Reduced dispersion energy  `[K]`
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Paramater (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m^3]`
## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Pair Parameter (`Float64`) - Mixed segment Diameter `[m]`
- `epsilon`: Pair Parameter (`Float64`) - Mixed reduced dispersion energy`[K]`
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume
## Input models
- `idealmodel`: Ideal Model
## Description
Simplified Perturbed-Chain SAFT (sPC-SAFT)
## References
1. von Solms, N., Michelsen, M. L., & Kontogeorgis, G. M. (2003). Computational and physical performance of a modified PC-SAFT equation of state for highly asymmetric and associating mixtures. Industrial & Engineering Chemistry Research, 42(5), 1098–1105. [doi:10.1021/ie020753p](https://doi.org/10.1021/ie020753p)
"""
sPCSAFT

function a_hc(model::sPCSAFTModel, V, T, z , _data = @f(data))
    _,_,_,_,η,m̄ = _data
    g_hs = (1-η/2)/(1-η)^3
    a_hs = (4η-3η^2)/(1-η)^2
    return m̄*a_hs - (m̄-1)*log(g_hs)
end

function g_hs(model::sPCSAFTModel, V, T, z,_data = @f(data))
    _,_,_,_,η,_ = _data
    return (1-η/2)/(1-η)^3
end

function a_hs(model::sPCSAFTModel, V, T, z)
    _,_,_,_,η,_ = _data
    return (4η-3η^2)/(1-η)^2
end

function Δ(model::sPCSAFTModel, V, T, z, i, j, a, b,_data = @f(data))
    ϵ_associjab = model.params.epsilon_assoc.values[i,j][a,b]
    κijab = model.params.bondvol.values[i,j][a,b]
    σij = model.params.sigma.values[i,j]
    g_hs_ = @f(g_hs,_data)
    return g_hs_*σij^3*(exp(ϵ_associjab/T)-1)*κijab
end

#custom method for sPCSAFT, we only calculate g_hs once
function  Δ(model::sPCSAFT, V, T, z,_data=@f(data))
    ϵ_assoc = model.params.epsilon_assoc.values
    κ = model.params.bondvol.values
    σ = model.params.sigma.values
    Δout = assoc_similar(κ,typeof(V+T+first(z)+one(eltype(model))))
    g_hs_ = @f(g_hs,_data)
    Δout.values .= false
    for (idx,(i,j),(a,b)) in indices(Δout)
        κijab = κ[idx]
        if κijab != 0
            Δout[idx] = g_hs_*σ[i,j]^3*(expm1(ϵ_assoc[i,j][a,b]/T))*κ[idx]
        end
    end
    return Δout
end