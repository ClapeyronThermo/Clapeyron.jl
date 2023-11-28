struct QCPRRuleParam <: EoSParam
    A::SingleParam{Float64}
    B::SingleParam{Float64}
    l::PairParam{Float64}
end

abstract type QCPRRuleModel <: MixingRule end

struct QCPRRule <: QCPRRuleModel
    components::Array{String,1}
    params::QCPRRuleParam
    references::Array{String,1}
end

"""
    QCPRRule <: MHV2RuleModel
    
    QCPRRule(components;
    activity = Wilson,
    userlocations=String[],
    activity_userlocations=String[],
    verbose::Bool=false)
## Input Parameters
None
## Input models 
- `activity`: Activity Model
## Description
Quantum-Corrected Mixing Rule, used by [`QCPR`](@ref) EoS:
```
aᵢⱼ = √(aᵢaⱼ)(1 - kᵢⱼ)
bᵢⱼ = (1 - lᵢⱼ)(bqᵢ + bqⱼ)/2
bqᵢ = bᵢβᵢ(T)
βᵢ(T) = (1 + Aᵢ/(T + Bᵢ))^3 / (1 + Aᵢ/(Tcᵢ + Bᵢ))^3
ā = ∑aᵢⱼxᵢxⱼ√(αᵢ(T)αⱼ(T))
b̄ = ∑bᵢⱼxᵢxⱼ
c̄ = ∑cᵢxᵢ
```
## References
1. Aasen, A., Hammer, M., Lasala, S., Jaubert, J.-N., & Wilhelmsen, Ø. (2020). Accurate quantum-corrected cubic equations of state for helium, neon, hydrogen, deuterium and their mixtures. Fluid Phase Equilibria, 524(112790), 112790. [doi:10.1016/j.fluid.2020.112790](https://doi.org/10.1016/j.fluid.2020.112790)
"""
QCPRRule


function QCPRRule(components; activity = nothing, userlocations=String[],activity_userlocations=String[], verbose::Bool=false)
    params = getparams(components, ["cubic/QCPR/QCPR_like.csv","cubic/QCPR/QCPR_unlike.csv"]; userlocations=userlocations, verbose=verbose)
    references = String["10.1016/j.fluid.2020.112790"]
    pkgparams = QCPRRuleParam(params["A"],params["B"],params["l"])
    model = QCPRRule(format_components(components), pkgparams ,references)
    return model
end

recombine_impl!(model::QCPRRule) = model

function ab_premixing(model::PRModel,mixing::QCPRRuleModel,kij,lij)
    Tc = model.params.Tc
    Pc = model.params.Pc
    Ωa, Ωb = ab_consts(model)
    comps = Tc.components
    n = length(Tc)
    a = model.params.a
    b = model.params.b
    aii,bii = diagvalues(a),diagvalues(b)
    @. aii = Ωa*R̄^2*Tc^2/Pc
    @. bii = Ωb*R̄*Tc/Pc
    epsilon_LorentzBerthelot!(a,kij)
    sigma_LorentzBerthelot!(b)
    if lij !== nothing
        mixing.params.l.values .= lij
    end
    return a,b
end

function mixing_rule(model::PRModel,V,T,z,mixing_model::QCPRRuleModel,α,a,b,c)
    n = sum(z)
    invn = (one(n)/n)
    invn2 = invn^2
    A = mixing_model.params.A.values
    B = mixing_model.params.B.values
    l = mixing_model.params.l.values
    Tc = model.params.Tc.values
    ā = zero(T+first(z))
    b̄ = zero(first(z))
    for i in 1:length(z)
        zi = z[i]
        αi = α[i]
        zi2 = zi^2
        Bi = B[i]
        Ai = A[i]
        βi = (1 + Ai/(T + Bi))^3 / (1 + Ai/(Tc[i] + Bi))^3
        bqi = βi*b[i,i]
        b̄ += bqi*zi2
        ā += a[i,i]*αi*zi2
        for j in 1:(i-1)
            zij = zi*z[j]
            Bj = B[j]
            Aj = A[j]
            βj = (1 + Aj/(T + Bj))^3 / (1 + Aj/(Tc[j] + Bj))^3
            bqj = βj*b[j,j]
            ā += 2*a[i,j]*sqrt(αi*α[j])*zij
            b̄ += zij*(bqi+bqj)*(1-l[i,j]) #2 * zij * 0.5(bi + bj)
        end
    end
    ā *= invn2
    b̄ *= invn2
    c̄ = dot(z,c)*invn
    #dot(z,Symmetric(a .* sqrt.(α*α')),z) * invn2
    return ā,b̄,c̄
end

function cubic_get_k(model::CubicModel,mixing::QCPRRuleModel,params)
    return get_k_geomean(params.a.values)
end

function cubic_get_l(model::CubicModel,mixing::QCPRRuleModel,params)
    return copy(mixing.params.l.values)
end

export QCPRRule
