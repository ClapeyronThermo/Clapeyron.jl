abstract type WSRuleModel <: MixingRule end

struct WSRule{γ} <: WSRuleModel
    components::Array{String,1}
    activity::γ
end

@registermodel WSRule
export WSRule
function WSRule(components::Vector{String}; activity = Wilson, userlocations::Vector{String}=String[],activity_userlocations::Vector{String}=String[], verbose::Bool=false)
    init_activity = activity(components;userlocations = activity_userlocations,verbose)
    
    model = WSRule(components, init_activity)
    return model
end

WS_λ(::PRModel) = 1/(2*√(2))*log((2+√(2))/(2-√(2))) #definitely check this
WS_λ(::RKModel) = log(2)

function mixing_rule(model::Union{RKModel,PRModel},V,T,z,mixing_model::WSRuleModel,α,a,b,c)
    λ = WS_λ(model)
    n = sum(z)
    invn = (one(n)/n)
    invn2 = invn^2
    Σab = sum(z[i]*a[i,i]*α[i]/b[i,i]/(R̄*T) for i ∈ @comps)*invn
    num = Σab
    den = 1 - (Σab-excess_gibbs_free_energy(mixing_model.activity,1e5,T,z)/(n*R̄*T)/λ)
    c̄ = dot(z,c)/n
    b̄  = num/den
    ā  = R̄*T*(b̄-num)
    return ā,b̄,c̄
end

