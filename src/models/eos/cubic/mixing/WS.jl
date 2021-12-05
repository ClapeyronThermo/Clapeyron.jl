abstract type WSRuleModel <: MixingRule end

struct WSRule{γ} <: WSRuleModel
    components::Array{String,1}
    activity::γ
    references::Array{String,1}
end

@registermodel WSRule
export WSRule
function WSRule(components::Vector{String}; activity = Wilson, userlocations::Vector{String}=String[],activity_userlocations::Vector{String}=String[], verbose::Bool=false)
    init_activity = activity(components;userlocations = activity_userlocations,verbose)
    
    references = ["10.1002/aic.690380505"]
    model = WSRule(components, init_activity,references)
    return model
end

WS_λ(::PRModel) = 0.6232252401402305 #1/(2*√(2))*log((2+√(2))/(2-√(2)))
WS_λ(::RKModel) = 0.6931471805599453#log(2)

function mixing_rule(model::Union{RKModel,PRModel},V,T,z,mixing_model::WSRuleModel,α,a,b,c)
    λ = WS_λ(model)
    n = sum(z)
    invn = (one(n)/n)
    Σab = sum(z[i]*a[i,i]*α[i]/b[i,i]/(n*R̄*T) for i ∈ @comps)
    num = sum(z[i]*(b[i,i]-a[i,i]*α[i]/(R̄*T)) for i ∈ @comps)*invn
    gE = excess_gibbs_free_energy(mixing_model.activity,1e5,T,z)
    den = 1 - (Σab-gE/(n*R̄*T)/λ)
    c̄ = dot(z,c)*invn
    b̄  = num/den
    ā  = R̄*T*(b̄-num)
    return ā,b̄,c̄
end

