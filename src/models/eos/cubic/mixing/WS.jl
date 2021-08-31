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

function mixing_rule(model::RKModel,V,T,z,mixing_model::WSRuleModel,α,a,b,c)
    n = sum(z)
    invn2 = (one(n)/n)^2
    num = sum(z[i]*(b[i,i]-a[i,i]*α[i]/(R̄*T)) for i ∈ @comps)/n
    den = 1 - (sum(z[i]*a[i,i]*α[i]/b[i,i]/(n*R̄*T) for i ∈ @comps)-excess_gibbs_free_energy(mixing_model.activity,1e5,T,z)/(n*R̄*T)/log(2))
    c̄ = dot(z,c)/n
    b̄  = num/den
    ā  = R̄*T*(b̄-num)
    return ā,b̄,c̄
end

function mixing_rule(model::PRModel,V,T,z,mixing_model::WSRuleModel,α,a,b,c)
    λ = 1/(2*√(2))*log((2+√(2))/(2-√(2)))
    n = sum(z)
    invn2 = (one(n)/n)^2
    num = sum(z[i]*(b[i,i]-a[i,i]*α[i]/(R̄*T)) for i ∈ @comps)/n
    den = 1 - (sum(z[i]*a[i,i]*α[i]/b[i,i]/(n*R̄*T) for i ∈ @comps)-excess_gibbs_free_energy(mixing_model.activity,1e5,T,z)/(n*R̄*T)/λ)
    c̄ = dot(z,c)/n
    b̄  = num/den
    ā  = R̄*T*(b̄-num)
    return ā,b̄,c̄
end

