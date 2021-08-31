abstract type PSRKRuleModel <: MixingRule end

struct PSRKRule{γ} <: PSRKRuleModel
    components::Array{String,1}
    activity::γ
end

@registermodel PSRKRule
export PSRKRule
function PSRKRule(components::Vector{String}; activity = Wilson, userlocations::Vector{String}=String[],activity_userlocations::Vector{String}=String[], verbose::Bool=false)
    init_activity = activity(components;userlocations = activity_userlocations,verbose)
    
    model = PSRKRule(components, init_activity)
    return model
end

function mixing_rule(model::RKModel,V,T,z,mixing_model::PSRKRuleModel,α,a,b,c)
    n = sum(z)
    x = z./n
    invn2 = (one(n)/n)^2
    g_E = excess_gibbs_free_energy(mixing_model.activity,1e5,T,z) / n
    b̄ = dot(z,Symmetric(b),z) * invn2
    c̄ = dot(z,c)/n
    ā = b̄*R̄*T*(sum(x[i]*a[i,i]*α[i]/b[i,i]/(R̄*T) for i ∈ @comps)-1/0.64663*(g_E/(R̄*T)+sum(x[i]*log(b̄/b[i,i]) for i ∈ @comps)))
    return ā,b̄,c̄
end