abstract type PSRKRuleModel <: MixingRule end

struct PSRKRule{γ} <: PSRKRuleModel
    components::Array{String,1}
    activity::γ
    references::Array{String,1}
end

@registermodel PSRKRule
export PSRKRule
function PSRKRule(components::Vector{String}; activity = Wilson, userlocations::Vector{String}=String[],activity_userlocations::Vector{String}=String[], verbose::Bool=false)
    init_activity = activity(components;userlocations = activity_userlocations,verbose)
    
    references = String["10.1016/0378-3812(91)85038-V"]
    model = PSRKRule(components, init_activity,references)
    return model
end

function mixing_rule(model::RKModel,V,T,z,mixing_model::PSRKRuleModel,α,a,b,c)
    n = sum(z)
    invn = (one(n)/n)
    invn2 = invn^2
    g_E = excess_gibbs_free_energy(mixing_model.activity,1e5,T,z) / n
    b̄ = dot(z,Symmetric(b),z) * invn2
    c̄ = dot(z,c)/n
    Σlogb = sum(z[i]*log(b̄/b[i,i]) for i ∈ @comps)*invn
    Σab = sum(z[i]*a[i,i]*α[i]/b[i,i]/(R̄*T) for i ∈ @comps)*invn
    ā = b̄*R̄*T*(Σab-1/0.64663*(g_E/(R̄*T)+Σlogb))
    return ā,b̄,c̄
end