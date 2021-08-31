abstract type HVRuleModel <: MixingRule end

struct HVRule{γ} <: HVRuleModel
    components::Array{String,1}
    activity::γ
end

@registermodel HVRule
export HVRule
function HVRule(components::Vector{String}; activity = Wilson, userlocations::Vector{String}=String[],activity_userlocations::Vector{String}=String[], verbose::Bool=false)
    init_activity = activity(components;userlocations = activity_userlocations,verbose)
    
    model = HVRule(components, init_activity)
    return model
end

function mixing_rule(model::RKModel,V,T,z,mixing_model::HVRuleModel,α,a,b,c)
    n = sum(z)
    invn2 = (one(n)/n)^2
    b̄ = dot(z,Symmetric(b),z) * invn2
    c̄ = dot(z,c)/n
    ā = b̄*((sum(z[i]*a[i,i]*α[i]/b[i,i] for i ∈ @comps))/n-excess_gibbs_free_energy(mixing_model.activity,1e5,T,z)/n/log(2))
    return ā,b̄,c̄
end

function mixing_rule(model::PRModel,V,T,z,mixing_model::HVRuleModel,α,a,b,c)
    λ = 1/(2*√(2))*log((2+√(2))/(2-√(2)))
    n = sum(z)
    invn2 = (one(n)/n)^2
    b̄ = dot(z,Symmetric(b),z) * invn2
    c̄ = dot(z,c)/n
    ā = b̄*((sum(z[i]*a[i,i]*α[i]/b[i,i] for i ∈ @comps))/n-excess_gibbs_free_energy(mixing_model.activity,1e5,T,z)/n/λ)
    return ā,b̄,c̄
end

