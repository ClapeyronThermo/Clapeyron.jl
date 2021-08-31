abstract type UMRRuleModel <: MixingRule end

struct UMRRule{γ} <: UMRRuleModel
    components::Array{String,1}
    activity::γ
end

@registermodel UMRRule
export UMRRule
function UMRRule(components::Vector{String}; activity = Wilson, userlocations::Vector{String}=String[],activity_userlocations::Vector{String}=String[], verbose::Bool=false)
    init_activity = activity(components;userlocations = activity_userlocations,verbose)
    
    model = UMRRule(components, init_activity)
    return model
end

function mixing_rule(model::PRModel,V,T,z,mixing_model::UMRRuleModel,α,a,b,c)
    n = sum(z)
    x = z./n
    invn2 = (one(n)/n)^2
    lnγ_SG_  = lnγ_SG(mixing_model.activity,1e5,T,z)
    lnγ_res_ = lnγ_res(mixing_model.activity,1e5,T,z)
    g_E = sum(x[i]*R̄*T*(lnγ_res_[i]+lnγ_SG_[i]) for i ∈ @comps)
    b = Diagonal(b).diag
    b = ((b.^(1/2).+b'.^(1/2))/2).^2
    b̄ = dot(z,Symmetric(b),z) * invn2
    c̄ = dot(z,c)/n
    ā = b̄*R̄*T*(sum(x[i]*a[i,i]*α[i]/b[i,i]/(R̄*T) for i ∈ @comps)-1/0.53*g_E/(R̄*T))
    return ā,b̄,c̄
end