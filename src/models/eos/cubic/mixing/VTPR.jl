abstract type VTPRRuleModel <: MixingRule end

struct VTPRRule{γ} <: VTPRRuleModel
    components::Array{String,1}
    activity::γ
end

@registermodel VTPRRule
export VTPRRule
function VTPRRule(components::Vector{String}; activity = Wilson, userlocations::Vector{String}=String[],activity_userlocations::Vector{String}=String[], verbose::Bool=false)
    init_activity = activity(components;userlocations = activity_userlocations,verbose)
    
    model = VTPRRule(components, init_activity)
    return model
end

function mixing_rule(model::PRModel,V,T,z,mixing_model::VTPRRuleModel,α,a,b,c)
    n = sum(z)
    x = z./n
    invn2 = (one(n)/n)^2
    lnγ_res_ = lnγ_res(mixing_model.activity,1e5,T,z)
    g_E_res = sum(x[i]*R̄*T*lnγ_res_[i] for i ∈ @comps)
    b = Diagonal(b).diag
    b = ((b.^(3/4).+b'.^(3/4))/2).^(4/3)
    b̄ = dot(z,Symmetric(b),z) * invn2
    c̄ = dot(z,c)/n
    ā = b̄*R̄*T*(sum(x[i]*a[i,i]*α[i]/b[i,i]/(R̄*T) for i ∈ @comps)-1/0.53087*(g_E_res/(R̄*T)))
    return ā,b̄,c̄
end