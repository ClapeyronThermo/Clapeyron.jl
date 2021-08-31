abstract type LCVMRuleModel <: MixingRule end

struct LCVMRule{γ} <: LCVMRuleModel
    components::Array{String,1}
    activity::γ
end

@registermodel LCVMRule
export LCVMRule
function LCVMRule(components::Vector{String}; activity = Wilson, userlocations::Vector{String}=String[],activity_userlocations::Vector{String}=String[], verbose::Bool=false)
    init_activity = activity(components;userlocations = activity_userlocations,verbose)
    
    model = LCVMRule(components, init_activity)
    return model
end


function mixing_rule(model::PRModel,V,T,z,mixing_model::LCVMRuleModel,α,a,b,c)
    n = sum(z)
    x = z./n
    invn2 = (one(n)/n)^2
    g_E = excess_gibbs_free_energy(mixing_model.activity,1e5,T,z) / n
    b̄ = dot(z,Symmetric(b),z) * invn2
    c̄ = dot(z,c)/n

    ᾱ = a.*sqrt.(α.*α')./(b*R̄*T)

    λ  = 0.7
    AV = -0.52
    AM = -0.623
    C1 = λ/AV+(1-λ)/AM

    ā = b̄*R̄*T*(C1*(g_E/(R̄*T)-0.3*sum(x[i]*log(b̄/b[i,i]) for i ∈ @comps))+sum(x[i]*ᾱ[i,i] for i ∈ @comps))
    return ā,b̄,c̄
end