abstract type LCVMRuleModel <: MixingRule end

struct LCVMRule{γ} <: LCVMRuleModel
    components::Array{String,1}
    activity::γ
    references::Array{String,1}
end

@registermodel LCVMRule
export LCVMRule
function LCVMRule(components::Vector{String}; activity = Wilson, userlocations::Vector{String}=String[],activity_userlocations::Vector{String}=String[], verbose::Bool=false)
    init_activity = activity(components;userlocations = activity_userlocations,verbose)
    references = ["10.1016/0378-3812(94)80043-X"]
    model = LCVMRule(components, init_activity,references)
    return model
end


function mixing_rule(model::PRModel,V,T,z,mixing_model::LCVMRuleModel,α,a,b,c)
    n = sum(z)
    #x = z./n
    invn = (one(n)/n)
    invn2 = invn^2
    g_E = excess_gibbs_free_energy(mixing_model.activity,1e5,T,z) / n
    b̄ = dot(z,Symmetric(b),z) * invn2
    c̄ = dot(z,c)/n
    #ᾱ  = a.*sqrt.(α.*α')./(b*R̄*T)
    Σxᾱ  = sum(α[i]*a[i,i]*z[i]/b[i,i] for i ∈ @comps)*invn/(R̄*T)
    λ  = 0.7
    AV = -0.52
    AM = -0.623
    C1 = -1.8276947771329795 #λ/AV+(1-λ)/AM,is this ok?, that C1 is independent of the input conditions
    Σlogb = sum(z[i]*log(b̄/b[i,i]) for i ∈ @comps)*invn #sum(x[i]*log(b̄/b[i,i]) for i ∈ @comps)
    ā = b̄*R̄*T*(C1*(g_E/(R̄*T)-0.3*Σlogb)+Σxᾱ )
    return ā,b̄,c̄
end