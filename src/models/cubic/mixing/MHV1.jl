abstract type MHV1RuleModel <: MixingRule end

struct MHV1Rule{γ} <: MHV1RuleModel
    components::Array{String,1}
    activity::γ
    references::Array{String,1}
end

@registermodel MHV1Rule
export MHV1Rule
function MHV1Rule(components::Vector{String}; activity = Wilson, userlocations::Vector{String}=String[],activity_userlocations::Vector{String}=String[], verbose::Bool=false)
    init_activity = activity(components;userlocations = activity_userlocations,verbose)
    references = ["10.1016/0378-3812(90)85053-D"]
    model = MHV1Rule(components, init_activity,references)
    return model
end

MHV1q(::PRModel) = 0.53
MHV1q(::RKModel) = 0.593
function mixing_rule(model::Union{RKModel,PRModel},V,T,z,mixing_model::MHV1RuleModel,α,a,b,c)
    n = sum(z)
    invn = (one(n)/n)
    invn2 = invn^2
    g_E = excess_gibbs_free_energy(mixing_model.activity,1e5,T,z) / n
    b̄ = dot(z,Symmetric(b),z) * invn2
    c̄ = dot(z,c)/n
    q = MHV1q(model)
    Σlogb = sum(z[i]*log(b̄/b[i,i]) for i ∈ @comps)*invn
    Σab = invn*sum(z[i]*a[i,i]*α[i]/b[i,i]/(R̄*T) for i ∈ @comps)
    ā = b̄*R̄*T*(Σab-1/q*(g_E/(R̄*T)+Σlogb))
    return ā,b̄,c̄
end