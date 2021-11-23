abstract type MHV2RuleModel <: MixingRule end

struct MHV2Rule{γ} <: MHV2RuleModel
    components::Array{String,1}
    activity::γ
end

@registermodel MHV2Rule
export MHV2Rule
function MHV2Rule(components::Vector{String}; activity = Wilson, userlocations::Vector{String}=String[],activity_userlocations::Vector{String}=String[], verbose::Bool=false)
    init_activity = activity(components;userlocations = activity_userlocations,verbose)
    
    model = MHV2Rule(components, init_activity)
    return model
end

MHV2q(::PRModel) = (-0.4347,-0.003654)
MHV2q(::RKModel) = (-0.4783,-0.0047 )
function mixing_rule(model::Union{PRModel,RKModel},V,T,z,mixing_model::MHV2RuleModel,α,a,b,c)
    n = sum(z)
    #x = z./n
    invn = (one(n)/n)
    invn2 = invn^2
    g_E = excess_gibbs_free_energy(mixing_model.activity,1e5,T,z) / n
    b̄ = dot(z,Symmetric(b),z) * invn2
    c̄ = dot(z,c)/n
    #ᾱ = a.*sqrt.(α.*α')./(b*R̄*T)
    q1,q2 = MHV2q(model)
    Σlogb = sum(z[i]*log(b̄/b[i,i]) for i ∈ @comps)
    Σab = sum(z[i]*a[i,i]*α[i]/b[i,i]/(R̄*T) for i ∈ @comps)*invn
    Σab2 = sum(z[i]*a[i,i]*(α[i]^2)/b[i,i]/(R̄*T) for i ∈ @comps)*invn
    c  = -q1*Σab-q2*Σab2-g_E/(R̄*T)-Σlogb
    ā = b̄*R̄*T*(-q1-sqrt(q1^2-4*q2*c))/(2*q2)
    return ā,b̄,c̄
end
#=
function mixing_rule(model::PRModel,V,T,z,mixing_model::MHV2RuleModel,α,a,b,c)
    n = sum(z)
    x = z./n
    invn2 = (one(n)/n)^2
    g_E = excess_gibbs_free_energy(mixing_model.activity,1e5,T,z) / n
    b̄ = dot(z,Symmetric(b),z) * invn2
    c̄ = dot(z,c)/n

    ᾱ = a.*sqrt.(α.*α')./(b*R̄*T)

    q1 = -0.4347
    q2 = -0.003654
    c  = -q1*sum(x[i]*ᾱ[i,i] for i ∈ @comps)-q2*sum(x[i]*ᾱ[i,i]^2 for i ∈ @comps)-g_E/(R̄*T)-sum(x[i]*log(b̄/b[i,i]) for i ∈ @comps)

    ā = b̄*R̄*T*(-q1-sqrt(q1^2-4*q2*c))/(2*q2)
    return ā,b̄,c̄
end=#