abstract type VTPRRuleModel <: MixingRule end

struct VTPRRule{γ} <: VTPRRuleModel
    components::Array{String,1}
    activity::γ
    references::Array{String,1}
end

@registermodel VTPRRule
export VTPRRule
function VTPRRule(components::Vector{String}; activity = UNIFAC, userlocations::Vector{String}=String[],activity_userlocations::Vector{String}=String[], verbose::Bool=false)
    init_activity = activity(components;userlocations = activity_userlocations,verbose)
    
    references = ["10.1016/S0378-3812(01)00626-4"]
    model = VTPRRule(components, init_activity,references)
    return model
end

function ab_premixing(::Type{PR},mixing::VTPRRule,Tc,pc,kij)
    Ωa, Ωb = ab_consts(PR)
    _Tc = Tc.values
    _pc = pc.values
    a = epsilon_LorentzBerthelot(SingleParam(pc, @. Ωa*R̄^2*_Tc^2/_pc),kij)
    bi = @. Ωb*R̄*_Tc/_pc
    bij = ((bi.^(3/4).+bi'.^(3/4))/2).^(4/3)
    b = PairParam("b",Tc.components,bij)
    return a,b
end

function mixing_rule(model::PRModel,V,T,z,mixing_model::VTPRRuleModel,α,a,b,c)
    n = sum(z)
    invn = (one(n)/n)
    invn2 = invn^2
    lnγ_res_ = lnγ_res(mixing_model.activity,1e5,T,z)
    g_E_res = sum(z[i]*R̄*T*lnγ_res_[i] for i ∈ @comps)*invn
    b̄ = dot(z,Symmetric(b),z) * invn2
    c̄ = dot(z,c)/n
    Σab = invn*sum(z[i]*a[i,i]*α[i]/b[i,i]/(R̄*T) for i ∈ @comps)
    ā = b̄*R̄*T*(Σab-1/0.53087*(g_E_res/(R̄*T)))
    return ā,b̄,c̄
end