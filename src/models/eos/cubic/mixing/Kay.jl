abstract type KayRuleModel <: MixingRule end

struct KayRuleParam <: EoSParam
end

@newmodelsimple KayRule KayRuleModel KayRuleParam

export KayRule

function KayRule(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false)
    #params = getparams(components, ["properties/critical.csv"]; userlocations=userlocations, verbose=verbose)
    #acentricfactor = SingleParam(params["w"],"acentric factor")
    packagedparams = KayRuleParam()
    model = KayRule(packagedparams, verbose=verbose)
    return model
end

function mixing_rule(model::ABCubicModel,V,T,z,mixing_model::KayRuleModel,α,a,b)
    n = sum(z)
    invn2 = (one(n)/n)^2
    b̄ = (dot(z,Symmetric(b.^(1/3)),z) * invn2).^3
    ā = √(dot(z,Symmetric(a .* (α*α') ./ b).^2,z)) * invn2 * b̄
    return ā,b̄
end

is_splittable(::KayRule) = false
