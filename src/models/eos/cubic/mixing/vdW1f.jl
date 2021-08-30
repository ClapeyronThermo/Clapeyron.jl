abstract type vdW1fRuleModel <: MixingRule end

struct vdW1fRuleParam <: EoSParam
end

@newmodelsimple vdW1fRule vdW1fRuleModel vdW1fRuleParam

export vdW1fRule
function vdW1fRule(components::Vector{String}; activity=nothing, userlocations::Vector{String}=String[], activity_userlocations::Vector{String}=String[], verbose::Bool=false)
    params = getparams(components, ["properties/critical.csv"]; userlocations=userlocations, verbose=verbose)
    acentricfactor = SingleParam(params["w"],"acentric factor")
    packagedparams = vdW1fRuleParam()
    model = vdW1fRule(packagedparams, verbose=verbose)
    return model
end

function mixing_rule(model::ABCubicModel,V,T,z,mixing_model::vdW1fRuleModel,α,a,b,c)
    n = sum(z)
    invn2 = (one(n)/n)^2
    b̄ = dot(z,Symmetric(b),z) * invn2
    c̄ = dot(z,c)/n
    ā = dot(z,Symmetric(a .* sqrt.(α*α')),z) * invn2
    return ā,b̄,c̄
end

is_splittable(::vdW1fRule) = false
