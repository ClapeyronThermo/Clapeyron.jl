abstract type KayRuleModel <: MixingRule end

struct KayRuleParam <: EoSParam
end

@newmodelsimple KayRule KayRuleModel KayRuleParam

export KayRule

function KayRule(components::Vector{String}; activity=nothing, userlocations::Vector{String}=String[], activity_userlocations::Vector{String}=String[], verbose::Bool=false)
    #params = getparams(components, ["properties/critical.csv"]; userlocations=userlocations, verbose=verbose)
    #acentricfactor = SingleParam(params["w"],"acentric factor")
    packagedparams = KayRuleParam()
    model = KayRule(packagedparams, verbose=verbose)
    return model
end

function mixing_rule(model::ABCubicModel,V,T,z,mixing_model::KayRuleModel,α,a,b,c)
    n = sum(z)
    invn2 = (one(n)/n)^2
    c̄ = dot(z,c)/n
    #b̄ = (dot(z,Symmetric(b.^(1/3)),z) * invn2).^3
    #ā = √(dot(z,Symmetric(a .* sqrt.(α*α') ./ b).^2,z)* invn2)  * b̄
    
    b̄ = zero(eltype(z))
    for i in 1:length(z)
        zi = z[i]
        b̄ += (cbrt(b[i,i])*invn2*zi^2)^3
        for j in 1:(i-1)
            b̄ += 2*(cbrt(b[i,j])*invn2*zi*z[j])^3
        end
    end

    ā = zero(eltype(α))+zero(eltype(z))
    for i in 1:length(z)
        zi = z[i]
        αi = α[i]
        ā += a[i,i]*αi/b[i,i]*zi*b̄
        for j in 1:(i-1)
            ā += 2*sqrt(zi*z[j]*invn2*(a[i,j]*sqrt(αi*α[j])/b[i,j])^2)*b̄
        end
    end
    ā
    return ā,b̄,c̄
end

is_splittable(::KayRule) = false
