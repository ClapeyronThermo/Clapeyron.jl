abstract type KayRuleModel <: MixingRule end

struct KayRuleParam <: EoSParam
end

@newmodelsimple KayRule KayRuleModel KayRuleParam

"""
    KayRule <: KayRuleModel
    
    KayRule(components::Vector{String};
    userlocations::Vector{String}=String[],
    verbose::Bool=false)

## Input Parameters

None

## Description

Kay mixing rule for cubic parameters:

```
aᵢⱼ = √(aᵢaⱼ)(1-kᵢⱼ)
bᵢⱼ = (bᵢ + bⱼ)/2
ā = b̄*(∑[aᵢⱼxᵢxⱼ√(αᵢ(T)αⱼ(T))/bᵢⱼ])^2
b̄ = (∑∛(bᵢⱼ)xᵢxⱼ)^3
c̄ = ∑cᵢxᵢ
```
"""
KayRule

export KayRule

function KayRule(components::Vector{String}; activity=nothing, userlocations::Vector{String}=String[], activity_userlocations::Vector{String}=String[], verbose::Bool=false, kwargs...)
    #params = getparams(components, ["properties/critical.csv"]; userlocations=userlocations, verbose=verbose)
    #acentricfactor = params["acentricfactor"]
    packagedparams = KayRuleParam()
    model = KayRule(packagedparams, verbose=verbose)
    return model
end

function mixing_rule(model::ABCubicModel,V,T,z,mixing_model::KayRuleModel,α,a,b,c)
    n = sum(z)
    invn = (one(n)/n)
    invn2 = invn^2
    #b̄ = (dot(z,Symmetric(b.^(1/3)),z) * invn2).^3
    #ā = √(dot(z,Symmetric(a .* sqrt.(α*α') ./ b).^2,z)* invn2)  * b̄
    b̄ = zero(eltype(z))
    ab2 = zero(T+first(z))
    for i in 1:length(z)
        zi = z[i]
        zi2 = zi*zi
        αi = α[i]
        b̄ += cbrt(b[i,i])*zi2
        abij2 = (a[i,i]*αi/b[i,i])^2
        ab2 += abij2*zi2
        for j in 1:(i-1)
            zij =zi*z[j]
            abij2 = (a[i,j]*sqrt(αi*α[j])/b[i,j])^2
            ab2 += 2*abij2*zij
            b̄ += 2*cbrt(b[i,j])*zij
        end
    end    
    b̄ = (b̄*invn2)^3
    ā = sqrt(ab2*invn2)*b̄
    c̄ = dot(z,c)*invn
    return ā,b̄,c̄
end

is_splittable(::KayRule) = false
