abstract type vdW1fRuleModel <: MixingRule end

struct vdW1fRuleParam <: EoSParam
end

@newmodelsimple vdW1fRule vdW1fRuleModel vdW1fRuleParam
export vdW1fRule

"""
    vdW1fRule <: vdW1fRuleModel
    
    vdW1fRule(components::Vector{String};
    userlocations::Vector{String}=String[],
    verbose::Bool=false)
## Input Parameters
None
## Description
van der Wals One-Fluid mixing rule for cubic parameters:
```
aᵢⱼ = √(aᵢaⱼ)(1-kᵢⱼ)
bᵢⱼ = (bᵢ + bⱼ)/2
ā = ∑aᵢⱼxᵢxⱼ√(αᵢ(T)αⱼ(T))
b̄ = ∑bᵢⱼxᵢxⱼ
c̄ = ∑cᵢxᵢ
```
"""
vdW1fRule

function vdW1fRule(components::Vector{String}; activity=nothing, userlocations::Vector{String}=String[], activity_userlocations::Vector{String}=String[], verbose::Bool=false)
    packagedparams = vdW1fRuleParam()
    model = vdW1fRule(packagedparams, verbose=verbose)
    return model
end

vdW1fRule() = vdW1fRule(vdW1fRuleParam())

function mixing_rule(model::ABCubicModel,V,T,z,mixing_model::vdW1fRuleModel,α,a,b,c)
    n = sum(z)
    invn = (one(n)/n)
    invn2 = invn^2
    #b̄ = dot(z,Symmetric(b),z) * invn2
    ā = zero(T+first(z))
    b̄ = zero(first(z))
    for i in 1:length(z)
        zi = z[i]
        αi = α[i]
        zi2 = zi^2
        b̄ += b[i,i]*zi2
        ā += a[i,i]*αi*zi^2
        for j in 1:(i-1)
            zij = zi*z[j]
            ā += 2*a[i,j]*sqrt(αi*α[j])*zij
            b̄ += 2*b[i,j]*zij
        end
    end
    ā *= invn2
    b̄ *= invn2
    c̄ = dot(z,c)*invn
    #dot(z,Symmetric(a .* sqrt.(α*α')),z) * invn2
    return ā,b̄,c̄
end



is_splittable(::vdW1fRule) = false