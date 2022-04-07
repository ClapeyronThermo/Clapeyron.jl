abstract type PPR78RuleModel <: MixingRule end

struct PPR78Param <: EoSParam
    A::PairParam{Float64}
    B::PairParam{Float64}
end

struct PPR78Rule <: PPR78RuleModel
    groups::GroupParam
    components::Vector{String}
    params::PPR78Param
    references::Vector{String}
end
@registermodel PPR78Rule

"""
    PPR78Rule <: PPR78RuleModel
    
    PPR78Rule(components;
    userlocations::Vector{String}=String[],
    verbose::Bool=false)

## Input Parameters

- `A`: Pair Parameter (`Float64`) - Fitted Parameter `[K]`
- `B`: Pair Parameter (`Float64`) - Fitted Parameter `[K]`

## Description

PPR78 Mixing Rule 
```
aᵢⱼ = √(aᵢaⱼ)
bᵢⱼ = (bᵢ +bⱼ)/2
b̄ = ∑bᵢⱼxᵢxⱼ
c̄ = ∑cᵢxᵢ
ā = b̄(∑[xᵢaᵢᵢαᵢ/(bᵢᵢ)] - ∑xᵢxⱼbᵢbⱼEᵢⱼ/2b̄)
Eᵢⱼ = ∑(z̄ᵢₖ - z̄ⱼₖ)(z̄ᵢₗ - z̄ⱼₗ) × Aₖₗ × (298.15/T)^(Aₖₗ/Bₖₗ - 1)
```

## References
1. Jaubert, J.-N., Privat, R., & Mutelet, F. (2010). Predicting the phase equilibria of synthetic petroleum fluids with the PPR78 approach. AIChE Journal. American Institute of Chemical Engineers, 56(12), 3225–3235. doi:10.1002/aic.12232
2. Jaubert, J.-N., Qian, J.-W., Lasala, S., & Privat, R. (2022). The impressive impact of including enthalpy and heat capacity of mixing data when parameterising equations of state. Application to the development of the E-PPR78 (Enhanced-Predictive-Peng-Robinson-78) model. Fluid Phase Equilibria, (113456), 113456. doi:10.1016/j.fluid.2022.113456

"""
PPR78Rule

export PPR78Rule

function PPR78Rule(components; activity = nothing, userlocations::Vector{String}=String[],activity_userlocations::Vector{String}=String[], verbose::Bool=false)
    groups = GroupParam(components,["cubic/EPPR78/EPPR78_groups.csv"]; verbose=verbose)
    params = getparams(groups, ["cubic/EPPR78/EPPR78_unlike.csv"]; userlocations=userlocations)
    pkgparams = PPR78Param(params["A"],params["B"])
    references = ["10.1002/aic.12232","10.1016/j.fluid.2022.113456"]
    model = PPR78Rule(groups,groups.components,pkgparams,references)
    return model
end


function mixing_rule(model::CubicModel,V,T,z,mixing_model::PPR78Rule,α,a,b,c)
    n = sum(z)
    invn = 1/n
    invn2 = invn*invn
    T̄ = 298.15/T
    b̄ = dot(z,Symmetric(b),z) * invn2
    c̄ = dot(z,c)/n
    _0 = zero(T+first(z))
    gᴱ = _0

    A = mixing_model.params.A.values
    B = mixing_model.params.B.values
    groups = mixing_model.groups 
    gc = groups.i_flattenedgroups
    z̄n = groups.n_groups_cache
   
    for i ∈ @comps
        zni = z̄n[i]
        ∑zni⁻¹ = 1/sum(zni)
        bi = b[i,i]
        for j in 1:i-1
            znj = z̄n[j]
            ∑znj⁻¹ = 1/sum(znj)
            Eij = _0
            for k in gc
                αik =zni[k]* ∑zni⁻¹
                αjk =znj[k]* ∑znj⁻¹
                Δαk = (αik - αjk)
                for l in 1:k-1 #Δαk*Δαl = Δαl*Δαk
                    αil =zni[l]* ∑zni⁻¹
                    αjl =znj[l]* ∑znj⁻¹
                    Δαl = (αil - αjl)
                    Akl = A[k,l]
                    Bkl = B[k,l]
                    if !iszero(Akl)
                        Eij -= Δαk*Δαl*Akl*T̄^(Bkl/Akl - 1) # -1/2 * 2
                    end
                end
            end
            gᴱ += bi*b[j,j]*z[i]*z[j]*Eij
        end
    end
    gᴱ = 0.5*gᴱ*invn2/b̄
    ∑ab = sum(z[i]*a[i,i]*α[i]/b[i,i] for i ∈ @comps)*invn
    ā = b̄*(∑ab-gᴱ)
    return ā,b̄,c̄
end

