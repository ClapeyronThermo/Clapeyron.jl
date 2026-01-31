abstract type PPR78RuleModel <: MixingRule end

struct PPR78Param <: EoSParam
    A::PairParam{Float64}
    B::PairParam{Float64}
end

struct PPR78Rule{G} <: PPR78RuleModel
    groups::GroupParam{G}
    components::Vector{String}
    params::PPR78Param
    references::Vector{String}
end

"""
    PPR78Rule <: PPR78RuleModel
    
    PPR78Rule(components;
    userlocations = String[],
    group_userlocations = String[]
    verbose::Bool=false)
## Input Parameters
- `A`: Pair Parameter (`Float64`) - Fitted Parameter `[K]`
- `B`: Pair Parameter (`Float64`) - Fitted Parameter `[K]`
## Description
PPR78 Mixing Rule, Uses E-PPR78 Group params. Default for [`EPPR78`](@ref) EoS.
```
aᵢⱼ = √(aᵢaⱼ)
bᵢⱼ = (bᵢ +bⱼ)/2
b̄ = ∑bᵢⱼxᵢxⱼ
c̄ = ∑cᵢxᵢ
ā = b̄(∑[xᵢaᵢαᵢ/(bᵢᵢ)] - ∑xᵢxⱼbᵢbⱼEᵢⱼ/2b̄)
Eᵢⱼ = ∑(z̄ᵢₖ - z̄ⱼₖ)(z̄ᵢₗ - z̄ⱼₗ) × Aₖₗ × (298.15/T)^(Aₖₗ/Bₖₗ - 1)
```
## References
1. Jaubert, J.-N., Privat, R., & Mutelet, F. (2010). Predicting the phase equilibria of synthetic petroleum fluids with the PPR78 approach. AIChE Journal. American Institute of Chemical Engineers, 56(12), 3225–3235. [doi:10.1002/aic.12232](https://doi.org/10.1002/aic.12232)
2. Jaubert, J.-N., Qian, J.-W., Lasala, S., & Privat, R. (2022). The impressive impact of including enthalpy and heat capacity of mixing data when parameterising equations of state. Application to the development of the E-PPR78 (Enhanced-Predictive-Peng-Robinson-78) model. Fluid Phase Equilibria, (113456), 113456. [doi:10.1016/j.fluid.2022.113456](https://doi.org/10.1016/j.fluid.2022.113456)
"""
PPR78Rule

export PPR78Rule

function PPR78Rule(components;
    activity = nothing,
    userlocations = String[],
    group_userlocations = String[],
    activity_userlocations = String[],
    verbose::Bool=false)
    
    _components = format_gccomponents(components)
    groups = GroupParam(_components,["cubic/EPPR78/EPPR78_groups.csv"]; group_userlocations = group_userlocations,verbose = verbose)
    params = getparams(groups, ["cubic/EPPR78/EPPR78_unlike.csv"]; userlocations = userlocations, verbose = verbose, ignore_missing_singleparams=["A","B"])
    pkgparams = PPR78Param(params["A"],params["B"])
    references = ["10.1002/aic.12232","10.1016/j.fluid.2022.113456"]
    model = PPR78Rule(groups,groups.components,pkgparams,references)
    return model
end

recombine_impl!(model::PPR78Rule) = model

function mixing_rule(model::CubicModel,V,T,z,mixing_model::PPR78Rule,α,a,b)
    n = sum(z)
    invn = 1/n
    invn2 = invn*invn
    T̄ = 298.15/T
    b̄ = dot(z,diagvalues(b)) * invn
    c̄ = translation2(model,V,T,z,model.translation,a,b,α)*invn
    _0 = zero(T+first(z))
    gᴱ = _0

    A = mixing_model.params.A.values
    B = mixing_model.params.B.values
    groups = mixing_model.groups 
    gc = 1:length(groups.flattenedgroups)
    z̄n = groups.n_flattenedgroups
   
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
            gᴱ += bi*b[j,j]*z[i]*z[j]*Eij #(0.5 * 2)
        end
    end
    gᴱ = gᴱ*invn2/b̄*1e6
    ∑ab = sum(z[i]*a[i,i]*α[i]/b[i,i] for i ∈ @comps)*invn
    ā = b̄*(∑ab-gᴱ)
    return ā,b̄,c̄
end

#=
groups:
["CH3",
"CH2",
"CH",
"C",
"CH4",
"C2H6",
"CH aro",
"C aro",
"C fused aromatic rings",
"CH2 cyclic",
"CH cyclic~|~C cyclic",
"CO2",
"N2",
"H2S",
"SH",
"H2O",
"C2H4",
"CH2 alkenic~|~CH alkenic",
"C alkenic",
"CH cycloalkenic~|~C cycloalkenic",
"H2",
"C2F6",
"CF3",
"CF2",
"CF2 double bond~|~CF double bond",
"C2H4F2",
"C2H2F4",
"CO",
"He",
"Ar",
"SO2",
"O2",
"NO",
"COS",
"NH3",
"NO2~|~N2O4",
"N2O",
"C2H2",
"HC=-C-",
"-C=-C-"]
=#