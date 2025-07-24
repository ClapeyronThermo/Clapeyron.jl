abstract type YFRModel <: DeltaCubicModel end

#=

Mixing Rule

=#

struct YFR1fRule <: vdW1fRuleModel end

function YFR1fRule(components; activity = nothing, userlocations = String[],activity_userlocations = String[], verbose::Bool=false)
    YFR1fRule()
end

function mixing_rule(model::YFRModel,V,T,z,mixing_model::YFR1fRule,α,a,b,c)
    n = sum(z)
    invn = (one(n)/n)
    invn2 = invn^2
    Tc = model.params.Tc.values
    Vc = model.params.Vc.values
    #b̄ = dot(z,Symmetric(b),z) * invn2
    ā = zero(Base.promote_eltype(model,T,z))
    b̄ = zero(Base.promote_eltype(model,primalval(T),z))
    for i in 1:length(model)
        ai,bi = a[i,i],b[i,i]
        zi,αi,Tri,Zci = z[i],α[i],T/Tc[i],Vc[i]/bi
        zi2 = zi^2
        Ωai,Ωbi = YFR_Ωab(Zci,Tri)
        b̄ += Ωbi*bi*zi2
        ā += Ωai*ai*αi*zi2
        for j in 1:(i-1)
            bj = b[j,j]
            zj,αj,Trj,Zcj = z[j],α[j],T/Tc[j],Vc[j]/bj
            zij = zi*zj
            Ωaj,Ωbj = YFR_Ωab(Zcj,Trj)
            ā += 2*a[i,j]*sqrt(Ωai*Ωaj*αi*α[j])*zij
            b̄ += bj*Ωbj + Ωbi*bi
        end
    end
    ā *= invn2
    b̄ *= invn2
    c̄ = dot(z,c)*invn
    #dot(z,Symmetric(a .* sqrt.(α*α')),z) * invn2
    return ā,b̄,c̄
end

function mixing_rule1(model::YFRModel,V,T,z,mixing_model::YFR1fRule,α,a,b,c)
    Tc = model.params.Tc.values[1]
    Vc = model.params.Vc.values[1]
    Zc = Vc/b[1,1]
    Tr = T/Tc
    _1 = oneunit(z[1])
    Ωa,Ωb = YFR_Ωab(Zc,Tr)
    ā = Ωa*a[1,1]*α[1]*_1
    b̄ = Ωb*b[1,1]*_1
    c̄ = c[1]*_1
    return ā,b̄,c̄
end

function cubic_lb_volume(model, T, z, mixing::YFR1fRule)
    n = sum(z)
    invn = (one(n)/n)
    Tc = model.params.Tc.values
    Vc = model.params.Vc.values
    b = model.params.b.values
    b̄ = zero(Base.promote_eltype(model,primalval(T),z))
    for i in 1:length(model)
        bi = b[i,i]
        zi,Tri,Zci = z[i],T/Tc[i],Vc[i]/bi
        Ωbi = YFR_Ωb(Zci,Tri)
        b̄ += Ωbi*bi*zi
    end
    b̄ *= invn
    return b̄
end

#=

Alpha function

=#

abstract type YFRAlphaModel <: GeneralizedSuaveAlphaModel end

const YFRAlphaParam = SimpleAlphaParam

@newmodelsimple YFRAlpha YFRAlphaModel YFRAlphaParam
export YFRAlpha

@inline function α_m(model::YFRModel,alpha_model::YFRAlpha,i)
    Tc = model.params.Tc.values[i]
    Pc = model.params.Pc.values[i]
    Vc = model.params.Vc.values[i]
    Zc = Vc*Pc/(Rgas(model)*Tc)
    ω = alpha_model.params.acentricfactor.values[i]
    return 2.779200*Zc + 5.208803*Zc*ω − 0.314477
end

struct YFR{T <: IdealModel,α,c,M} <: YFRModel
    components::Array{String,1}
    alpha::α
    mixing::M
    translation::c
    params::ABCCubicParam
    idealmodel::T
    references::Array{String,1}
end

export YFR

"""
    YFR(components;
    idealmodel = BasicIdeal,
    alpha = YFRAlpha,
    mixing = YFRfRule,
    activity = nothing,
    translation = NoTranslation,
    userlocations = String[],
    ideal_userlocations = String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    reference_state = nothing,
    verbose = false)

## Input parameters
- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Pressure `[Pa]`
- `Vc`: Single Parameter (`Float64`) - Critical Volume `[m3/mol]`
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `k`: Pair Parameter (`Float64`) (optional)
- `l`: Pair Parameter (`Float64`) (optional)

## Model Parameters
- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Pressure `[Pa]`
- `Vc`: Single Parameter (`Float64`) - Critical Volume `[m3/mol]`
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `a`: Pair Parameter (`Float64`)
- `b`: Pair Parameter (`Float64`)
- `c`: Pair Parameter (`Float64`)


## Input models
- `idealmodel`: Ideal Model
- `alpha`: Alpha model
- `mixing`: Mixing model
- `activity`: Activity Model, used in the creation of the mixing model.
- `translation`: Translation Model

## Description

Yang-Frotscher-Richter Equation of state.

```
P = RT/(v-b) + a•α(T)/((v - Δ₁b)*(v - Δ₂b))
aᵢᵢ = Ωaᵢ(R²Tcᵢ²/Pcᵢ)
bᵢᵢ = Ωbᵢ(R²Tcᵢ/Pcᵢ)
cᵢ = Ωcᵢ(R²Tcᵢ/Pcᵢ)
Zcᵢ = Pcᵢ*Vcᵢ/(R*Tcᵢ)
Trᵢ = T/Tcᵢ
ξcᵢ = 0.144894*exp(-Trᵢ^4) − 0.129835*exp(-Trᵢ^3) + 0.957454*Zcᵢ + 0.036884
Ωaᵢ = 3ξcᵢ² + 3(1 - 2ξcᵢ)Ωbᵢ + Ωbᵢ² + 1 - 3ξcᵢ
Ωbᵢ: maximum real solution of 0 = -ξcᵢᵢ³ + (3ξcᵢ²)*Ωbᵢ + (2 - 3ξcᵢ)*Ωbᵢ² + Ωbᵢ³
Ωcᵢ = 1 - 3ξcᵢ

γ = ∑cᵢxᵢ/∑bᵢxᵢ
δ = 1 + 6γ + γ²
ϵ = 1 + γ

Δ₁ = -(ϵ + √δ)/2
Δ₂ = -(ϵ - √δ)/2
```

## Model Construction Examples
```julia
# Using the default database
model = YFR("water") #single input
model = YFR(["water","ethanol"]) #multiple components
model = YFR(["water","ethanol"], idealmodel = ReidIdeal) #modifying ideal model
model = YFR(["water","ethanol"],alpha = SoaveAlpha) #modifying alpha function
model = YFR(["water","ethanol"],translation = RackettTranslation) #modifying translation
model = YFR(["water","ethanol"],mixing = KayRule) #using another mixing rule
model = YFR(["water","ethanol"],mixing = WSRule, activity = NRTL) #using advanced EoS+gᴱ mixing rule

# Passing a prebuilt model

my_alpha = SoaveAlpha(["ethane","butane"],userlocations = Dict(:acentricfactor => [0.1,0.2]))
model = YFR(["ethane","butane"],alpha = my_alpha)

# User-provided parameters, passing files or folders

model = YFR(["neon","hydrogen"]; userlocations = ["path/to/my/db","cubic/my_k_values.csv"])

# User-provided parameters, passing parameters directly

model = YFR(["neon","hydrogen"];
        userlocations = (;Tc = [44.492,33.19],
                        Pc = [2679000, 1296400],
                        Mw = [20.17, 2.],
                        acentricfactor = [-0.03,-0.21]
                        k = [0. 0.18; 0.18 0.], #k,l can be ommited in single-component models.
                        l = [0. 0.01; 0.01 0.])
                    )
```

## Notes

The original paper also uses correlations for Ωaᵢ and Ωbᵢ. We opt to use the full solution instead to make the EoS consistent at the critical point.

## References

1. Yang, X., Frotscher, O., & Richter, M. (2025). Symbolic-regression aided development of a new cubic equation of state for improved liquid phase density calculation at pressures up to 100 MPa. International Journal of Thermophysics, 46(2). doi:10.1007/s10765-024-03490-5

"""
YFR

function YFR(components;
    idealmodel = BasicIdeal,
    alpha = YFRAlpha,
    mixing = YFR1fRule,
    activity = nothing,
    translation = NoTranslation,
    userlocations = String[],
    ideal_userlocations = String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    reference_state = nothing,
    verbose = false)

    formatted_components = format_components(components)

    params = getparams(formatted_components, ["properties/critical.csv", "properties/molarmass.csv"];
                        userlocations = userlocations,
                        verbose = verbose,
                        ignore_missing_singleparams = __ignored_crit_params(alpha))

    model = CubicModel(YFR,params,formatted_components;
                        idealmodel,alpha,mixing,activity,translation,
                        userlocations,ideal_userlocations,alpha_userlocations,activity_userlocations,mixing_userlocations,translation_userlocations,
                        reference_state, verbose)
    
    k = get(params,"k",nothing)
    l = get(params,"l",nothing)
    recombine_cubic!(model,k,l)
    set_reference_state!(model,reference_state;verbose)
    return model
end


#premixing: done at eos evaluation time
function ab_premixing(model::YFRModel,mixing::MixingRule,k,l)
    _Tc = model.params.Tc
    _pc = model.params.Pc
    a = model.params.a
    b = model.params.b
    for i in 1:length(model)
        Tci,pci = _Tc[i],_pc[i]
        a[i] = R̄^2*Tci^2/pci
        b[i] = R̄*Tci/pci
    end
    epsilon_LorentzBerthelot!(a,k)
    sigma_LorentzBerthelot!(b,l)
    return a,b
end

function c_premixing(model::YFRModel)
    _Tc = model.params.Tc
    _pc = model.params.Pc
    c = model.params.c
    c .= R̄ .*_Tc ./ _pc
    return c
end

function cubic_ΔT(model::YFRModel,T,z)
    b = diagvalues(model.params.b.values)
    c = diagvalues(model.params.c.values)
    b̄ = zero(Base.promote_eltype(b,z))
    c̄ = zero(Base.promote_eltype(c,z))
    Pc = model.params.Pc.values
    Vc = model.params.Vc.values
    Tc = model.params.Tc.values
    for i in 1:length(model)
        Zc = Vc[i]/b[i]
        Tr = T/Tc[i]
        Ωbi,Ωci = YFR_Ωbc(Zc,Tr)
        b̄ += z[i]*b[i]*Ωbi
        c̄ += z[i]*c[i]*Ωci
    end
    γ = c̄/b̄
    #=
    1 + 6g + g2 = 0
    γ > -3 +- 2*sqrt(2) = -0.1715728752538097 for this to work
    on YFR, undecane at 333K fails at this.
    =#
    δ = sqrt(evalpoly(γ,(1,6,1)))
    ϵ = 1 + γ
    return (-0.5*(ϵ + δ), -0.5*(ϵ - δ))
end

#technically we do dX/dTr = 0 instead of dX/dT = 0, but it is equivalent for eos evaluation purposes
#this is not correct if you are actually fitting the parameters.
function YFR_ξc(Zc,_Tr)
    Tr = primalval(_Tr) 
    return 0.144894*exp(-Tr^4) − 0.129835*exp(-Tr^3) + 0.957454*Zc + 0.036884
end

function YFR_Ωc(Zc,Tr)
    ξc = YFR_ξc(Zc,Tr)
    return 1 - 3ξc
end

#we opt for the exact solution of Ωb and Ωa instead of solving for the root
function YFR_Ωb(Zc,Tr,ξc = YFR_ξc(Zc,Tr))
    Ωb = PatelTeja_Ωb(ξc)
    #Ωb = 0.048371*exp(-Tr^4) −0.043334*exp(-Tr^3) + 0.319103*Zc −0.012341
end

function YFR_Ωab(Zc,Tr)
    ξc = YFR_ξc(Zc,Tr)
    Ωa,Ωb = PatelTeja_Ωab(ξc)
    #Ωa = −0.174696*exp(-Tr^4) + 0.156625*exp(-Tr^3) − 1.158565*Zc + 0.784751
    return Ωa,Ωb
end

function YFR_Ωbc(Zc,Tr)
    ξc = YFR_ξc(Zc,Tr)
    Ωb = YFR_Ωb(Zc,Tr,ξc)
    return Ωb,1 - 3ξc
end

function crit_pure(model::YFRModel)
    single_component_check(crit_pure,model)
    Tc = model.params.Tc.values[1]
    Pc = model.params.Pc.values[1]
    b = cubic_lb_volume(model,Tc,SA[1.0])
    Δ1,Δ2 = cubic_ΔT(model,Tc,SA[1.0])
    RT = Rgas(model)*Tc
    RTp = RT/Pc
    Vc0 = real((RTp + (Δ1 + Δ2 + 1)*b)/3)
    c = translation(model,Vc0,Tc,SA[1.0])
    Vc = Vc0 - c[1]
    return Tc,Pc,Vc
end
