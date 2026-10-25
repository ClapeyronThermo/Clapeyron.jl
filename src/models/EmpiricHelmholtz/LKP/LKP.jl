abstract type LKPModel <: EmpiricHelmholtzModel end

abstract type LKPReferenceModel <: EoSModel end 
#=
all instances need to define:
- lkp_Zc(model::LKPReferenceModel,ω)
- a_res_lkp(model::LKPReferenceModel,V,T,(T̄,V̄,ω̄))
=#
is_splittable(::LKPReferenceModel) = false

struct LKPMBWR <: LKPReferenceModel
    params_simple::NTuple{13,Float64}
    params_reference::NTuple{13,Float64}
    type::Symbol
end

show_info_lkp(io,model::LKPMBWR) = print(io,("MBWR ($(model.type))"))

lkp_Zc(::LKPMBWR,ω) = 0.2905 - 0.085*ω

function lkp0()
    #"b1", "b2", "b3", "b4", "c1", "c2", "c4", "c3", "d1", "d2", "𝛽", "𝛾", "ω"
    params_simple = (0.1181193, 0.265728, 0.15479, 0.030323, 0.0236744, 0.0186984, 0.0, 0.042724, 1.55428e-5, 6.23689e-5, 0.65392, 0.060167, 0.0)
    params_reference = (0.2026579, 0.331511, 0.027655, 0.203488, 0.0313385, 0.0503618, 0.016901, 0.041577, 4.8736e-5, 7.40336e-6, 1.226, 0.03754, 0.3978)
    return LKPMBWR(params_simple,params_reference,:original)
end

function a_res_lkp(model::LKPMBWR,V,T,data)
    T̄,V̄,ω̄ = data
    Zc = lkp_Zc(model,ω̄)
    δ = V̄/V
    τ = T̄/T
    αr_0    = reduced_a_res_lkp(model.params_simple,δ,τ,Zc)
    αr_ref  = reduced_a_res_lkp(model.params_reference,δ,τ,Zc)
    ω0,ωref = last(model.params_simple),last(model.params_reference)
    ωᵣ = (ω̄ - ω0)/(ωref - ω0)
    return (1 - ωᵣ)*αr_0 + ωᵣ*αr_ref
end

function reduced_a_res_lkp(params,δ,τ,Zc)
    b1,b2,b3,b4,c1,c2,c3,c4,d1,d2,β,γ = params
    B = evalpoly(τ,(b1,-b2,-b3,-b4))
    C = evalpoly(τ,(c1,-c2,zero(c2),c3))
    D = d1 + d2*τ
    c4τ = (c4/(2*γ))*τ^3
    δr = δ/Zc
    δr2 = δr*δr
    return δr*B + 0.5*C*δr2 + 0.2*(D*δr2*δr2*δr) - c4τ*(γ*δr2 + β + 1) * exp(-γ*δr2) + c4τ*(β + 1)
end

struct LKPParam <: EoSParam
    Mw::SingleParam{Float64}
    Tc::SingleParam{Float64}
    Pc::SingleParam{Float64}
    Vc_fit::SingleParam{Float64}
    k::PairParam{Float64}
    acentricfactor::SingleParam{Float64}
end

struct LKP{T,R} <: LKPModel
    components::Array{String,1}
    params::LKPParam
    idealmodel::T
    references::Array{String,1}
    reference_model::R
end

function show_info(io,model::LKP) 
    ref = model.reference_model
    println(io)
    print(io,"Reference model: ")
    show_info_lkp(io,ref)
end

function LKP(components;
    idealmodel = BasicIdeal,
    reference_model = lkp0(),
    userlocations = String[],
    ideal_userlocations = String[],
    reference_state = nothing,
    verbose = false)

    _components = format_components(components)
    params1 = getparams(_components, default_locations(LKP); userlocations = userlocations, verbose = verbose)
    params2 = transform_params_lkp(reference_model,params1,_components)
    packagedparams = build_eosparam(LKPParam,params2)
    init_idealmodel = init_model(idealmodel,_components,ideal_userlocations,verbose)
    references = default_references(LKP)
    model = LKP(_components,packagedparams,init_idealmodel, references, reference_model)
    set_reference_state!(model,reference_state;verbose)
    return model
end

"""
    LKP <: EmpiricHelmholtzModel
    LKP(components;
        idealmodel=BasicIdeal,
        verbose=false)

## Input parameters
- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Pressure `[Pa]`
- `Vc_fit`: Single Parameter (`Float64`) - Fitted Critical Volume `[m³·mol⁻¹]`
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `acentricfactor`: Single Parameter (`Float64`) - Acentric Factor (no units)
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Parameter (no units)

## Input models
- `idealmodel`: Ideal Model

## Description
Lee-Kesler-Plöker equation of state. Corresponding states using interpolation between a simple, spherical fluid (methane, `∅`)  and a reference fluid (n-octane, `ref`):
```
αᵣ = (1 - ωᵣ)*αᵣ(δr,τ,params(∅)) + ωᵣ*αᵣ(δr,τ,params(ref))
τ = Tr/T
δr = Vr/V/Zr
Zr = (0.2905 - 0.085*ω̄)
ωᵣ = (ω̄ - ω(∅))/(ω(ref) - ω(∅))
ω̄ = ∑xᵢωᵢ
Tr = ∑xᵢ*xⱼ*Tcᵢⱼ*Vcᵢⱼ^η * (1-kᵢⱼ)
Vr = ∑xᵢ*xⱼ*Tcᵢⱼ*Vcᵢⱼ
Tcᵢⱼ = √Tcᵢ*Tcⱼ
Vcᵢⱼ = 0.125*(∛Vcᵢ + ∛Vcⱼ)^3
η = 0.25
```

## Model Construction Examples
```julia
# Using the default database
model = LKP("water") #single input
model = LKP(["water","ethanol"]) #multiple components
model = LKP(["water","ethanol"], idealmodel = ReidIdeal) #modifying ideal model

# Passing a prebuilt model

my_idealmodel = MonomerIdeal(["neon","hydrogen"];userlocations = (;Mw = [20.17, 2.]))
model =  LKP(["neon","hydrogen"],idealmodel = my_idealmodel)

# User-provided parameters, passing files or folders
model = LKP(["neon","hydrogen"]; userlocations = ["path/to/my/db","lkp/my_k_values.csv"])

# User-provided parameters, passing parameters directly

model = LKP(["neon","hydrogen"];
        userlocations = (;Tc = [44.492,33.19],
                        Pc = [2679000, 1296400],
                        Vc = [4.25e-5, 6.43e-5],
                        Mw = [20.17, 2.],
                        acentricfactor = [-0.03,-0.21]
                        k = [0. 0.18; 0.18 0.]) #k,l can be ommited in single-component models.
                    )
```

## References
1. Plöcker, U., Knapp, H., & Prausnitz, J. (1978). Calculation of high-pressure vapor-liquid equilibria from a corresponding-states correlation with emphasis on asymmetric mixtures. Industrial & Engineering Chemistry Process Design and Development, 17(3), 324–332. [doi:10.1021/i260067a020](https://doi.org/10.1021/i260067a020)
"""
LKP

default_references(::Type{LKP}) = ["10.1021/i260067a020"]
default_locations(::Type{LKP}) = ["properties/critical.csv","properties/molarmass.csv","Empiric/LKP/LKP_unlike.csv"]
default_ignore_missing_singleparams(::Type{LKP}) = ["Vc"]


function get_k(model::LKPModel)   
    return copy(model.params.k.values)
end

function set_k!(model::LKPModel,k)
    check_arraysize(model,k)
    model.params.k.values .= k
end

function transform_params_lkp(refmodel::LKPReferenceModel,params,components)
    get!(params,"k") do
        PairParam("k",components)
    end

    Vc = get!(params,"Vc_fit") do
        SingleParam("Vc_fit",components)
    end

    Tc,Pc,ω = params["Tc"],params["Pc"],params["acentricfactor"]
    for i in eachindex(Vc.values)
        if Vc.ismissingvalues[i]
            Zci   = lkp_Zc(refmodel,ω[i])
            Vc[i] = Zci*Rgas()*Tc[i]/Pc[i]
        end
    end
    return params
end

function a_res(model::LKPModel,V,T,z)
    _data = @f(data)
    return a_res_lkp(model.reference_model,V/sum(z),T,_data)
end


data(model::LKPModel,V,T,z) = lkp_data(model,z)

function lkp_data(model::LKPModel,z)
    ω = model.params.acentricfactor.values
    T̄,V̄ = TV_scale(model,z)
    ω̄ = dot(ω,z)/sum(z)
    return T̄,V̄,ω̄
end

function V_scale(model::LKPModel,z)
    Vc = model.params.Vc_fit.values
    nc = length(model)
    V̄ = zero(Base.promote_eltype(model,z))
    ∑z = sum(z)
    for i in 1:nc
        Vci,zi = Vc[i],z[i]
        V̄ += zi*zi*Vci
        Vci_cbrt = cbrt(Vci) 
        for j in 1:(i-1)
            Vcj,zj = Vc[j],z[j]
            Vcij = 0.125*(Vci_cbrt + cbrt(Vcj))^3
            V̄ += 2*zi*zj*Vcij
        end
    end
    V̄ = V̄/∑z/∑z
    return V̄
end

function TV_scale(model::LKPModel,z)
    Vc = model.params.Vc_fit.values
    Tc = model.params.Tc.values
    nc = length(model)
    T̄ = zero(Base.promote_eltype(model,z))
    V̄η,V̄ = zero(T̄),zero(T̄)
    ∑z = sum(z)
    k = model.params.k.values
    for i in 1:nc
        Vci,Tci,zi = Vc[i],Tc[i],z[i]
        Vciη = sqrt(sqrt(Vci)) #Vci^0.25
        T̄ += zi*zi*Tci*Vciη
        V̄ += zi*zi*Vci
        Vci_cbrt = cbrt(Vci) 
        for j in 1:(i-1)
            Vcj,Tcj,zj = Vc[j],Tc[j],z[j]
            Vcij = 0.125*(Vci_cbrt + cbrt(Vcj))^3
            Vcijη = sqrt(sqrt(Vcij)) #Vci^0.25
            Tcij = sqrt(Tci*Tcj)*(1 - k[i,j])
            T̄ += 2*zi*zj*Tcij*Vcijη
            V̄ += 2*zi*zj*Vcij
        end
    end
    V̄ = V̄/∑z/∑z
    V̄η = V̄^0.25
    T̄ = T̄/∑z/∑z/V̄η
    return T̄,V̄
end

#=
function crit_pure(model::LKPModel)
    ω = only(model.params.acentricfactor.values)
    tc = only(model.params.Tc.values)
    vc = only(model.params.Vc.values)
    pc = (0.2905 - 0.085*ω)*Rgas(model)*tc/vc
    return tc,pc,vc
end =#

function x0_crit_pure(model::LKPModel,z)
    V̄ = V_scale(model,z)
    return one(V̄), log10(V̄)
end

has_fast_crit_pure(::LKPModel) = false

function x0_sat_pure(model::LKPModel,T)
    ω = only(model.params.acentricfactor.values)
    crit = crit_pure(model)
    _,vl,vv = x0_sat_pure_lk(model,T,crit,ω)
    return vl,vv
end

function x0_psat(model::LKPModel,T)
    ω = only(model.params.acentricfactor.values)
    crit = crit_pure(model)
    p,_,_ = x0_sat_pure_lk(model,T,crit,ω)
    return p
end

function lb_volume(model::LKPModel,T,z)
    V̄ = V_scale(model,z)
    return sum(z)*V̄*0.25
end

function T_scale(model::LKPModel,z)
    T̄,_ = TV_scale(model,z)
    return T̄
end

export LKP
