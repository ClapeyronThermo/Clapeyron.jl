struct LKPParam <: EoSParam
    Mw::SingleParam{Float64}
    Tc::SingleParam{Float64}
    Pc::SingleParam{Float64}
    Vc::SingleParam{Float64}
    k::PairParam{Float64}
    acentricfactor::SingleParam{Float64}
end

abstract type LKPModel <: EmpiricHelmholtzModel end
@newmodel LKP LKPModel LKPParam false

"""
    LKP <: EmpiricHelmholtzModel
    LKP(components;
        idealmodel=BasicIdeal,
        verbose=false)

## Input parameters
- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Pressure `[Pa]`
- `Vc`: Single Parameter (`Float64`) (optional) - Critical Volume `[m³]`
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
Zr = Pr*Vr/(R*Tr)
Pr = (0.2905 - 0.085*ω̄)*R*Tr/Vr
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

function transform_params(::Type{LKP},params,components)
    k = get(params,"k",nothing)
    if k === nothing
        nc = length(components)
        params["k"] = PairParam("k",components)
    end
    _Vc = get(params,"Vc",nothing)
    if _Vc === nothing
        Vc = SingleParam("Vc",components)
        params["Vc"] = Vc
    else
        Vc = _Vc
    end
    Tc,Pc,ω = params["Tc"],params["Pc"],params["acentricfactor"]
    for i in 1:length(Vc)
        if Vc.ismissingvalues[i]
            Vc[i] = (0.2905 - 0.085*ω[i])*Rgas()*Tc[i]/Pc[i]
        end
    end
    return params
end

function get_k(model::LKPModel)   
    return copy(model.params.k.values)
end

function set_k!(model::LKPModel,k)
    check_arraysize(model,k)
    model.params.k.values .= k
end

function a_res(model::LKPModel,V,T,z = SA[1.0])
    Tr,Pr,Vr,ω̄ = @f(data)
    Zr = Pr*Vr/(Rgas(model)*Tr)
    δ = sum(z)*Vr/V
    τ = Tr/T
    δr = δ/Zr
    params_simple = lkp_params_simple(model)
    params_reference = lkp_params_reference(model)
    ω0,ωref = last(params_simple),last(params_reference)
    αr_0 = reduced_a_res_lkp(model,δ,τ,δr,params_simple)
    αr_ref = reduced_a_res_lkp(model,δ,τ,δr,params_reference)
    ωᵣ = (ω̄ - ω0)/(ωref - ω0)
    return (1 - ωᵣ)*αr_0 + ωᵣ*αr_ref
end
                                     #"b1", "b2", "b3", "b4", "c1", "c2", "c4", "c3", "d1", "d2", "𝛽", "𝛾", "ω"
lkp_params_simple(model::LKPModel) = (0.1181193, 0.265728, 0.15479, 0.030323, 0.0236744, 0.0186984, 0.0, 0.042724, 1.55428e-5, 6.23689e-5, 0.65392, 0.060167, 0.0)
lkp_params_reference(model::LKPModel) = (0.2026579, 0.331511, 0.027655, 0.203488, 0.0313385, 0.0503618, 0.016901, 0.041577, 4.8736e-5, 7.40336e-6, 1.226, 0.03754, 0.3978)

function data(model::LKPModel,V,T,z)
    ω = model.params.acentricfactor.values
    Vc = model.params.Vc.values
    Tc = model.params.Tc.values
    nc = length(model)
    T̄ = zero(1. + first(z))
    V̄η = zero(T̄)
    V̄ = zero(T̄)
    ∑z = sum(z)
    k = model.params.k.values
    for i in 1:nc
        Vci,Tci,zi = Vc[i],Tc[i],z[i]
        Vciη = sqrt(sqrt(Vci)) #Vci^0.25
        T̄ += zi*zi*Tci*Vciη
        V̄ += zi*zi*Vci
        for j in 1:(i-1)
            Vcj,Tcj,zj = Vc[j],Tc[j],z[j]
            Vcij = 0.125*(cbrt(Vci) + cbrt(Vcj))^3
            Vcijη = sqrt(sqrt(Vcij)) #Vci^0.25
            Tcij = sqrt(Tci*Tcj)*(1 - k[i,j])
            T̄ += 2*zi*zj*Tcij*Vcijη
            V̄ += 2*zi*zj*Vcij
        end
    end
    V̄ = V̄/∑z/∑z
    V̄η = V̄^0.25
    T̄ = T̄/∑z/∑z/V̄η
    ω̄ = dot(ω,z)/∑z
    #vci = (0.2905 - 0.085w)*RTci/Pci
    Pc = (0.2905 - 0.085*ω̄)*Rgas(model)*T̄/V̄
    return T̄,Pc,V̄,ω̄
end

function reduced_a_res_lkp(model::LKPModel,δ,τ,δr,params)
    b1,b2,b3,b4,c1,c2,c3,c4,d1,d2,β,γ,ω = params
    B = evalpoly(τ,(b1,-b2,-b3,-b4))
    C = evalpoly(τ,(c1,-c2,0.,c3))
    D = d1 + d2*τ
    c4τ = (c4/(2*γ))*τ^3
    return δr*B + 0.5*C*δr^2 + 0.2*(D*δr^5) - 
    c4τ*(γ*δr^2 + β + 1) * exp(-γ*δr^2) +
    c4τ*(β + 1)
end

function crit_pure(model::LKPModel)
    ω = only(model.params.acentricfactor.values)
    tc = only(model.params.Tc.values)
    vc = only(model.params.Vc.values)
    pc = (0.2905 - 0.085*ω)*Rgas(model)*tc/vc
    return tc,pc,vc
end

has_fast_crit_pure(::LKPModel) = true

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

function lb_volume(model::LKPModel,z)
    return lb_volume(model,0.0,z)
end

function lb_volume(model::LKPModel,T,z)
    V = 0.0
    Tc,Pc,Vc,ω̄ = @f(data)
    return sum(z)*Vc/4 #?
end

function T_scale(model::LKPModel,z)
    V,T = 0.0,0.0
    Tc,Pc,Vc,ω̄ = @f(data)
    return Tc
end

function p_scale(model::LKPModel,z)
    V,T = 0.0,0.0
    Tc,Pc,Vc,ω̄ = @f(data)
    return Pc
end

export LKP
