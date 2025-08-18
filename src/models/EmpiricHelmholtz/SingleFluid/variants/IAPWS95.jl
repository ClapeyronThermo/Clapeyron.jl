"""
    IAPWS95 <: EmpiricHelmholtzModel
    IAPWS95()

## Input parameters

None

## Description

IAPWS95 (International Association for the Properties of Water and Steam) Pure water Model, 2018 update.

```
δ = ρ/ρc
τ = T/Tc
a⁰(δ,τ) = log(δ) + n⁰₁ + n⁰₂τ + n⁰₃log(τ) + ∑n⁰ᵢ(1-exp(-γ⁰ᵢτ)), i ∈ 4:8
aʳ(δ,τ)  = aʳ₁+ aʳ₂ + aʳ₃ + aʳ₄
aʳ₁(δ,τ)  = ∑nᵢδ^(dᵢ)τ^(tᵢ), i ∈ 1:7
aʳ₂(δ,τ)  = ∑nᵢexp(-δ^cᵢ)δ^(dᵢ)τ^(tᵢ), i ∈ 8:51
aʳ₃(δ,τ)  = ∑nᵢexp(-αᵢ(δ - εᵢ)^2 - βᵢ(τ - γᵢ)^2)δ^(dᵢ)τ^(tᵢ), i ∈ 52:54
aʳ₄(δ,τ) = ∑nᵢδΨΔ^(bᵢ), i ∈ 55:56
Δ = θ^2 + Bᵢ[(δ - 1)^2]^aᵢ
θ = (1 - τ) + Aᵢ[(δ - 1)^2]^(1/2βᵢ)
Ψ = exp(-Cᵢ(δ - 1)^2 - Dᵢ(τ - 1)^2)
```
Parameters `n⁰`,`γ⁰`,`n`,`t`,`d`,`c`,`α`,`β`,`γ`,`ε`,`A`,`B`,`C`,`D` where obtained via fitting.

## References

1. Wagner, W., & Pruß, A. (2002). The IAPWS formulation 1995 for the thermodynamic properties of ordinary water substance for general and scientific use. Journal of physical and chemical reference data, 31(2), 387–535. [doi:10.1063/1.1461829](https://doi.org/10.1063/1.1461829)
2. IAPWS R6-95 (2018). Revised Release on the IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary Water Substance for General and Scientific Use

"""
function IAPWS95()

    components = ["water"]

    Mw = 18.015268 #g·mol-1
    T_c = 647.096  #K
    P_c = 2.2064e7 #Pa
    rho_c= 17873.72799560906 # mol·m-3, calculated from rho_c = 322 kg·m-3
    Tr = T_c
    rhor = rho_c
    lb_volume = 1.4393788065379039e-5
    Ttp = 273.16 #K
    ptp = 611.6548008968684
    rhov_tp  = 0.2694716052752858
    rhol_tp = 55496.95513999978
    Rgas = 8.3143713575874
    acentric_factor = 0.3442920843

    properties = SingleFluidProperties(Mw,T_c,rho_c,lb_volume,T_c,P_c,rho_c,Ttp,ptp,rhov_tp,rhol_tp,acentric_factor,Rgas)

    a₁ = -8.3204464837497
    a₂ = 6.6832105275932
    c0 = 4.00632 - 1
    n_pe = [0.012436, 0.97315, 1.2795, 0.96956, 0.24873]
    t_pe = -[1.28728967, 3.53734222, 7.74073708, 9.24437796,27.5075105]
    ideal = SingleFluidIdealParam(a₁,a₂,c0,n_pe,t_pe)

    n = [0.12533547935523e-1, 0.78957634722828e1, -0.87803203303561e1,
    0.31802509345418, -0.26145533859358, -0.78199751687981e-2,
    0.88089493102134e-2,
    -0.66856572307965, 0.20433810950965, -0.66212605039687e-4,
    -0.19232721156002, -0.25709043003438, 0.16074868486251,
    -0.040092828925807, 0.39343422603254e-6, -0.75941377088144e-5,
    0.56250979351888e-3, -0.15608652257135e-4, 0.11537996422951e-8,
    .36582165144204e-6, -.13251180074668e-11, -.62639586912454e-9,
    -0.10793600908932, 0.17611491008752e-1, 0.22132295167546,
    -0.40247669763528, 0.58083399985759, 0.49969146990806e-2,
    -0.31358700712549e-1, -0.74315929710341, 0.47807329915480,
    0.20527940895948e-1, -0.13636435110343, 0.14180634400617e-1,
    0.83326504880713e-2, -0.29052336009585e-1, 0.38615085574206e-1,
    -0.20393486513704e-1, -0.16554050063734e-2, .19955571979541e-2,
    0.15870308324157e-3, -0.16388568342530e-4, 0.43613615723811e-1,
    0.34994005463765e-1, -0.76788197844621e-1, 0.22446277332006e-1,
    -0.62689710414685e-4, -0.55711118565645e-9, -0.19905718354408,
    0.31777497330738, -0.11841182425981,
    -0.31306260323435e2, 0.31546140237781e2, -0.25213154341695e4]

    t = [-0.5, 0.875, 1, 0.5, 0.75, 0.375, 1,
    4, 6, 12, 1, 5, 4, 2, 13, 9, 3, 4, 11, 4, 13, 1, 7, 1, 9, 10, 10, 3, 7,
    10, 10, 6, 10, 10, 1, 2, 3, 4, 8, 6, 9, 8, 16, 22, 23,23, 10, 50, 44, 46, 50,
    0,1,4]

    d = [1, 1, 1, 2, 2, 3, 4,
    1, 1, 1, 2, 2, 3, 4, 4, 5, 7, 9, 10, 11, 13, 15, 1, 2, 2, 2, 3,
    4, 4, 4, 5, 6, 6, 7, 9, 9, 9, 9, 9, 10, 10, 12, 3, 4, 4, 5, 14, 3, 6, 6, 6,
    3,3,3]

    l = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 6, 6,6, 6]

    η = [20.0,20.0,20.0]
    β = [150.0,150.0,250.0]
    γ = [1.21,1.21,1.25]
    ε = [1.,1.,1.]

    polexpgauss = PolExpGaussTerm(n,t,d,l,ones(length(l)),η,β,γ,ε)
    
    NA_A = [0.32,0.32]
    NA_B = [0.2,0.2]
    NA_C = [28,32]
    NA_D = [700,800]
    NA_a = [3.5,3.5]
    NA_b = [0.85,0.95]
    NA_beta = [0.3,0.3]
    NA_n = [-0.14874640856724,0.31806110878444]
    na = NonAnalyticTerm(NA_A,NA_B,NA_C,NA_D,NA_a,NA_b,NA_beta,NA_n)
    
    residual = SingleFluidResidualParam(;polexpgauss,na)

    ancillary_gas = GenericAncEvaluator([0.9791749335365787, -2.6190679042770215, -3.9166443712365235, -20.313306821636637, 16.497589490043744, -125.36580458432083],[0.21, 0.262, 0.701, 3.909, 4.076, 17.459],T_c,rho_c,:exp,true) |> PolExpVapour
    ancillary_liquid = GenericAncEvaluator([0.8157021355019343, 2.0434712177006693, -78.58278372496308, 1026.4273940070307, -2290.5642779377695, 8420.141408210317],[0.276, 0.455, 7.127, 9.846, 11.707, 17.805],T_c,rho_c,:noexp,false) |> PolExpLiquid
    ancillary_pressure = GenericAncEvaluator([-9.75639641045262, 3.3357600887120102, -1.10029278432831, 0.02037617155190105, -2.6668589845604367, 6.676721087238668],[1.018, 1.206, 2.327, 5.753, 4.215, 14.951],T_c,P_c,:exp,true) |> PolExpSat
    ancillaries = CompositeModel(components,gas = ancillary_gas,liquid = ancillary_liquid,saturation = ancillary_pressure)

    references = ["IAPWS R6-95(2018)"]

    return SingleFluid(components,properties,ancillaries,ideal,residual,references)
end

IAPWS95(components;idealmodel = nothing,userlocations = nothing,verbose = false,reference_state = nothing) = IAPWS95()


"""
    IAPWS95Ideal <: IdealModel
    IAPWS95Ideal(components;
    userlocations::Array{String,1}=String[],
    verbose = false)

    IAPWS95Ideal()

## Input parameters

None

## Description

IAPWS95 ideal helmholtz model for use in other models. Only valid for water. Check [`IAPWS95`](@ref) for more information.

## References

1. Wagner, W., & Pruß, A. (2002). The IAPWS formulation 1995 for the thermodynamic properties of ordinary water substance for general and scientific use. Journal of physical and chemical reference data, 31(2), 387–535. [doi:10.1063/1.1461829](https://doi.org/10.1063/1.1461829)
2. IAPWS R6-95 (2018). Revised Release on the IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary Water Substance for General and Scientific Use

"""
function IAPWS95Ideal()
    return idealmodel(IAPWS95())
end

export IAPWS95