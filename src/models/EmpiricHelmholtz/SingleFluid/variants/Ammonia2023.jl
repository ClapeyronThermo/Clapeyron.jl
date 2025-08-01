"""
    Ammonia2023 <: EmpiricHelmholtzModel
    Ammonia2023()

## Input parameters

None

## Description
Ammonia Reference Equation of State (2023)
```
δ = ρ/ρc
τ = T/Tc
a⁰(δ,τ) = log(δ) + n⁰₁ + n⁰₂τ + n⁰₃log(τ) + ∑n⁰ᵢ(1-exp(-γ⁰ᵢτ)), i ∈ 4:7
aʳ(δ,τ)  = aʳ₁+ aʳ₂ + aʳ₃
aʳ₁(δ,τ)  = ∑nᵢδ^(dᵢ)τ^(tᵢ), i ∈ 1:5
aʳ₂(δ,τ)  = ∑nᵢexp(-δ^cᵢ)δ^(dᵢ)τ^(tᵢ), i ∈ 6:8
aʳ₃(δ,τ)  = ∑nᵢexp(-ηᵢ(δ - εᵢ)^2 - βᵢ(τ - γᵢ)^2)δ^(dᵢ)τ^(tᵢ), i ∈ 9:18
aʳ₃(δ,τ)  = ∑nᵢexp(-ηᵢ(δ - εᵢ)^2 - 1/(βᵢ*(τ -γᵢ)^2 + bᵢ))δ^(dᵢ)τ^(tᵢ), i ∈ 19:20

```
Parameters  `n⁰`,`γ⁰`,`n`,`t`,`d`,`c`,`η`,`β`,`γ`,`ε` where obtained via fitting.

## References
1. Gao, K., Wu, J., Bell, I. H., Harvey, A. H., & Lemmon, E. W. (2023). A reference equation of state with an associating term for the thermodynamic properties of ammonia. Journal of Physical and Chemical Reference Data, 52(1), 013102. [doi:10.1063/5.0128269](https://doi.org/10.1063/5.0128269)
"""
Ammonia2023

function Ammonia2023()

    components = ["ammonia"]

    Mw = 17.03052 #g·mol-1
    T_c = 405.56    #K
    P_c = 11.3634e6 #Pa
    rho_c= 13696.0 # mol·m-3
    lb_volume = 1/53130
    Ttp = 195.49 #K
    ptp = 6.05339e3
    rhov_tp  = 0.003740e3
    rhol_tp = 43.091e3
    Rgas = 8.314462618
    acentric_factor = NaN

    properties = SingleFluidProperties(Mw,T_c,rho_c,lb_volume,T_c,P_c,rho_c,Ttp,ptp,rhov_tp,rhol_tp,acentric_factor,Rgas)

    a₁ = -6.59406093943886
    a₂ = 5.601011519879
    u = -[4.0585856593352405, 9.776605187888352, 17.829667620080876]
    v = [2.224,3.148,0.9579]
    c0 = 4.0 - 1
    ideal = SingleFluidIdealParam(a₁,a₂,c0,v,u)

    n = [0.006132232,1.7395866,-2.2261792,-0.30127553,0.08967023,-0.076387037,-0.84063963,-0.27026327,6.212578,-5.7844357,2.4817542,-2.3739168,0.01493697,-3.7749264,0.0006254348,
    -1.7359e-05,-0.13462033,0.07749072839]
    t = [1.0,0.382,1.0,1.0,0.677,2.915,3.51,1.063,0.655,1.3,3.1,1.4395,1.623,0.643,1.13,4.5,1.0,4.0]
    d = [4,1,1,2,3,3,2,3,1,1,1,2,2,1,3,3,1,1]
    l = [2,2,1]
    g = [1,1,1]
    η = [0.42776,0.6424,0.8175,0.7995,0.91,0.3574,1.21,4.14,22.56,22.68]
    β = [1.708,1.4865,2.0915,2.43,0.488,1.1,0.85,1.14,945.64,993.85]
    γ = [1.036,1.2777,1.083,1.2906,0.928,0.934,0.919,1.852,1.05897,1.05277]
    ε = [-0.0726,-0.1274,0.7527,0.57,2.2,-0.243,2.96,3.02,0.9574,0.9576]

    gaob_n = [-1.6909858,0.93739074]
    gaob_t = [4.3315,4.015]
    gaob_d = [1.,1.]
    gaob_eta = [-2.8452,-2.8342]
    gaob_beta = [0.3696,0.2962]
    gaob_gamma = [1.108,1.313]
    gaob_epsilon = [0.4478,0.44689]
    gaob_b = [1.244,0.6826]

    gao_b = GaoBTerm(gaob_n,gaob_t,gaob_d,gaob_eta,gaob_beta,gaob_gamma,gaob_epsilon,gaob_b)
    polexpgauss = PolExpGaussTerm(n,t,d,l,g,η,β,γ,ε)
    residual = SingleFluidResidualParam(;gao_b,polexpgauss)
    ancillary_gas = GenericAncEvaluator([-0.089966,-3.8722,-8.1183,-25.293,-54.279,-400.83],[0.112,0.473,1.5,3.875,8.0,20.0],T_c,rho_c,:exp,false) |> PolExpVapour
    ancillary_liquid = GenericAncEvaluator([0.051236,3.7925,-3.5929,4.6409,-1.9893,1.5978],[0.07,0.46,0.77,1.05,1.25,8.0],T_c,rho_c,:noexp,false) |> PolExpLiquid
    ancillary_pressure = GenericAncEvaluator([-7.3128,3.8888,-2.9908,-2.8636],[1.0,1.5,1.6,3.7],T_c,P_c,:exp,true) |> PolExpSat
    ancillaries = CompositeModel(components,gas = ancillary_gas,liquid = ancillary_liquid,saturation = ancillary_pressure)

    references = ["10.1063/5.0128269"]

    return SingleFluid(components,properties,ancillaries,ideal,residual,references)
end

export Ammonia2023
