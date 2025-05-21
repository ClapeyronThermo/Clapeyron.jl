function TholLJ()
    components = ["Leonard-Jones Fluid (unscaled)"]
    Mw = 1. #g·mol-1
    T_c = 1.32    #K
    P_c = 0.13006 #Pa
    rho_c = 0.31 # mol·m-3
    lb_volume = 1/(π/6)
    Ttp = 0.661 #K
    ptp = NaN
    rhov_tp  = NaN
    rhol_tp = NaN
    Rgas = 1.0
    acentric_factor = NaN

    properties = SingleFluidProperties(Mw,T_c,rho_c,lb_volume,T_c,P_c,rho_c,Ttp,ptp,rhov_tp,rhol_tp,acentric_factor,Rgas)

    a₁ = 6.262265814
    a₂ = -1.515151515
    u = Float64[]
    v = Float64[]
    c0 = 2.5

    ideal = SingleFluidIdealParam(a₁,a₂,c0,v,u)

    n = [0.005208073, 2.186252, -2.161016, 1.4527, -2.041792,
        0.18695286, -0.090988445, -0.4974561, 0.10901431, -0.80055922,
        -0.568839, -0.6208625, -1.4667177, 1.891469, -0.1383701,
        -0.3869645, 0.1265702, 0.605781, 1.179189, -0.47732679,
        -9.9218575, -0.5747932, 0.003772923]
    t = [1.0, 0.32, 0.505, 0.672, 0.843, 0.898, 1.294, 2.59, 1.786, 2.77,
        1.786, 1.205, 2.83, 2.548, 4.65, 1.385, 1.46, 1.351, 0.66, 1.496,
        1.83, 1.616, 4.97]
    d = [4, 1, 1, 2, 2, 3, 5, 2, 2, 3, 1, 1, 1, 1, 2, 3, 3, 2, 1, 2, 3, 1, 1]
    l = [1, 2, 1, 2, 2, 1]
    η = [2.067, 1.522, 8.82, 1.722, 0.679, 1.883, 3.925, 2.461, 28.2, 0.753, 0.82]
    β = [0.625, 0.638, 3.91, 0.156, 0.157, 0.153, 1.16, 1.73, 383.0, 0.112, 0.119]
    γ = [0.71, 0.86, 1.94, 1.48, 1.49, 1.945, 3.02, 1.11, 1.17, 1.33, 0.24]
    ε = [0.2053, 0.409, 0.6, 1.203, 1.829, 1.397, 1.39, 0.539, 0.934, 2.369, 2.43]

    polexpgauss = PolExpGaussTerm(n,t,d,l,ones(length(l)),η,β,γ,ε)

    residual = SingleFluidResidualParam(;polexpgauss)

    anc_gas_fn = GenericAncEvaluator([-0.69655e+1,-0.10331e+3,-0.20325e+1,-0.44481e+2,-0.18463e+2,-0.26070e+3],[1.320 ,19.24,0.360,8.780,4.040,41.60],T_c,rho_c,:exp,false)
    anc_liquid_fn = GenericAncEvaluator([0.1362e+1,0.2093e+1,-0.2110e+1,0.3290e0,0.1410e+1],[0.313 ,0.940,1.630,17.,2.4],T_c,rho_c,:noexp,false)
    anc_ps_fn = GenericAncEvaluator([0.54000e+1,0.44704e01,-0.18530e+1,0.19890e0,-0.11250e+1],[1.,1.5,4.7,2.5,21.4],T_c,P_c,:exp,true)

    ancillary_gas = PolExpVapour(anc_gas_fn)
    ancillary_liquid = PolExpLiquid(anc_liquid_fn)
    ancillary_pressure = PolExpSat(anc_ps_fn)
    ancillaries = CompositeModel(components,gas = ancillary_gas,liquid = ancillary_liquid,saturation = ancillary_pressure)

    references = ["10.1063/1.4945000"]

    return SingleFluid(components,properties,ancillaries,ideal,residual,references)
end


"""
    TholLJ()

Lennard-Jones Reference equation of state. valid from 0.5 < T/Tc < 7 and pressures up to p/pc = 500.
```
τᵢ = 1.32/T
δᵢ = n/0.31V
a⁰ᵢ(δ,τ) = log(δᵢ) + 1.5log(τᵢ) - 1.515151515τᵢ + 6.262265814
a⁰(δ,τ,z) = ∑xᵢ(a⁰ᵢ + log(xᵢ))
aʳ(δ,τ)  = aʳ₁+ aʳ₂ + aʳ₃ + aʳ₄
aʳ₁(δ,τ)  = ∑nᵢδ^(dᵢ)τ^(tᵢ), i ∈ 1:6
aʳ₂(δ,τ)  = ∑nᵢexp(-δ^cᵢ)δ^(dᵢ)τ^(tᵢ), i ∈ 7:12
aʳ₃(δ,τ)  = ∑nᵢexp(-ηᵢ(δ - εᵢ)^2 - βᵢ(τ - γᵢ)^2)δ^(dᵢ)τ^(tᵢ), i ∈ 13:23
"""
TholLJ