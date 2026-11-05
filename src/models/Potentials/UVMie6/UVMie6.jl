module UVMie6Utils

# Table I: c0..c7
const UVMIE_C = (26.454, 1.8045, 1.7997, 161.96, 11.605, 12.0, 0.4, 2.0)

# Table II: k0i, k1i, k2i, k3i — Eqs. (20), (22)
const UVMIE_K = (
    (-3.9806,  79.565,  0.5489,  5.3632),
    ( 1.42455, 57.292,  0.0,     1.0031),
    (-39.755, -81.213,  0.6987, 30.156 ),
    (-23.692, -85.006,  0.7762, 12.798 ),
)

# Table III: p_i, l_i, m_i — Eqs. (23), (24)
const UVMIE_PLM = (
    ( 0.808444,   4.9485,  0.118532),
    (-0.095409, -21.300,   0.078556),
    ( 0.475257,   7.0000, -0.550394),
    (-2.8328,     3.2162,  0.0091630),
)

# Table VII: c1..c9 used in Eq. (A8)
const UVMIE_C_MAT = (
    ( 1.92840364363978,   0.44316589626508,   0.520120816141764),
    ( 0.182526759234415,  0.01103199896599296, -7.97813995328348e-5),
    ( 0.01298851560872428, 6.41039871789327e-3, 1.85866741090323e-5),
)

# Table VIII: A_ij, B_ij used in Eqs. (A12), (A13)
const UVMIE_A = (
    ( 2.64043218,    -1.2184421,    -22.90786387,  0.96433414),
    (-16.75643936,  630.83929772,     3.08711814, -166.57701616),
    ( 19.53170162,  -88.87955657,   -76.51387192,  443.68942745),
    ( -3.77740877,   83.04694547,    21.62502721, -304.86431760),
)
const UVMIE_B = (
    ( 2.19821588, -20.45005484),
    (-13.47050687, 56.65701375),
    ( 12.90119266, -42.71680606),
)

# Table IX: c_i1..c_i6 used in Eq. (A19), 6 rows
const UVMIE_IU = (
    ( -0.2622378162,   0.6585817423,   5.5318022309,   0.6902354794,  -3.6825190645,   -1.7263213318),
    ( -0.1899241690,  -0.5555205158,   9.1361398949,   0.7966155658,  -6.1413017045,    4.9553415149),
    (  0.1169786415,  -0.2216804790,  -2.0470861617,  -0.3742261343,   0.9568416381,   10.1401796764),
    (  0.5852642702,   2.0795520346,  19.0711829725,  -2.3403594600,   2.5833371420,  432.3858674425),
    ( -0.6084232211,  -7.2376034572,  19.0412933614,   3.2388986513,  75.4442555789, -588.3837110653),
    (  0.0512327656,   6.6667943569,  47.1109947616,  -0.5011125797, -34.8918383146,  189.5498636006),
)

# Table X: c_i0..c_i3 used in Eq. (A21)
const UVMIE_RSAP = (
    (-0.06355098,  6.206829830, -37.4582955,  40.7284977),
    ( 1.51905341, 13.1498964,   85.3505867, 374.1906360),
    ( 0.69345622,  9.4599462,  -53.2898422, 315.8199084),
    ( 0.00749260,  0.5461712,    7.9795626, -119.6126395),
)
# ------------------------------------------------------------------
# Reference fluid a₀
# ------------------------------------------------------------------
# Eq. (A9)
function η_A(η, data)
    A = UVMIE_A
    νinv = data.νinv
    τ = data.τ
    A2 = τ*((A[2][1] + A[2][2]*νinv) + (A[2][3] + A[2][4]*νinv)*τ)
    A3 = τ*((A[3][1] + A[3][2]*νinv) + (A[3][3] + A[3][4]*νinv)*τ)
    A4 = τ*((A[4][1] + A[4][2]*νinv) + (A[4][3] + A[4][4]*νinv)*τ)
    η*evalpoly(η, (data.A1 + 1, A2, A3, A4))
end

# Eq. (A10)
function η_B(η, data)
    B = UVMIE_B
    τ = data.τ
    B2 = τ*(B[2][1] + B[2][2]*τ)
    B3 = τ*(B[3][1] + B[3][2]*τ)
    η*evalpoly(η, (data.B1 + 1, B2, B3))
end

# Eq. (A19)
C_i(c, τ, νinv) = evalpoly(νinv,(c[1],c[2],c[3])) + evalpoly(νinv,(c[4],c[5],c[6]))*τ

end #utils module

abstract type UVMieModel <: EoSModel end

struct UVMie{T} <: UVMieModel
    nu::T
end

function __data(model::UVMieModel, T̄)
    ν = model.nu
    # Eq. (2)
    pref = ν/(ν-6) * (ν/6)^(6/(ν-6))
    # Eq. (3)
    rm   = (ν/6)^(1/(ν-6))
    # Eq. (A4)
    c    = (ν/6)^(-1/(12-2ν)) - 1
    # Eq. (A3)
    sqrtT̄ = sqrt(T̄)
    d    = rm * (1/(1 + c*sqrtT̄))^(2/ν)

    # Eq. (A7) with coefficients from Eq. (A8)
    M   = UVMie6Utils.UVMIE_C_MAT
    v7  = ν - 7
    v72 = v7*v7
    v73 = v72*v7
    C1q = M[1][1]*v7 + M[1][2]*v72 + M[1][3]*v73
    C2q = M[2][1]*v7 + M[2][2]*v72 + M[2][3]*v73
    C3q = M[3][1]*v7 + M[3][2]*v72 + M[3][3]*v73
    bracket = 1 + sqrt(2π*ν/6)*sqrtT̄ + C1q*T̄ + C2q*T̄*sqrtT̄ + C3q*T̄*T̄
    q = rm * bracket^(-1/(2ν))

    # Eq. (A11)
    τ  = rm - d
    d3  = d*d*d
    q3  = q*q*q
    rm3 = rm*rm*rm

    # Eq. (A18)
    αmf = pref*(rm3/3 - rm^(3-ν)/(3-ν))

    νinv  = 1/ν

    # C1 from Eq. (A19) — reused in _I_u [Eq. (A17)] and _ΔB31 [Eq. (A26)]
    C1  = UVMie6Utils.C_i(UVMie6Utils.UVMIE_IU[1], τ, νinv)

    # A1 [Eq. (A12)] — reused in _η_A [Eq. (A9)] and _B30_tilde [Eq. (A23)]
    a1r = UVMie6Utils.UVMIE_A[1]
    A1  = τ*((a1r[1] + a1r[2]*νinv) + (a1r[3] + a1r[4]*νinv))
    # B1 [Eq. (A13)] — reused in _η_B [Eq. (A10)] and _B30_tilde [Eq. (A23)]
    b1r = UVMie6Utils.UVMIE_B[1]
    B1  = τ*(b1r[1] + b1r[2]*τ)

    return (;pref, d, q, τ, d3, q3, rm3, αmf,
              νinv, C1, A1, B1)
end

# Eqs. (A1), (A2)
function a0(model::UVMieModel, ρs, T̄, data = __data(model, T̄))
    η  = π/6 * ρs * data.d3
    η2 = η*η
    a_hs = (4η - 3η2)/(1-η)^2

    # Eq. (A6)
    ηA = UVMie6Utils.η_A(η, data)
    ηB = UVMie6Utils.η_B(η, data)
    termA = (1 - ηA/2)/(1-ηA)^3 * (data.rm3 - data.q3)
    termB = (1 - ηB/2)/(1-ηB)^3 * (data.rm3 - data.d3)
    Δa0 = -2/3 * π * ρs * (termA - termB)
    return a_hs + Δa0
end

# ------------------------------------------------------------------
# First-order perturbation Δa₁ᵘ
# ------------------------------------------------------------------
# Eq. (A17)
function I_u(model::UVMieModel, ρs, T, data = __data(model,T))
    τ = data.τ
    C2 = UVMie6Utils.C_i(UVMie6Utils.UVMIE_IU[2], τ, data.νinv)
    C3 = UVMie6Utils.C_i(UVMie6Utils.UVMIE_IU[3], τ, data.νinv)
    C4 = UVMie6Utils.C_i(UVMie6Utils.UVMIE_IU[4], τ, data.νinv)
    C5 = UVMie6Utils.C_i(UVMie6Utils.UVMIE_IU[5], τ, data.νinv)
    C6 = UVMie6Utils.C_i(UVMie6Utils.UVMIE_IU[6], τ, data.νinv)

    ρs2 = ρs*ρs
    ρs3 = ρs2*ρs
    num = data.C1*ρs + C2*ρs2 + C3*ρs3
    den = 1 + C4*ρs + C5*ρs2 + C6*ρs3
    (data.q3 - data.rm3)/3 - data.αmf + num/den
end


# ------------------------------------------------------------------
# Virial coefficients
# ------------------------------------------------------------------

# Eqs. (A20), (A21)
function ΔB2(model::UVMieModel, T̄, data = __data(model, T̄))
    βs   = 1/T̄
    c    = UVMie6Utils.UVMIE_RSAP
    νinv = data.νinv
    ν2inv = νinv*νinv
    ν3inv = ν2inv*νinv
    c1 = c[1][1] + c[1][2]*νinv + c[1][3]*ν2inv + c[1][4]*ν3inv
    c2 = c[2][1] + c[2][2]*νinv + c[2][3]*ν2inv + c[2][4]*ν3inv
    c3 = c[3][1] + c[3][2]*νinv + c[3][3]*ν2inv + c[3][4]*ν3inv
    c4 = c[4][1] + c[4][2]*νinv + c[4][3]*ν2inv + c[4][4]*ν3inv
    B20 = 2/3*π*data.q3
    Y   = expm1(βs)
    s   = sum(βs^i/(factorial(i)*i) for i in 1:15)
    (B20 - 2/3*π*data.rm3 - c1)*Y - c2*s - c3*βs - c4*βs^2
end


function ΔB3(model::UVMieModel,T̄, data = __data(model, T̄))
    ν = model.nu

    #B3_total, Eqs. (19), (21)–(24)
    βs   = 1/T̄
    q3,d3,rm3 = data.q3,data.d3,data.rm3
    B30  = 10*(π/6*q3)^2
    νm12 = ν - 12
    νm6  = ν - 6
    a1, a2, a3, a4 = ntuple(4) do i
        k0, k1, k2, k3 = UVMie6Utils.UVMIE_K[i]
        p, l, m        = UVMie6Utils.UVMIE_PLM[i]
        k1t = k1*(1 + p*νm12/(ν-l))
        k3t = k3*(1 + m*νm12/νm6)
        k0 + k1t/(k2 + T̄)^k3t
    end

    B3_total = evalpoly(βs, (B30, a1, a2, a3, a4))

    # Eq. (A23)
    A1,B1 = data.A1,data.B1
    B30_tilde = 10*(π/6*d3)^2 - 5/9*π^2*d3 * ((1 + A1)*(rm3 - q3) - (1 + B1)*(rm3 - d3))

    return B3_total - B30_tilde
end

function a_virial(model::UVMieModel,ρs,T,data = __data(model,T))
    ΔB₂  = ΔB2(model, T̄, data)
    ΔB₂₁ = 2π*(1/T̄)*((data.q3 - data.rm3)/3 - data.αmf) # Eq. (A24)
    ΔB₃  = ΔB3(model, T̄, data)
    ΔB₃₁ = π*(1/T̄) * data.pref * data.C1 # Eq. (A26)

    # Eq. (13)
    ν = model.nu
    c0, c1 = 26.454, 1.8045
    α = exp(-(c0/ν)*(1 + c1/T̄^2)*ρs)

    return (ΔB₂ - ΔB₂₁)*ρs + 0.5*α*(ΔB₃ - ΔB₃₁)*ρs^2
end

# ------------------------------------------------------------------
# u-fraction and damping
# ------------------------------------------------------------------

# Eq. (15)
function u_fraction(model::UVMieModel, ρs, T̄)
    ν = model.nu
    c2, c3, c4 = 1.7997, 161.96, 11.605

    # Eq. (16)
    c5, c6, c7 = 12.0, 0.4, 2.0
    ω = exp(-c5*T̄*(ρs-c6)^2) * (1/tanh(c7*T̄) - 1)^2

    1 - exp(-(c2 + c3/(c4+ν))*ρs^2 - ω)
end

# ------------------------------------------------------------------
# Eq. (14)
# ------------------------------------------------------------------
function a_potential(model::UVMieModel, V, T̄)
    ρs = 1/V
    data = __data(model, T̄)
    a0   = a0(model, ρs, T̄, data)
    Δa1  = 2π*ρs*(1/T̄)*I_u(model, ρs, T̄, data) # Eq. (A15)
    ϕu   = u_fraction(model, ρs, T̄)
    av = a_virial(model,ρs,T̄,data)
    a0 + Δa1 + (1 - ϕu)*av
end

# Eq. (17)
function Clapeyron.lb_volume(model::UVMieModel, T, z)
    ρ_max = 1.1 + 0.12 * sqrt(T)
    sum(z) / ρ_max
end

Clapeyron.x0_volume_liquid(model::UVMieModel, p, T, z) = 1.25 * lb_volume(model, T, z)

function second_virial_coefficient_impl(model::UVMieModel, T, z)
    data = __data(model, T)
    B20 = 2/3 * π * data.q3
    B2  = B20 + ΔB2(model, T, data)
    B2 * sum(z)
end