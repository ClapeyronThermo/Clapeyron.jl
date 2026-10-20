module SAFTVRMieUtils

using Clapeyron: SAFTVRMieconsts, SA, log1p

function aS1(λ,ζₓ)
    ζeff_ = ζeff(λ,ζₓ)
    return -1/(λ-3)*(1-ζeff_/2)/(1-ζeff_)^3
end

function ζeff(λ,ζₓ)
    A = SAFTVRMieconsts.A
    λ⁻¹ = one(λ)/λ
    Aλ⁻¹ = A * SA[one(λ); λ⁻¹; λ⁻¹*λ⁻¹; λ⁻¹*λ⁻¹*λ⁻¹]
    ζₓi = SA[ζₓ; ζₓ^2; ζₓ^3; ζₓ^4]
    return Aλ⁻¹[1]*ζₓi[1] + Aλ⁻¹[2]*ζₓi[2] + Aλ⁻¹[3]*ζₓi[3] + Aλ⁻¹[4]*ζₓi[4]
end

function B(λ, x_0, ζₓ)
    x_0_3λ = x_0^(3-λ)
    ζₓm13 = (1-ζₓ)^3
    I = (1-x_0_3λ)/(λ-3)
    J = (1-(λ-3)*x_0^(4-λ)+(λ-4)*x_0_3λ)/((λ-3)*(λ-4))
    return I*(1-ζₓ/2)/ζₓm13-9*J*ζₓ*(ζₓ+1)/(2*ζₓm13)
end

KHS(ζₓ) = (1-ζₓ)^4/evalpoly(ζₓ,(1,4,4,-4,1))

const _ϕ = SAFTVRMieconsts.ϕ
const _ϕnum = ntuple(m -> (_ϕ[1][m], _ϕ[2][m], _ϕ[3][m], _ϕ[4][m]), Val(6))
const _ϕden = ntuple(m -> (1 + _ϕ[5][m], _ϕ[6][m], _ϕ[7][m]), Val(6))

@inline f123456(α) = ntuple(m -> evalpoly(α, _ϕnum[m]) / evalpoly(α, _ϕden[m]), Val(6))

function g_HS(x_0ij,ζₓ)
    ζX3 = (1-ζₓ)^3
    #evalpoly(ζₓ,(0,42,-39,9,-2)) = (42ζₓ-39ζₓ^2+9ζₓ^3-2ζₓ^4)
    k_0 = -log1p(-ζₓ)+evalpoly(ζₓ,(0,42,-39,9,-2))/(6*ζX3)
    #evalpoly(ζₓ,(0,-12,6,0,1)) = (ζₓ^4+6*ζₓ^2-12*ζₓ)
    k_1 = evalpoly(ζₓ,(0,-12,6,0,1))/(2*ζX3)
    k_2 = -3*ζₓ^2/(8*(1-ζₓ)^2)
    #(-ζₓ^4+3*ζₓ^2+3*ζₓ) = evalpoly(ζₓ,(0,3,3,0,-1))
    k_3 = evalpoly(ζₓ,(0,3,3,0,-1))/(6*ζX3)
    return exp(evalpoly(x_0ij,(k_0,k_1,k_2,k_3)))
end

# ── ∂ζeff/∂ζ_X ─────────────────────────────────────────────────────────────
@inline function ∂ζeff╱∂ζ_X(λ, ζ_X)
    A    = SAFTVRMieconsts.A
    λ⁻¹  = one(λ)/λ
    # ζeff = A·(1·ζ + λ⁻¹·ζ² + λ⁻²·ζ³ + λ⁻³·ζ⁴)
    return A * evalpoly(ζ_X, (one(λ), 2λ⁻¹, 3λ⁻¹*λ⁻¹, 4λ⁻¹*λ⁻¹*λ⁻¹))
end



end #module