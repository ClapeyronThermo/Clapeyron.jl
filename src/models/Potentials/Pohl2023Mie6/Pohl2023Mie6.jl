
const Pohl2023Mie6consts = (
    c1_pol = SA[-0.0192944, 1.38,    -2.2653, 1.6291, -1.974, 0.40412 ],
    c1_exp = SA[ 0.1845,   -0.3227,   1.1351,  2.232,  -2.344, -0.4238],
    c1_gbs = SA[-4.367,     0.0371,   1.3895,  2.835                  ],
    c2_pol = SA[ 0.26021,  -5.525,    8.329,  -19.492, 25.8,   -3.8133],
    c2_exp = SA[-5.05,      2.7842,  -9.523,  -30.383, 17.902,  2.2264],
    c2_gbs = SA[48.445,    -5.506,  -11.643, -24.36                   ],
    t_pol  = SA[1.0,   0.236, 0.872, 0.313, 0.407, 0.703],
    t_exp  = SA[1.78,  2.99,  2.866, 1.20,  3.06,  1.073],
    t_gbs  = SA[1.50,  1.03,  4.02,  1.57 ],
    d_pol  = SA[4.0,   1.0,   1.0,   2.0,   2.0,   3.0  ],
    d_exp  = SA[1.0,   1.0,   3.0,   2.0,   2.0,   5.0  ],
    d_gbs  = SA[2.0,   3.0,   2.0,   2.0                ],
    l      = SA[1.0,   2.0,   2.0,   1.0,   2.0,   1.0  ],
    η      = SA[0.362, 0.313, 1.17,  0.957              ],
    β      = SA[0.0761, 0.143, 0.63, 1.32               ],
    γ      = SA[1.55,  -0.0826, 1.505, 1.07             ],
    ε      = SA[-1.0,  -1.0,   -0.195, -0.287           ],
)
#=
function PohlMie_constructor(λ)
    (11.0 ≤ lambda_r ≤ 13.0) ||@warn "λr = $lambda_r is outside the recommended range [11,13]"

    c1_pol = Pohl2023Mie6consts.c1_pol
    c1_exp = Pohl2023Mie6consts.c1_exp
    c1_gbs = Pohl2023Mie6consts.c1_gbs
    c2_pol = Pohl2023Mie6consts.c2_pol
    c2_exp = Pohl2023Mie6consts.c2_exp
    c2_gbs = Pohl2023Mie6consts.c2_gbs

    t_pol = [1.0,   0.236, 0.872, 0.313, 0.407, 0.703]
    t_exp = [1.78,  2.99,  2.866, 1.20,  3.06,  1.073]
    t_gbs = [1.50,  1.03,  4.02,  1.57 ]
    d_pol = [4.0,   1.0,   1.0,   2.0,   2.0,   3.0  ]
    d_exp = [1.0,   1.0,   3.0,   2.0,   2.0,   5.0  ]
    d_gbs = [2.0,   3.0,   2.0,   2.0                ]
    p     = [1.0,   2.0,   2.0,   1.0,   2.0,   1.0  ]
    η     = [0.362, 0.313, 1.17,  0.957              ]
    β     = [0.0761, 0.143, 0.63, 1.32               ]
    γ     = [1.55,  -0.0826, 1.505, 1.07             ]
    ε     = [-1.0,  -1.0,   -0.195, -0.287           ]

    cp0_R = 2.5
    cv0_R = cp0_R - 1.0

    # ---------- λ-dependent coefficients ----------
    n_pol = c1_pol .+ c2_pol ./ lambda_r
    n_exp = c1_exp .+ c2_exp ./ lambda_r
    n_gbs = c1_gbs .+ c2_gbs ./ lambda_r

    # ---------- empirical critical point (Pohl) ----------
    Tc   = 0.668 + 6.84/lambda_r + 145/lambda_r^3
    rhoc = 0.2516 + 0.049*log10(lambda_r)

    # ---------- residual term (PolExpGauss) ----------
    n = vcat(n_pol, n_exp, n_gbs)
    t = vcat(t_pol, t_exp, t_gbs)
    d = vcat(d_pol, d_exp, d_gbs)
    l = copy(p)                    # exponent δ^l for the exp terms
    e = ones(length(p))            # coefficients in exp(-e δ^l) = 1
    polexpgauss = PolExpGaussTerm(n, t, d, l, e, η, β, γ, ε)
    residual = SingleFluidResidualParam(; polexpgauss = polexpgauss)

    # ---------- ideal-gas coefficients ----------
    # Reference state used by calc_ref_stat in main.cxx:
    #   t_ref = Tc * (0.8/1.32),  p_ref = rhoc * (0.00125/0.31) * t_ref
    # τ_ref = Tc / t_ref = 1.65,  δ_ref = 0.00125/0.31   (λ-independent)
    t_ref     = Tc * (0.8 / 1.32)
    p_ref     = rhoc * (0.00125 / 0.31) * t_ref
    rho_ref   = p_ref / t_ref
    tauref    = Tc / t_ref
    delta_ref = rho_ref / rhoc

    a2 = -cp0_R / tauref                                    # ≡ c1
    a1 = cv0_R - log(delta_ref) - cv0_R * log(tauref)       # ≡ c2

    ideal = SingleFluidIdealParam(a1, a2, cp0_R, Float64[], Float64[])

    # ---------- critical pressure from the EOS ----------
    #   Pc = ρc·Tc·(1 + A01(τ=1, δ=1))
    A01_crit = sum(n_pol .* d_pol) +
               sum(n_exp .* (d_exp .- p) .* exp(-1.0)) +
               sum((d_gbs .- 2.0 .* η .* (1.0 .- ε)) .* n_gbs .*
                   exp.(-η .* (1.0 .- ε).^2 .- β .* (1.0 .- γ).^2))
    Pc = (1.0 + A01_crit) * Tc * rhoc

    # ---------- fluid properties ----------
    Mw         = 1.0
    lb_volume  = 1 / 0.6
    Ttp        = 0.661  * Tc   / 1.32    # scaled from LJ reference
    ptp        = 0.0012 * Pc   / 0.13006
    rhov_tp    = 0.0017 * rhoc / 0.31
    rhol_tp    = 0.0845 * rhoc / 0.31
    acentric   = NaN
    Rgas       = 1.0

    properties = SingleFluidProperties(Mw, Tc, rhoc, lb_volume,
                                       Tc, Pc, rhoc, Ttp, ptp,
                                       rhov_tp, rhol_tp, acentric, Rgas)

    components = ["Mie($(lambda_r),$(lambda_a)) fluid (Pohl)"]

    # Ancillaries are not provided by the correlation; a dummy
    # CompositeModel is used so that the type is valid.
    if isnothing(ancillaries)
        ancillaries = CompositeModel(components;
                                     gas=nothing, liquid=nothing,
                                     saturation=nothing)
    end

    references = ["Pohl, Fingerhut, Thol, Vrabec, Span, J. Chem. Phys. (2021) " *
                  "- Mie(λr,6) EOS, λr ∈ [11,13]"]

    return SingleFluid(components, properties, ancillaries, ideal, residual, references)
end
=#
struct Pohl2023Mie6{T} <: Mie6PotentialModel
    lambda_r::T
end

function a_potential(model::Pohl2023Mie6,V,T)
    λr = model.lambda_r
    λrinv = 1/λr
    Tr   = 0.668 + λrinv*(6.84 + 145*λrinv*λrinv)
    rhor = 0.2516 + 0.049*log10(λr)
    τ   = Tr / T
    δ   =  1/(V*rhor)
    reduced_a_res_pohl(model::Pohl2023Mie6,δ,τ,λrinv)
end

function reduced_a_res_pohl(model::Pohl2023Mie6,δ,τ,λrinv)
    lnδ = log(δ)
    lnτ = log(τ)
    
    _0 = zero(lnδ+lnτ+λrinv)
    αᵣ = _0
    npol = Pohl2023Mie6consts.c1_pol .+ λrinv .* Pohl2023Mie6consts.c2_pol
    αᵣ += term_ar_pol(δ,τ,lnδ,lnτ,_0,npol,Pohl2023Mie6consts.t_pol,Pohl2023Mie6consts.d_pol)

    nexp = Pohl2023Mie6consts.c1_exp .+ λrinv .* Pohl2023Mie6consts.c2_exp
    αᵣ += term_ar_exp(δ,τ,lnδ,lnτ,_0,nexp,Pohl2023Mie6consts.t_exp,Pohl2023Mie6consts.d_exp,Pohl2023Mie6consts.l)

    ngbs = Pohl2023Mie6consts.c1_gbs .+ λrinv .* Pohl2023Mie6consts.c2_gbs
    αᵣ += term_ar_gauss(δ,τ,lnδ,lnτ,_0,ngbs,Pohl2023Mie6consts.t_gbs,Pohl2023Mie6consts.d_gbs,Pohl2023Mie6consts.η,Pohl2023Mie6consts.β,Pohl2023Mie6consts.γ,Pohl2023Mie6consts.ε)
    return αᵣ
end

lb_volume(model::Pohl2023Mie6,T,z) = sum(z)

function pohl_vt_crit(model::Pohl2023Mie6)
    λr = model.lambda_r
    λrinv = 1/λr
    Tr   = 0.668 + λrinv*(6.84 + 145*λrinv*λrinv)
    rhor = 0.2516 + 0.049*log10(λr)
    return Tr,rhor
end

function x0_crit_pure(model::Pohl2023Mie6,z)
    Tr,rhor = pohl_vt_crit(model)
    return Tr,-log10(rhor)
end

model_showname(::Pohl2023Mie6) = "Pohl2023Mie6"
default_references(::Pohl2023Mie6) = ["10.1080/00268979300100411"]

"""
    Pohl2023Mie6 <: MiePotentialModel
    Pohl2023Mie6(λr)

## Description

Multiparameter Equation of state for the ν-6 Mie fluid, with ν (c) between 11 and 13.


## References
Pohl, S., Fingerhut, R., Thol, M., Vrabec, J., & Span, R. (2023). Equation of state for the Mie (λr,6) fluid with a repulsive exponent from 11 to 13. The Journal of Chemical Physics, 158(8), 084506. [doi:10.1063/5.0133412](https://doi.org/10.1063/5.0133412)
"""
Pohl2023Mie6

export Pohl2023Mie6