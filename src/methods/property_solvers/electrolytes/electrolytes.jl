include("mT.jl")

"""
    critical_salt_fraction(model::EoSModel, T; s_lo=0.01, s_hi=0.999, max_iters=60)

Critical salt fraction `s_c` for a net-neutral species + symmetric 1:1 salt system: (`x(s) = [1-s, s/2, s/2]`, `s` the total salt mole fraction) at a fixed temperature.

The critical salt fraction is found by doing bisection over [`crit_mix(model,x(s))`](@ref Clapeyron.crit_mix) until `T_c ≈ T`.

## Inputs:
 - `model`, electrolyte model
 - `T`, temperature `[K]`
 - `s_lo`, minimum molar salt fraction
 - `s_hi`, maximum molar salt fraction
 - `max_iters`, maximum iterations

## Outputs:
 - `s_c`, critical salt fraction (molar fraction of the salt)
 - `P_c`, pressure at the critical salt fraction point `[Pa]`
 - `V_c`, volume at the critical salt fraction point `[m³]`
"""
function critical_salt_fraction(model::EoSModel, T::Real; s_lo::Float64=0.01, s_hi::Float64=0.999, max_iters::Int=60)
    _1 = Base.promote_eltype(model,T,s_lo,s_hi)
    
    crit_mix_at(s) = crit_mix(model, [1 - s, s / 2, s / 2])

    Tc_lo, = crit_mix_at(s_lo)
    Tc_hi, = crit_mix_at(s_hi)
    
    _is_positive(Tc_lo) || Tc_lo > T || throw(ArgumentError("no coexistence found even at s_lo=$s_lo (crit_mix Tc=$Tc_lo < T=$T); lower it"))
    _is_positive(Tc_hi) || Tc_hi < T || throw(ArgumentError("still coexisting at s_hi=$s_hi (crit_mix Tc=$Tc_hi > T=$T); raise it"))

    lo, hi = s_lo*_1, s_hi*_1
    sc, Pc, Vc = NaN*_1, NaN*_1, NaN*_1
    for _ in 1:max_iters
        mid = 0.5 * (lo + hi)
        Tc_mid, Pc_mid, Vc_mid = crit_mix_at(mid)
        #=
        `crit_mix` becomes numerically fragile very close to `s→1` (pure salt, a genuinely degenerate limit) 

        - `Tc_mix(s) > T` means still unstable (coexistence exists somewhere reachable from `x(s)` at this `T`);
        - `Tc_mix(s) < T` means stable.
        - a `NaN` return during bisection is treated as "still unstable" (the fragility is on the far side of where `s_c` is expected for any solvent/polymer-containing system).

        =#
        if isnan(Tc_mid) || Tc_mid < T
            hi = mid
        else
            lo = mid
        end
        sc, Pc, Vc = mid, Pc_mid, Vc_mid
    end
    return (sc, Vc, Pc)
end

export critical_salt_fraction