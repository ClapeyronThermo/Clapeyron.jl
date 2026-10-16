"""
    donnan_psi_bulk(model::ElectrolyteModel, T, ρbulk_dense::AbstractVector, ρbulk_dilute::AbstractVector)

Read off the Donnan potential `Ψ`:

```julia
`μᵢ(dense) - μᵢ(dilute) = Zᵢ·Rgas·T·Ψ`
```


(dense-minus-dilute convention, matching `donnan_bubble_residual!`'s own `μᵢ(dense) - μᵢ(dilute) = Zᵢ·Rgas·T·Ψ`)
between two ALREADY-KNOWN coexisting bulk states, given directly asper-component densities rather than `(V,x)` pairs -- since chemical
potential only depends on density (`VT_chemical_potential(model,1,T,ρ)` and
`VT_chemical_potential(model,V,T,ρ*V)` agree for any `V`, confirmed
directly), no volume/mole-fraction bookkeeping is needed from the caller.

No iterative solve: unlike `bubble_pressure(...,DonnanBubblePressure())`,
which solves the FULL coupled equilibrium (`Ψ` along with both phases'
volumes/composition), this is a direct readout for a caller that already
has a converged coexisting pair of bulk states in hand (e.g. a DFT solver's
own `structure.ρbulk`/`structure.topology.ρbulk2`, built earlier from a
`bubble_pressure` result) and only needs `Ψ` itself, not to re-solve the
equilibrium that produced it. Averaged over every charged component
(`Zᵢ≠0`) for numerical robustness -- at a genuine coexistence point these
should all agree; returns exactly `0.0` if `model` has no charged
components (nothing to average).
"""
function donnan_psi_bulk(model::ElectrolyteModel, T, ρbulk_dense::AbstractVector, ρbulk_dilute::AbstractVector)
    Z = component_charges(model)
    Ψ = zero(Base.promote_eltype(model,T,ρbulk_dense,ρbulk_dilute))
    μd = VT_chemical_potential(model, 1.0, T, ρbulk_dense)
    μl = VT_chemical_potential(model, 1.0, T, ρbulk_dilute)
    RT = Rgas(model) * T
    ncharged = 0
    for i in 1:length(model)
        Zi = Z[i]
        if Zi != 0
            ncharged += 1
            Ψ += (μd[i] - μl[i]) / (RT * Zi)
        end
    end
    iszero(ncharged) && return Ψ
    return Ψ/ncharged
end

export donnan_psi_bulk
