
"""
krichevskii_parameter(model::EoSModel, T, crit = nothing)

Calculates the krichevskii parameter,defined as:
```
∂p/∂x₂ |T → Tc₁,V → Vc₁, x₂ → 0
```
where the first component is the solvent and second is the solute.
"""
function krichevskii_parameter(model,crit = nothing)
    binary_component_check(krichevskii_parameter,model)
    solvent,solute = split_model(model)
    if crit === nothing
        crit_solvent = crit_pure(solvent)
    else
        crit_solvent = crit
    end
    Tc,Pc,Vc = crit_solvent
    z = [1 - 1e-30, 1e-30]
    ∂p = VT_molar_gradient(model,Vc,Tc,z,pressure)
    return ∂p[2]
end