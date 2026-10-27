struct Thol2016LJ{T} <: LJPotentialModel
    model::T
end

Thol2016LJ() = Thol2016LJ(TholLJ_constructor())
Base.show(io::IO,model::Thol2016LJ) = print(io,model_showname(model),"()")
Base.show(io::IO,::MIME"text/plain",model::Thol2016LJ) = print(io,model_showname(model),"()")

model_showname(::Thol2016LJ) = "Thol2016LJ"

a_potential(model::Thol2016LJ,V,T) = a_res(model.model,V,T,SA[1.0])
x0_sat_pure(model::Thol2016LJ,T) = x0_sat_pure(model.model,T)
x0_psat(model::Thol2016LJ,T) = x0_psat(model.model,T)
has_fast_crit_pure(::Thol2016LJ) = true
crit_pure(model::Thol2016LJ) = crit_pure(model.model)
x0_volume_liquid(model::Thol2016LJ,p,T,z) = x0_volume_liquid(model.model,p,T,z)
lb_volume(model::Thol2016LJ,T,z) = lb_volume(model.model,T,z)


"""
    Thol2016LJ <: LJPotentialModel
    Thol2016LJ()

Lennard-Jones Reference equation of state. Valid from 0.5 < T/Tc < 7 and pressures up to p/pc = 500.
```
τᵢ = 1.32/T
δᵢ = 1/0.31V
a⁰ᵢ(δ,τ) = log(δᵢ) + 1.5log(τᵢ) - 1.515151515τᵢ + 6.262265814
a⁰(δ,τ,z) = ∑xᵢ(a⁰ᵢ + log(xᵢ))
τ = 1.32ϵ/T
δ = n(Nₐσ^3)/0.31V
aʳ(δ,τ)  = aʳ₁+ aʳ₂ + aʳ₃ + aʳ₄
aʳ₁(δ,τ)  = ∑nᵢδ^(dᵢ)τ^(tᵢ), i ∈ 1:6
aʳ₂(δ,τ)  = ∑nᵢexp(-δ^cᵢ)δ^(dᵢ)τ^(tᵢ), i ∈ 7:12
aʳ₃(δ,τ)  = ∑nᵢexp(-ηᵢ(δ - εᵢ)^2 - βᵢ(τ - γᵢ)^2)δ^(dᵢ)τ^(tᵢ), i ∈ 13:23
```
## References
1. Thol, M., Rutkai, G., Köster, A., Lustig, R., Span, R., & Vrabec, J. (2016). Equation of state for the Lennard-Jones fluid. Journal of physical and chemical reference data, 45(2), 023101. [doi:10.1063/1.4945000](https://doi.org/10.1063/1.4945000)
"""
Thol2016LJ

export Thol2016LJ
