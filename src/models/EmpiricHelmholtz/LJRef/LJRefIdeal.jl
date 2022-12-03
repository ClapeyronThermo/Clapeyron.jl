struct LJRefIdeal <: IdealModel
    components::Vector{String}
    lj::LJRef
end

function LJRefIdeal(components;userlocations=String[], verbose=false)
    lj = LJRef(components;userlocations,verbose)
    return LJRefIdeal(components,lj)
end

function a_ideal(model::LJRefIdeal,V,T,z)
    return a_ideal(model.lj,V,T,z)
end

mw(model::LJRefIdeal) = model.lj.params.Mw.values

@registermodel LJRefIdeal

idealmodel(model::LJRef) = LJRefIdeal(model.components,model)

"""
    LJRefIdeal <: IdealModel
    LJRef(components;
    userlocations=String[],
    verbose=false)

## Input parameters

- `sigma`: Single Parameter (`Float64`) - particle size [Å]
- `epsilon`: Single Parameter (`Float64`) - dispersion energy [`K`]
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`

## Description

Lennard-Jones Reference equation of state. Ideal Part. valid from 0.5 < T/Tc < 7 and pressures up to p/pc = 500.


```
τᵢ = 1.32ϵᵢ/T
δᵢ = n(Nₐσᵢ^3)/0.31V
a⁰ᵢ(δ,τ) = log(δᵢ) + 1.5log(τᵢ) + 1.515151515τᵢ + 6.262265814 
a⁰(δ,τ,z) = ∑xᵢ(a⁰ᵢ + log(xᵢ))

```

`LJRefIdeal` acts as a wrapper of `LJRef` model, you can access it with `LJRef(model::LJRefIdeal)`.

!!! warning "Mutiple component warning"

    The original model was done with only one component in mind. to support multiple components, a VDW 1-fluid mixing rule (shown above) is implemented, but it is not tested.

## References

1. Thol, M., Rutkai, G., Köster, A., Lustig, R., Span, R., & Vrabec, J. (2016). Equation of state for the Lennard-Jones fluid. Journal of physical and chemical reference data, 45(2), 023101. [doi:10.1063/1.4945000](https://doi.org/10.1063/1.4945000)

"""
LJRefIdeal

LJRef(model::LJRefIdeal) = model.lj

export LJRefIdeal

