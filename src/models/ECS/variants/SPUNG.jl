"""
    SPUNGScaling(components, reference_component, model, refmodel;
    shape_userlocations = nothing,
    shaperef_userlocations = nothing,
    verbose = false)

    SPUNGScaling(shape_model, shaperef_model)

## Input models

- `model`: cubic equation of state that describes the fluid (or mixture) of interest, the *shape equation* (for example `SRK` or `PR`). A cubic built beforehand can also be passed; it is then used as it is.
- `refmodel`: cubic equation of state that represents the reference fluid. It must be the **same equation** as `model`: same functional form, same `α(T)` function and same mixing and translation options. Only the substance, and the parameters describing it, may differ (see the notes below).

## Input parameters

- `components`: substances of the fluid or mixture of interest.
- `reference_component`: the reference substance, as named in the cubic-equation database.
- `shape_userlocations`: alternative parameter source (files or values) for the cubic of the fluid of interest. Replaces database critical constants, acentric factors or binary interaction parameters, and therefore changes the scale factors.
- `shaperef_userlocations`: the same, for the cubic of the reference fluid.
- `verbose`: prints information about the parameters used while the model is built. It does not change any result.

## Description

Scaling model of the [`SPUNG`](@ref) equation of state, to be used with [`ECS`](@ref). 
Its role is to say *which state of the reference fluid corresponds to a given state of the fluid of interest*: a state `(T, v, x)` is mapped to the reference state

```
T₀ = T/f
v₀ = v/h
```

where `f` and `h` are the temperature and volume scale factors. They are obtained from the cubic equation of state, written as

```
P = RT/(v-b) - a(T)/((v+δ₁b)(v+δ₂b))
```

Two states are considered *corresponding* when the two dimensionless groups of the cubic coincide,

```
Γ = a(T)/(b·R·T)       β = b/v
```

For a cubic these two groups fully determine the reduced residual Helmholtz energy, so equality of `Γ` and `β` means equal reduced residual Helmholtz energy. This gives

```
h = v/v₀ = b/b₀
f = T/T₀ = [a(T)/b] / [a₀(T₀)/b₀]
```

where `a, b` belong to the fluid of interest and `a₀, b₀` to the reference fluid, both from the same cubic equation.
`f` appears on both sides (`T₀ = T/f`), so it is the solution of an implicit equation. 
For pure fluids described by Soave-type `α(T)` functions (SRK, PR) the solution is available in closed form; in all other cases it is found iteratively. 
Both routes give the same value.

For mixtures, `a` and `b` are the mixture parameters given by the mixing rule of the cubic (van der Waals one-fluid rules by default), so the composition dependence of the scale factors, and the effect of the binary interaction parameters `kᵢⱼ` and `lᵢⱼ`, come entirely from the cubic. 
The scaling model supplies only the mapping; the accuracy of the final Helmholtz energy comes from the reference equation.

## Notes

- If the fluid of interest is the reference fluid, and both cubics use the same parameters, `f = h = 1` and the reference equation is recovered unchanged.
- The two cubics must be the same equation because equality of `Γ` and `β` implies equal Helmholtz energy only if both fluids obey the same functional form (same `δ₁`, `δ₂`, `α(T)` and mixing rule). 
Using different cubic equations for the fluid of interest and the reference fluid is not currently supported.

## References

1. Mollerup, J. (1998). Unification of the two-parameter equation of state and the principle of corresponding states. Fluid Phase Equilibria, 148(1–2), 1–19. [doi:10.1016/s0378-3812(98)00230-1](https://doi.org/10.1016/s0378-3812(98)00230-1)
2. Michelsen, M. L., & Mollerup, J. M. (2007). Thermodynamic Models: Fundamentals and Computational Aspects, 2nd ed. Tie-Line Publications.
"""
struct SPUNGScaling{C <: DeltaCubicModel} <: EoSModel
    shape_model::C
    shape_ref::C
end

Base.show(io::IO,mime::MIME"text/plain",model::SPUNGScaling) = eosshow(io,mime,model)
Base.show(io::IO,model::SPUNGScaling) = eosshow(io,model)
component_list(model::SPUNGScaling) = component_list(model.shape_model)
Base.length(model::SPUNGScaling) = length(model.shape_model)
mw(model::SPUNGScaling) = mw(model.shape_model)
molecular_weight(model::SPUNGScaling,z) = molecular_weight(model.shape_model,z)
function each_split_model(scaling::SPUNGScaling, I)
    # split the fluid(s), keep the single-component reference cubic whole
    return SPUNGScaling(each_split_model(scaling.shape_model, I), scaling.shape_ref)
end

function split_pure_model(scaling::SPUNGScaling, I)
    pure_scaling = each_split_model(scaling.shape_model, I)
    return map(x -> SPUNGScaling(x,scaling.shape_ref),pure_scaling) # reference is pure, never split
end

function shape_factors(scaling::SPUNGScaling,reference::EoSModel,T,z)
    shape_model = scaling.shape_model

    if has_fast_crit_pure(shape_model) && shape_model.alpha isa GeneralizedSuaveAlphaModel && length(z) == 1
        _1 = one(eltype(z))
        f0,h0 =  shape_factors_cubic1_soave(scaling, reference, T)
        return f0*_1,h0*_1
    end
    shape_ref = scaling.shape_ref

    a,b = cubic_ab(shape_model,1.0,T,z)
    n = sum(z)
    #v = V/n
    # initial point
    amix = dot(z,shape_model.params.a.values,z)/(n*n)
    a00 = shape_ref.params.a.values[1,1]
    b00 = shape_ref.params.b.values[1,1]

    fT0 = one(T)*b00*amix/a00/b

    prob = Roots.ZeroProblem(__cubic_shape_factors_prob,fT0)
    prob_data = (shape_ref,T,a,b)
    f = roots_solve_ad(prob,Roots.Order0(),prob_data)
    _,b0 = cubic_ab(shape_ref,1.0,T/f)
    h = b/b0
    return f,h
end

function __cubic_shape_factors_prob(f,data)
    shape_ref,T,a,b = data
    a0f,b0f = cubic_ab(shape_ref,1.0,T/f)
    return f - (a/a0f)*(b0f/b)
end

function shape_factors_cubic1_soave(scaling::SPUNGScaling, reference::EoSModel, T)
    shape_model, shape_ref = scaling.shape_model, scaling.shape_ref
    ac,  bc  = shape_model.params.a.values[1,1], shape_model.params.b.values[1,1]
    a0c, b0c = shape_ref.params.a.values[1,1],   shape_ref.params.b.values[1,1]
    m,  m0   = α_m(shape_model,shape_model.alpha,1), α_m(shape_ref,shape_ref.alpha,1)
    Tc, Tc0  = shape_model.params.Tc.values[1], shape_ref.params.Tc.values[1]
    Tr, Tr0 = T/Tc, T/Tc0
    α = (1 + m*(1 - sqrt(Tr)))^2                   # Tr = T/Tc
    C = ac*b0c*α/(bc*a0c)
    u = (m0*sqrt(Tr0) + sqrt(C))/(1 + m0)  
    return u*u, bc/b0c
end

function inverse_shape_factors(scaling::SPUNGScaling, reference, T0, z)
    shape_model = scaling.shape_model

    if has_fast_crit_pure(shape_model) && shape_model.alpha isa GeneralizedSuaveAlphaModel && length(z) == 1
        _1 = one(eltype(z))
        f0, h0 = inverse_shape_factors_cubic1_soave(scaling, reference, T0)
        return f0*_1, h0*_1
    end

    shape_ref = scaling.shape_ref

    # T₀ is fixed and the reference is pure ⇒ a₀, b₀ are constants of the problem
    a0, b0 = cubic_ab(shape_ref, 1.0, T0)

    # initial guess (f ≈ 1 ⇒ T ≈ T₀)
    amix, b = cubic_ab(shape_model, 1.0, T0, z)
    f_init = one(T0) * b0 * amix / (a0 * b)

    prob = Roots.ZeroProblem(__cubic_inverse_shape_factors_prob, f_init)
    prob_data = (shape_model, T0, z, a0, b0)
    f = roots_solve_ad(prob, Roots.Order0(), prob_data)

    T = T0 * f
    _, bT = cubic_ab(shape_model, 1.0, T, z)
    h = bT / b0
    return f, h
end

function __cubic_inverse_shape_factors_prob(f, data)
    shape_model, T0, z, a0, b0 = data
    T = T0 * f
    a, b = cubic_ab(shape_model, 1.0, T, z)
    return f - (a/a0)*(b0/b)
end

function inverse_shape_factors_cubic1_soave(scaling::SPUNGScaling, reference::EoSModel, T0)
    shape_model, shape_ref = scaling.shape_model, scaling.shape_ref
    ac,  bc  = shape_model.params.a.values[1,1], shape_model.params.b.values[1,1]
    a0c, b0c = shape_ref.params.a.values[1,1],   shape_ref.params.b.values[1,1]
    m,  m0   = α_m(shape_model, shape_model.alpha, 1), α_m(shape_ref, shape_ref.alpha, 1)
    Tc, Tc0  = shape_model.params.Tc.values[1],       shape_ref.params.Tc.values[1]
    Tr0 = T0/Tc0
    α0  = (1 + m0*(1 - sqrt(Tr0)))^2
    K   = ac*b0c/(bc*a0c)
    L   = α0/K
    u   = (1 + m) / (sqrt(L) + m*sqrt(T0/Tc))
    return u*u, bc/b0c
end

function SPUNGScaling(components,reference_component,model,refmodel;shape_userlocations = nothing,shaperef_userlocations = nothing,verbose = false)
    shapemodel = init_model(model,components,shape_userlocations,verbose)
    shaperef   = init_model(refmodel,reference_component,shaperef_userlocations,verbose)
    return SPUNGScaling(shapemodel,shaperef)
end

"""
    SPUNG(components;
    reference_model = PropaneRef(),
    shape_model = SRK,
    shaperef_model = shape_model,
    idealmodel = BasicIdeal,
    reference_component = only(Clapeyron.component_list(reference_model)),
    ideal_userlocations = String[],
    shape_userlocations = String[],
    shaperef_userlocations = String[],
    verbose = false,
    reference_state = nothing)

## Input models
- `reference_model`: accurate pure-fluid equation of state (normally a multiparameter Helmholtz equation) used as athe reference model. It must be already created.
- `shape_model`: cubic equation of state that defines the mapping between the fluid of interest and the reference fluid (for example `SRK` or `PR`).
- `shaperef_model`: cubic equation that represents the reference fluid in the mapping. Default: the same as `shape_model`. It must be the **same model type as `shape_model`**. It is created with components equal to `reference_component`.

## Description

SPUNG: _State Research Program for Utilization of Natural Gas_

[`ECS`](@ref) method. The residual Helmholtz energy of a pure fluid or mixture is taken from a very accurate equation of state for a pure reference fluid, evaluated at a corresponding state:

```
Aʳᵉˢ(T,v,x)/RT = Aʳᵉˢ₀(T₀,v₀)/R₀T₀        T₀ = T/f ,  v₀ = v/h
```

The scale factors `f` and `h` depend on temperature and composition, and are calculated from a cubic equation of state (by default [`SRK`](@ref)) with the [`SPUNGScaling`](@ref) model: `h` is the ratio of the covolumes `b/b₀`, and `f` follows from equating `a(T)/(bRT)` for the fluid of interest and for the reference fluid. 
In mixtures the parameters `a` and `b` come from the cubic's mixing rules. 
If the reference equation is itself the same cubic, SPUNG reduces to that cubic.

The ideal-gas part is independent of the scaling and is given by `idealmodel`.

## Model Construction Examples

```julia
# Default: SRK shape equation, propane reference fluid
model = SPUNG(["methane","carbon dioxide"])

# Another reference fluid (any pure-fluid Helmholtz equation)
model = SPUNG(["carbon dioxide","nitrogen"], reference_model = SingleFluid("methane"))

# Another shape equation (used for the mixture and for the reference fluid)
model = SPUNG(["methane","ethane"], shape_model = PR)

# Realistic ideal-gas heat capacities, for caloric properties
model = SPUNG(["methane","ethane"], idealmodel = ReidIdeal)

# Custom parameters (e.g. SPUNG-specific interaction parameters) for the mixture's cubic
model = SPUNG(["methane","ethane"], shape_userlocations = ["path/to/my/db"])
```

## References

1. Mollerup, J. (1998). Unification of the two-parameter equation of state and the principle of corresponding states. Fluid Phase Equilibria, 148(1–2), 1–19. [doi:10.1016/s0378-3812(98)00230-1](https://doi.org/10.1016/s0378-3812(98)00230-1)
1. Michelsen, M. L., & Mollerup, J. M. (2007). Thermodynamic Models: Fundamentals and Computational Aspects, 2nd ed. Tie-Line Publications.
"""
function SPUNG(components;
        reference_model::EoSModel = PropaneRef(),
        shape_model               = SRK,
        shaperef_model            = shape_model,
        idealmodel                = BasicIdeal,
        reference_component       = only(Clapeyron.component_list(reference_model)),
        ideal_userlocations       = String[],
        shape_userlocations       = String[],
        shaperef_userlocations    = String[],
        verbose                   = false,
        reference_state           = nothing)

    _components = format_components(components)
    _reference_component = format_components(reference_component)
    scaling = SPUNGScaling(_components,_reference_component,shape_model,shaperef_model;shape_userlocations,shaperef_userlocations,verbose)
    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    model = ECS(scaling,reference_model,init_idealmodel)
    set_reference_state!(model,reference_state;verbose)
    return model
end

export SPUNG, SPUNGScaling