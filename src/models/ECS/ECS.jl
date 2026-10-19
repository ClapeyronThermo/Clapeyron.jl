abstract type CorrespondingStatesModel <: EoSModel end

scaling_model(model::CorrespondingStatesModel) = model.scaling
reference_model(model::CorrespondingStatesModel) = model.reference

struct ExtendedCorrespondingStates{S,R,I} <: CorrespondingStatesModel
    scaling::S
    reference::R
    idealmodel::I
end

"""
    ExtendedCorrespondingStates <: CorrespondingStatesModel

    function ECS(scaling,reference,idealmodel)

## Input Models

- `scaling`: determines the scales of temperature and volume for the model.
- `reference`: determines the shape of the helmholtz surface (V,T).
- `idealmodel`: determines the ideal gas model.
## Description

A Extended Corresponding states method.

The idea is to use a "shape model" that provides a corresponding states parameters
and a "reference model" that implements a Helmholtz energy function, so that:

```
eos(shape_model,v,T,x)/RT = eos(model_ref,v₀,T₀)/RT₀
```

where:
```
T₀ = T/f
v₀ = v/h
f,h = shape_factors(model::ECS,shape_ref::EoSModel,V,T,z)
```
[`shape_factors`](@ref) can be used to create custom Extended Corresponding state models.

## References

1. Mollerup, J. (1998). Unification of the two-parameter equation of state and the principle of corresponding states. Fluid Phase Equilibria, 148(1–2), 1–19. [doi:10.1016/s0378-3812(98)00230-1](https://doi.org/10.1016/s0378-3812(98)00230-1)
"""
ECS

const ECS = ExtendedCorrespondingStates

Base.length(model::ECS) = length(scaling_model(model))


#=
function eos_res(model::ECS,V,T,z=SA[1.0])
    f,h = shape_factors(model,T,z)
    n = sum(z)
    T0 = T/f
    V0 = V/h/n
    return n*eos_res(reference_model(model),V0,T0)*f
end =#

function a_res(model::CorrespondingStatesModel,V,T,z)
    reference = reference_model(model)
    scaling = scaling_model(model)
    return a_res_ecs(reference,scaling,V,T,z)
end

#we use this dispatch to allow potential models to define their own a_res dispatch on data.
function a_res_ecs(refmodel,scalingmodel,V,T,z)
    f,h = shape_factors(scalingmodel,refmodel,T,z)
    T0 = T/f
    V0 = V/h
    return a_res(refmodel,V0,T0,z)
end


rgas_ratio(model::CorrespondingStatesModel) = rgas_ratio(scaling_model(model),reference_model(model))
rgas_ratio(model1,model2) = Rgas(model1)/Rgas(model2)

"""
    shape_factors(model::ECS,V,T,z=SA[1.0])
    shape_factors(model::ECS,shape_ref::ABCubicModel,V,T,z=SA[1.0])
    shape_factors(model::ECS,shape_ref::EoSModel,V,T,z=SA[1.0])

Returns `f` and `h` scaling factors, used by the [`ECS`](@ref) Equation of state.
```
eos(shape_model,v,T,x)/RT = eos(model_ref,v₀,T₀)/RT₀
```

where:
```
T₀ = T/f
v₀ = v/h
```

For cubics, a general procedure is defined in [1]:

```
h = b/b₀
fh = a(T)/a₀(T₀)

```

!!! info "General Shape Factors?"
    For general EoS, there is no existent publications on how to obtain shape factors. However, we can "map" any EoS to a cubic with:
    ```
    b ≈ lb_volume(model,T,z)
    a ≈ RT*(b - B)
    B = second_virial_coefficient(model,T)
    ```
    This is not tested extensively and it is considered an Experimental feature, subject to future changes.

## References

1. Mollerup, J. (1998). Unification of the two-parameter equation of state and the principle of corresponding states. Fluid Phase Equilibria, 148(1–2), 1–19. [doi:10.1016/s0378-3812(98)00230-1](https://doi.org/10.1016/s0378-3812(98)00230-1)

"""
shape_factors

shape_factors(model,T) = shape_factors(model,T,SA[1.0])
shape_factors(model::EoSModel,T,z) = shape_factors(scaling_model(model),reference_model(model),T,z)
shape_factors(model1::EoSModel,model2::EoSModel,T) = shape_factors(model1,model2,T,SA[1.0])

Rgas(model::ECS) = Rgas(scaling_model(model))
mw(model::ECS) = mw(scaling_model(model))
molecular_weight(model::ECS,z) = molecular_weight(scaling_model(model),z)
component_list(model::ECS) = component_list(scaling_model(model))

function Base.show(io::IO,mime::MIME"text/plain",model::ECS)
    println(io,"Extended Corresponding States model")
    println(io," reference model: ",reference_model(model))
    print(io," scaling model: ",scaling_model(model))
end

function Base.show(io::IO,model::ECS)
    print(io,string(typeof(model)),scaling_model(model).components)
end

function lb_volume(model::CorrespondingStatesModel,T,z)
    lb_v0 = lb_volume(reference_model(model),T,z)
    _,h = shape_factors(model,T,z)
    return lb_v0*h
end

function x0_volume_liquid(model::CorrespondingStatesModel, p, T, z)
    f, h = shape_factors(model, T, z)
    r = rgas_ratio(model)
    p0 = p*h/(r*f)
    return h*x0_volume_liquid(reference_model(model), p0, T/f, z)
end

function x0_volume_gas(model::CorrespondingStatesModel,p,T,z)
    f,h = shape_factors(model,T,z)
    r = rgas_ratio(model)
    T0 = T/f
    p0 = p*h/(f*r)
    v0v = x0_volume_gas(reference_model(model),p0,T0,z)
    return v0v*h
end

function T_scale(model::CorrespondingStatesModel,z)
    T0 = T_scale(reference_model(model),SA[1.0])
    f,h = shape_factors(model,T0,z)
    return T0*f
end

function p_scale(model::CorrespondingStatesModel, z)
    T0 = T_scale(reference_model(model), z)
    p0 = p_scale(reference_model(model), z)
    f, h = shape_factors(model, T0, z)
    r = rgas_ratio(model)
    return r*f/h*p0
end

function x0_sat_pure(model::CorrespondingStatesModel,T,crit=nothing)
    f,h = shape_factors(model,T)
    T0 = T/f
    if crit === nothing
        v0l,v0v = x0_sat_pure(reference_model(model),T0)
    else
        Tc,Vc,Pc = crit
        # f, h should be evaluated at Tc for the critical-point scaling
        f_c,h_c = shape_factors(model,Tc)
        r = rgas_ratio(model)
        Tc0, Vc0 = Pc*h_c/(r*f_c)
        Pc0 = Pc*h_c/(r*f_c)                        # from p = r·f/h·p0
        crit0 = (Tc0, Vc0, Pc0)
        v0l,v0v = x0_sat_pure(reference_model(model),T0,crit0)
    end
    return (v0l*h,v0v*h)
end

function each_split_model(model::ECS,I)
    shape_model = each_split_model(scaling_model(model),I)
    imodel = each_split_model(model.idealmodel,I)
    shape_ref,model_ref = model.shape_ref, reference_model(model)
    ECS(shape_model,shape_ref,model_ref,imodel)
end

#=
Basic Scaling models

=#

@newmodelsingleton NoECSScaling EoSModel

shape_factors(scaling::NoECSScaling,reference::EoSModel,V,T,z) = 1.0,1.0
rgas_ratio(scaling::NoECSScaling,reference::EoSModel) = 1.0

export ExtendedCorrespondingStates
export ECS, shape_factors, NoECSScaling
