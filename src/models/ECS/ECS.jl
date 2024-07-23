struct ExtendedCorrespondingStates{S<:EoSModel,E<:EoSModel} <: EoSModel
    shape_model::S
    shape_ref::S
    model_ref::E
end

"""
    ExtendedCorrespondingStates <: EoSModel

    function ECS(components,
        refmodel=PropaneRef(),
        shapemodel=SRK(components),
        shaperef = SRK(refmodel.components))
    
## Input Models

- `shape_model`: shape model
- `shape_ref`:  shape reference. is the same type of EoS that `shape_model`
- `model_ref`: Reference model 

## Description

A Extended Corresponding states method.

The idea is to use a "shape model" that provides a corresponding states parameters
and a "reference model" that implements a helmholtz energy function, so that:

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

.1 Mollerup, J. (1998). Unification of the two-parameter equation of state and the principle of corresponding states. Fluid Phase Equilibria, 148(1–2), 1–19. [doi:10.1016/s0378-3812(98)00230-1](https://doi.org/10.1016/s0378-3812(98)00230-1)
"""
const ECS = ExtendedCorrespondingStates

Base.length(model::ECS) = length(model.shape_model)

function eos_impl(model::ECS,V,T,z)
    f,h = shape_factors(model,V,T,z)
    n = sum(z)
    T0 = T/f
    V0 = V/h/n
    #eos(V,T)/RT = eos0(V0,T0)/RT0
    #eos(V,T) = eos0(V0,T0)*T/T0
    return n*eos(model.model_ref,V0,T0)*f
end

function eos_res(model::ECS,V,T,z=SA[1.0])
    f,h = shape_factors(model,V,T,z)
    n = sum(z)
    T0 = T/f
    V0 = V/h/n
    return n*eos_res(model.model_ref,V0,T0)*f
end

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
    b ≈ lb_volume(model,z)
    a ≈ RT*(b - B)
    B = second_virial_coefficient(model,T)
    ```
    This is not tested extensively and it is considered an Experimental feature, subject to future changes.  

## References

1. Mollerup, J. (1998). Unification of the two-parameter equation of state and the principle of corresponding states. Fluid Phase Equilibria, 148(1–2), 1–19. [doi:10.1016/s0378-3812(98)00230-1](https://doi.org/10.1016/s0378-3812(98)00230-1)
    
"""
function shape_factors end
shape_factors(model::ECS,V,T,z=SA[1.0]) = shape_factors(model,model.shape_ref,V,T,z)

function shape_factors(model::ECS,shape_ref::ABCubicModel,V,T,z=SA[1.0])
    a,b = cubic_ab(model.shape_model,V,T,z)
    n = sum(z)
    v = V/n
    # initial point
    amix = dot(z,model.shape_model.params.a.values,z)/(n*n)
    a00 = shape_ref.params.a.values[1,1]
    b00 = shape_ref.params.b.values[1,1]
    
    fT0 = one(T)*b00*amix/a00/b
    function f_0(f)
        a0f,b0f = cubic_ab(shape_ref,v,T/f)
        return f*b/b0f - a/a0f
    end

    prob = Roots.ZeroProblem(f_0,fT0)
    f = Roots.solve(prob)
    _,b0 = cubic_ab(shape_ref,v,T/f)
    h = b/b0
    return f,h
end

mw(model::ECS) = mw(model.shape_model)
molecular_weight(model::ECS,z=SA[1.0]) = molecular_weight(model.shape_model,z)

function Base.show(io::IO,mime::MIME"text/plain",model::ECS)
    println(io,"Extended Corresponding States model")
    println(io," reference model: ",model.model_ref)
    print(io," shape model: ",model.shape_model)
end

function Base.show(io::IO,model::ECS)
    print(io,string(typeof(model)),model.shape_model.components)
end

function lb_volume(model::ECS,z)
    lb_v0 = lb_volume(model.model_ref,z)
    T0 = T_scale(model.model_ref)
    f,h = shape_factors(model,lb_v0,T0,z) #h normaly should be independent of temperature
    return lb_v0*h
end

function x0_volume(model::ECS,T,z)
    lb_v0 = lb_volume(model.model_ref,T,z)
    f,h = shape_factors(model,lb_v0,T,z)
    return lb_v0*h
end

function x0_volume_liquid(model::ECS,p,T,z)
    lb_v0 = lb_volume(model.model_ref,T,z)
    f,h = shape_factors(model,lb_v0,T,z)
    T0 = T/f
    v0l = x0_volume_liquid(model.model_ref,p,T0,z)
    return v0l*h
end

function T_scale(model::ECS,z)
    lb_v0 = lb_volume(model.model_ref)
    T0 = T_scale(model.model_ref)
    f,h = shape_factors(model,lb_v0,T0,z) #h normaly should be independent of temperature
    return T0*f
end

function p_scale(model::ECS,z)
     lb_v0 = lb_volume(model.model_ref)
     T0 = T_scale(model.model_ref)
     p0 = p_scale(model.model_ref)
     f,h = shape_factors(model,lb_v0,T0,z) #h normaly should be independent of temperature
     ps = p0*f/h
     return ps
end

function x0_sat_pure(model::ECS,T)
    f,h = shape_factors(model,zero(T),T) 
    T0 = T/f
    v0l,v0v = x0_sat_pure(model.model_ref,T0)
    return (v0l*h,v0v*h) 
end

function split_model(model::ECS,splitter)
    shape_model_vec = split_model(model.shape_model,splitter)
    shape_ref,model_ref = model.shape_ref, model.model_ref
    return [ECS(shape_modeli,shape_ref,model_ref) for shape_modeli in shape_model_vec]
end

#==
Experimental way of trying to make general shape factors
WIP

It will try to make extended corresponding states by the same criteria as the cubic:
f*h = a(T)/a0(T0)
where h = lb_volume(model)/lb_volume(model0)
a(T) = RT*(b-B(T))
==#
function shape_factors(model::ECS,shape_ref::EoSModel,V,T,z=SA[1.0])
    n = sum(z)
    shape_ref = model.shape_ref
    RT = R̄*T
    b = lb_volume(model.shape_model,T,z)
    b0 = lb_volume(shape_ref,T,SA[1.0])
    n = sum(z)
    v = V/n
    B = second_virial_coefficient(model.shape_model,T,z)
    B00 = second_virial_coefficient(shape_ref,T)
    #B = b-a/RT
    #a/RT = b-B
    a = R̄*T*(b-B)
    a00 = R̄*T*(b0-B00) #f = 1
    h = b/b0
    fT0 = one(T)*b0*a/a00/b
    function f_0(f)
        T0 = T/f
        B0f = second_virial_coefficient(shape_ref,T0)
        a0f = R̄*T0*(b0-B0f)
        return f*h - a/a0f
    end
    prob = Roots.ZeroProblem(f_0,fT0)
    f = Roots.solve(prob)
    return f,h
end

export ExtendedCorrespondingStates
export ECS, shape_factors
