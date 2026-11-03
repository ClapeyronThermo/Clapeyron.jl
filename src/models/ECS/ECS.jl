"""
abstract type CorrespondingStatesModel <: EoSModel end

Abstract type containing corresponding states (CS) models. 
corresponding states models are defined as:

```
Aʳᵉˢ(T,v,x)/RT = Aʳᵉˢ₀(T₀,v₀)/R₀T₀        T₀ = T/f ,  v₀ = v/h

```

In particular, the code implementation is just two methods.
```
function a_res(model::CorrespondingStatesModel,V,T,z)
    reference = reference_model(model)
    scaling = scaling_model(model)
    return a_res_ecs(reference,scaling,V,T,z)
end

function a_res_ecs(refmodel,scalingmodel,V,T,z)
    f,h = shape_factors(scalingmodel,refmodel,T,z)
    T0 = T/f
    V0 = V/h
    return a_res(refmodel,V0,T0,z)
end
```

The following functions are required:
- [`reference_model`](@ref) : provides a single-component, helmholtz equation of state, used as reference for the shape their `a_res` V-T surface.
- [`scaling_model`](@ref)   : provides a struct that defines the scaling between the real fluid and the reference fluid.
- [`shape_factors`](@ref)   : the scaling model needs to define `f,h = shape_factors(scaling::MyModel,reference,T,z)`.


Some reference models may carry its own scaling (or no scaling at all). For those cases, you may need to overload `a_res_ecs`
"""
abstract type CorrespondingStatesModel <: EoSModel end

"""
    shape_factors(model::CorrespondingStatesModel,T,z=SA[1.0])
    shape_factors(scaling_model,reference_model,T,z)

returns `f`,`h` scaling factors, used by subtypes of [`CorrespondingStatesModel`](@ref). 
`T` is the scaled temperature. To calculate the shape factors given a reference temperature (`T0 = T/f`), use [`inverse_scaling_factors`](@ref)
"""
function shape_factors end

"""
    inverse_shape_factors(model::CorrespondingStatesModel,T0,z=SA[1.0])
    inverse_shape_factors(scaling_model,reference_model,T0,z)

returns `f`,`h` scaling factors, used by subtypes of [`CorrespondingStatesModel`](@ref). 
`T` is the temperature used by the reference model. To calculate the shape factors at a scaled temperature (`T = T0*f`), use [`scaling_factors`](@ref)
"""
function inverse_shape_factors end

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

shape_factors(model,T) = shape_factors(model,T,SA[1.0])
shape_factors(model,T,z::AbstractVector) = shape_factors(scaling_model(model),reference_model(model),T,z)
shape_factors(model1,model2,T::Number) = shape_factors(model1,model2,T,SA[1.0])

inverse_shape_factors(model,T::Number) = inverse_shape_factors(model,T,SA[1.0])
inverse_shape_factors(model,T,z::AbstractVector) = inverse_shape_factors(scaling_model(model),reference_model(model),T,z)
inverse_shape_factors(model1,model2,T::Number) = inverse_shape_factors(model1,model2,T,SA[1.0])


rgas_ratio(model::CorrespondingStatesModel) = rgas_ratio(scaling_model(model),reference_model(model))
rgas_ratio(model1,model2) = Rgas(model1)/Rgas(model2)

Rgas(model::CorrespondingStatesModel) = Rgas(scaling_model(model))
mw(model::CorrespondingStatesModel) = mw(scaling_model(model))
molecular_weight(model::CorrespondingStatesModel,z) = molecular_weight(scaling_model(model),z)
component_list(model::CorrespondingStatesModel) = component_list(scaling_model(model))

function Base.show(io::IO,mime::MIME"text/plain",model::CorrespondingStatesModel)
    print(io, model_showname(model))
    length(model) == 1 && print(io, " with 1 component:")
    length(model) > 1 && print(io, " with ", length(model), " components:")
    println(io)
    show_pairs(io,component_list(model))
    print(io,'\n'," Reference Model: ", parameterless_type(reference_model(model)))
    #print(io,"(\"",only(component_list(reference_model(model))),"\")")
    print(io,'\n'," Fluid Model: ",typeof(scaling_model(model)))
    show_info(io,model)
    show_params(io,model)
    show_reference_state(io,model)
end

Base.show(io::IO,model::CorrespondingStatesModel) = eosshow(io,model)

function lb_volume(model::CorrespondingStatesModel,T,z)
    lb_v0 = lb_volume(reference_model(model),T,SA[1.0])
    _,h = shape_factors(model,T,z)
    return lb_v0*h*sum(z)
end

function x0_volume_liquid(model::CorrespondingStatesModel, p, T, z)
    f, h = shape_factors(model, T, z)
    r = rgas_ratio(model)
    T0 = T/f
    p0 = p*h/(r*f)
    v0l = x0_volume_liquid(reference_model(model),p0,T0,SA[1.0])
    return h*v0l*sum(z)
end

function x0_volume_gas(model::CorrespondingStatesModel,p,T,z)
    f,h = shape_factors(model,T,z)
    r = rgas_ratio(model)
    T0 = T/f
    p0 = p*h/(f*r)
    v0v = x0_volume_gas(reference_model(model),p0,T0,SA[1.0])
    return v0v*h*sum(z)
end

function volume_impl(model::CorrespondingStatesModel,p,T,z,phase,threaded,_vol0)
    f,h = shape_factors(model,T,z)
    r = rgas_ratio(model)
    T0 = T/f
    p0 = p*h/(f*r)
    n = sum(z)
    if _vol0 !== nothing && !isnan(_vol0)
        vol0 = _vol0/(h*n)
        v0 = volume(reference_model(model),p0,T0,SA[1.0];phase,threaded,vol0)
    else
        v0 = volume(reference_model(model),p0,T0,SA[1.0];phase,threaded)
    end
    return v0*h*n
end

function T_scale(model::CorrespondingStatesModel,z)
    T0 = T_scale(reference_model(model),SA[1.0])
    f, _ = shape_factors(model,T0,z)
    return T0*f
end

function p_scale(model::CorrespondingStatesModel, z)
    T0 = T_scale(reference_model(model), z)
    p0 = p_scale(reference_model(model), z)
    f, h = shape_factors(model, T0, z)
    r = rgas_ratio(model)
    return r*f/h*p0
end

function x0_sat_pure(model::CorrespondingStatesModel, T, crit = nothing)
    f, h = shape_factors(model, T)
    T0 = T/f
    if crit === nothing
        v0l, v0v = x0_sat_pure(reference_model(model), T0)
    else
        Tc, Pc, Vc = crit                    # Clapeyron order is (Tc, Pc, Vc)
        f_c, h_c = shape_factors(model, Tc)  # scale factors at the critical temperature
        r  = rgas_ratio(model)
        Tc0 = Tc/f_c
        Vc0 = Vc/h_c
        Pc0 = Pc*h_c/(r*f_c)                 # same map as p0 = p·h/(r·f) in volume_impl
        v0l, v0v = x0_sat_pure(reference_model(model), T0, (Tc0, Pc0, Vc0))
    end
    return (v0l*h, v0v*h)
end

function x0_psat(model::CorrespondingStatesModel, T)
    f, h = shape_factors(model, T)
    p0 = x0_psat(reference_model(model), T/f)
    r  = rgas_ratio(model)
    return r*f/h*p0
end

function second_virial_coefficient_impl(model::CorrespondingStatesModel,T,z)
    f, h = shape_factors(model, T, z)
    B0 = second_virial_coefficient(reference_model(model),T/f,SA[1.0])
    return B0*h*sum(z)
end

function crit_pure(model::CorrespondingStatesModel)
    Tc0,Pc0,Vc0 = crit_pure(reference_model(model))
    fc,hc = inverse_shape_factors(model,Tc0)
    r = rgas_ratio(model)
    return Tc0*fc, r*fc/hc*Pc0, Vc0*hc
end

struct ECSSaturation{M} <: SaturationMethod
    method::M
end

function init_preferred_method(method::typeof(saturation_pressure),model::CorrespondingStatesModel,kwargs)
    inner = init_preferred_method(method,reference_model(model),kwargs)
    return ECSSaturation(inner)
end

#=
function init_preferred_method(method::typeof(saturation_temperature),model::CorrespondingStatesModel,kwargs)
    inner = init_preferred_method(method,reference_model(model),kwargs)
    return ECSSaturation(inner)
end =#

function saturation_pressure_impl(model::CorrespondingStatesModel,T,methodwrapper::ECSSaturation{M}) where M
    method = methodwrapper.method
    f,h = shape_factors(model,T)
    T0 = T/f
    psat0,vl0,vv0 = saturation_pressure(reference_model(model),T0,method)
    r = rgas_ratio(model)
    return r*f/h*psat0,vl0*h,vv0*h
end

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

Generic instance of a [`CorrespondingStatesModel`](@ref) EoS. It can be used to mix and match scaling and reference models.
"""
ECS

const ECS = ExtendedCorrespondingStates

Base.length(model::ECS) = length(scaling_model(model))
scaling_model(model::ECS) = model.scaling
reference_model(model::ECS) = model.reference

function each_split_model(model::ECS, I)
    scaling = each_split_model(scaling_model(model), I)
    ideal   = each_split_model(model.idealmodel, I)
    return ECS(scaling, reference_model(model), ideal)   # reference is pure, never split
end

function split_pure_model(model::ECS, I)
    scaling = split_pure_model(scaling_model(model), I)
    ideal   = split_pure_model(model.idealmodel, I)
    ref = reference_model(model)
    return map(x -> ECS(x[1],ref,x[2]),zip(scaling,ideal)) # reference is pure, never split
end


#=
Basic Scaling models

=#

@newmodelsingleton NoECSScaling EoSModel

"""
    NoECSScaling <: EoSModel

## Input parameters
None

## Description
Null scaling model, to be used along a reference model. temperature and shape factors are equal to one.
"""
NoECSScaling

shape_factors(scaling::NoECSScaling,reference,T,z) = 1.0,1.0
inverse_shape_factors(scaling::NoECSScaling,reference,T,z) = 1.0,1.0
rgas_ratio(scaling::NoECSScaling,reference) = 1.0

export ExtendedCorrespondingStates
export ECS, shape_factors, inverse_shape_factors, NoECSScaling
