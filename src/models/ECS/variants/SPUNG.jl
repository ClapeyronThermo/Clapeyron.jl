"""

    function function SPUNG(components,
        refmodel=PropaneRef(),
        shapemodel=SRK(components),
        shaperef = SRK(refmodel.components))
    
## Description

SPUNG: _State Research Program for Utilization of Natural Gas_

[`ECS`](@ref) method. It uses [`SRK`](@ref) as the shape model and [`PropaneRef`](@ref) as the reference model.

## References

1. Wilhelmsen, Ø., Skaugen, G., Jørstad, O., & Li, H. (2012). Evaluation of SPUNG* and other equations of state for use in carbon capture and storage modelling. Energy Procedia, 23, 236–245. [doi:10.1016/j.egypro.2012.06.024](https://doi.org/10.1016/j.egypro.2012.06.024)
"""
function SPUNG(components,
    refmodel=PropaneRef(),
    shapemodel = SRK(components),
    shaperef = SRK(component_list(refmodel));idealmodel = BasicIdeal)

    scaling = CubicECSScaling(shapemodel,shaperef)

    model = ECS(scaling,refmodel,init_idealmodel)
    return model
end


struct CubicECSScaling{C <: DeltaCubicModel} <: EoSModel
    shape_model::C
    shape_ref::C
end

component_list(model::CubicECSScaling) = component_list(model.shape_model)
Base.length(model::CubicECSScaling) = length(model.shape_model)

function shape_factors(scaling::CubicECSScaling,reference::EoSModel,V,T,z)
    shape_model = scaling.shape_model
    shape_ref = scaling.shape_ref
    a,b = cubic_ab(shape_model,V,T,z)
    r_shape = _rgas_ratio(model)
    n = sum(z)
    v = V/n
    # initial point
    amix = dot(z,shape_model.params.a.values,z)/(n*n)
    a00 = shape_ref.params.a.values[1,1]
    b00 = shape_ref.params.b.values[1,1]

    fT0 = one(T)*b00*amix/a00/b

    prob = Roots.ZeroProblem(__cubic_shape_factors_prob,fT0)
    prob_data = (shape_ref,v,T,a,b,r_shape)
    f = roots_solve_ad(prob,Roots.Order0(),prob_data)
    _,b0 = cubic_ab(shape_ref,v,T/f)
    h = b/b0
    return f,h
end

function __cubic_shape_factors_prob(f,data)
    shape_ref,v,T,a,b,r_shape = data
    a0f,b0f = cubic_ab(shape_ref,v,T/f)
    return f*b/b0f - a*r_shape/a0f
end

export SPUNG
