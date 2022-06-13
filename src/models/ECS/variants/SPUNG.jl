

"""

    function function SPUNG(components::Vector{String},
        refmodel=PropaneRef(),
        shapemodel=SRK(components),
        shaperef = SRK(refmodel.components))
    
## Description

SPUNG: State Research Program for Utilization of Natural Gas

[`ECS`](@ref) method. It uses [`SRK`] as the shape model and (`PropaneRef`) as the reference model.

## References

1. Wilhelmsen, Ø., Skaugen, G., Jørstad, O., & Li, H. (2012). Evaluation of SPUNG* and other equations of state for use in carbon capture and storage modelling. Energy Procedia, 23, 236–245. doi:10.1016/j.egypro.2012.06.024
"""
function SPUNG(components::Vector{String},
    refmodel=PropaneRef(),
    shapemodel=SRK(components),
    shaperef = SRK(refmodel.components))
    model = ECS(shapemodel,shaperef,refmodel)
    return model
end

export SPUNG