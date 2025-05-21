abstract type ZuoFurstModel <: RSPModel end

struct ZuoFurst <: ZuoFurstModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    references::Array{String,1}
end

@registermodel ZuoFurst
export ZuoFurst

"""
    ZuoFurst(solvents::Array{String,1},
         ions::Array{String,1};
         userlocations::Vector{String}=[],
         verbose::Bool=false)

## Description
This function is used to create a Zuo-First model. The Zuo-First expression estimates the dielectric constant of water for a range of temperatures.

## References
1. Zuo, Y-X., Fürst, W. (1997). Prediction of vapor pressure for nonaqueous electrolyte solutions using an electrolyte equation of state, Fluid Phase Equilibria, 138(1-2), 87-104.
"""
function ZuoFurst(solvents,ions; userlocations::Vector{String}=String[], verbose::Bool=false)
    components = deepcopy(solvents)
    append!(components,ions)
    icomponents = 1:length(components)

    references = String[]

    model = ZuoFurst(components, icomponents ,references)
    return model
end

function dielectric_constant(model::ZuoFurstModel,V,T,z,_data=nothing)
    return -19.2905+29814.5/T-0.019678*T+1.318e-4*T^2-3.1144e-7*T^3
end

is_splittable(::ZuoFurst) = false