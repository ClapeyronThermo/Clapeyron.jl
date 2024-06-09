abstract type ConstRSPModel <: RSPModel end

struct ConstRSP <: ConstRSPModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    value::Float64
    references::Array{String,1}
end

@registermodel ConstRSP
export ConstRSP

"""
    ConstRSP(solvents::Array{String,1}, 
         ions::Array{String,1}; 
         userlocations::Vector{String}=[], 
         value::Float64 = 78.38484961, 
         verbose::Bool=false)

## Input parameters
- `value::Float64`: Constant Relative Static Permittivity `[-]`

## Description
This function is used to create a constant Relative Static Permittivity model, given by `value`.
"""
function ConstRSP(solvents,ions; userlocations::Vector{String}=String[], value =  78.38484961, verbose::Bool=false)
    components = deepcopy(solvents)
    append!(components,ions)
    icomponents = 1:length(components)

    references = String[]
    
    model = ConstRSP(components, icomponents, value ,references)
    return model
end

function dielectric_constant(model::ConstRSPModel,V,T,z,_data=nothing)
    return model.value
end

is_splittable(::ConstRSP) = false