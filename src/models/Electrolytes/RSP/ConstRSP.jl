abstract type ConstRSPModel <: RSPModel end

struct ConstRSP <: ConstRSPModel
    value::Float64
end

export ConstRSP

"""
    ConstRSP(solvents::Array{String,1},
        ions::Array{String,1};
        userlocations::Vector{String}=[],
        value::Float64 = 78.38484961,
        verbose::Bool=false)

    ConstRSP(val::Float64)

## Input parameters
- `value::Float64`: Constant Relative Static Permittivity `[-]`

## Description
This function is used to create a constant Relative Static Permittivity model, given by `value`.
"""
function ConstRSP(solvents,ions; userlocations=nothing, value=78.38484961, verbose::Bool=false)
    return ConstRSP(value)
end

ConstRSP() = ConstRSP(78.38484961)

function dielectric_constant(model::ConstRSPModel,V,T,z,Z = nothing)
    return model.value
end

is_splittable(::ConstRSP) = false