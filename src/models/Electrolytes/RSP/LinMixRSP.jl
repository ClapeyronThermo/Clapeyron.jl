abstract type LinMixRSPModel <: RSPModel end
#
struct LinMixRSPParam <: EoSParam
    dielectric_constant::SingleParam{Float64}
end

struct LinMixRSP <: LinMixRSPModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    params::LinMixRSPParam
    references::Array{String,1}
end

@registermodel LinMixRSP
export LinMixRSP

"""
    LinMixRSP(solvents::Array{String,1},
         ions::Array{String,1};
         userlocations::Vector{String}=[],
         verbose::Bool=false)

## Input parameters
- `dielectric_constant::Float64`: Constant Relative Static Permittivity `[-]`

## Description
This function is used to create a Linear Mixing-Rule Relative Static Permittivity model, for a mixture of solvents, where each solvent has a `dielectric_constant`.
"""
function LinMixRSP(solvents,ions; userlocations::Vector{String}=String[],assoc_userlocations::Vector{String}=String[], verbose::Bool=false)

    components = deepcopy(solvents)
    append!(components,ions)
    icomponents = 1:length(components)

    params = getparams(components, ["Electrolytes/RSP/dielectric.csv","SAFT/PCSAFT/ePCSAFTAdv/dielectric.csv"]; userlocations=userlocations, verbose=verbose)
    e_r = params["dielectric"]
    packagedparams = LinMixRSPParam(e_r)

    references = String[]

    model = LinMixRSP(components, icomponents, packagedparams,references)
    return model
end

function dielectric_constant(model::LinMixRSPModel, V, T, z,_data = nothing)
    ϵᵣ = model.params.dielectric_constant.values

    res = zero(eltype(z))
    ∑z = sum(z)
    for i ∈ @comps
        res += z[i]*ϵᵣ[i]
    end

    return res/∑z
end

is_splittable(::LinMixRSP) = true