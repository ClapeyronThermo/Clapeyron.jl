abstract type LinMixRSPModel <: RSPModel end
#
struct LinMixRSPParam <: EoSParam
    dielectric_constant::SingleParam{Float64}
end

struct LinMixRSP <: LinMixRSPModel
    components::Array{String,1}
    params::LinMixRSPParam
    references::Array{String,1}
end

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
function LinMixRSP(solvents,ions; userlocations=String[], verbose::Bool=false)

    solvents = format_components(solvents)
    ions = format_components(ions)
    components = vcat(solvents, ions)

    userlocations = normalize_userlocations(userlocations)

    params = getparams(components, ["Electrolytes/RSP/dielectric.csv","SAFT/PCSAFT/ePCSAFTAdv/dielectric.csv"]; userlocations=userlocations, verbose=verbose)
    e_r = params["dielectric"]
    packagedparams = LinMixRSPParam(e_r)

    references = String[]

    model = LinMixRSP(components, packagedparams, references)
    return model
end

function dielectric_constant(model::LinMixRSPModel, V, T, z, Z = nothing)
    ϵᵣᵢ = model.params.dielectric_constant.values
    ϵᵣ = dot(ϵᵣᵢ,z)/sum(z)
end

is_splittable(::LinMixRSP) = true