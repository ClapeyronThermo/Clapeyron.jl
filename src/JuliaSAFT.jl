module JuliaSAFT

using NLopt, NLsolve, DiffResults, ForwardDiff, LinearAlgebra
include("solvers/Solvers.jl")
using .Solvers

using Combinatorics

include("utils/database.jl")
include("methods/constants.jl")
include("methods/param_structs.jl")
include("methods/model_structs.jl")
include("methods/eos/combining_rules.jl")
include("methods/import_params.jl")
include("utils/misc.jl")
include("methods/eos/ideal.jl")

include("methods/eos/SAFT/PCSAFT.jl")
include("methods/eos/SAFT/sPCSAFT.jl")
include("methods/eos/SAFT/SAFTVRMie.jl")
include("methods/eos/SAFT/ogSAFT.jl")

include("methods/eos/eos.jl")

export system

export get_volume, get_Psat, get_Pcrit, get_enthalpy_vap, get_pressure, get_entropy, get_chemical_potential, get_internal_energy, get_enthalpy, get_Gibbs_free_energy, get_Helmholtz_free_energy, get_isochoric_heat_capacity, get_isobaric_heat_capacity, get_thermal_compressibility, get_isentropic_compressibility, get_speed_of_sould, get_isobaric_expansitivity, get_Joule_themson_coefficient


function system(components::Array{String,1}, method::String; kwargs...)
    # possible kwargs... are filepaths for
    # customdatabase_like, customdatabase_unlike, customdatabase_assoc
    raw_params = retrieveparams(components, method)
    set_components = [Set([components[i]]) for i in 1:length(components)]
    if method == "PCSAFT"
        model = PCSAFT(set_components, create_PCSAFTParams(raw_params; kwargs...))
    elseif method == "SAFTVRMie"
        model = SAFTVRMie(set_components, create_SAFTVRMieParams(raw_params; kwargs...))
    elseif method == "ogSAFT"
        model = ogSAFT(set_components, create_ogSAFTParams(raw_params; kwargs...))
    elseif method == "sPCSAFT"
        model = sPCSAFT(set_components, create_sPCSAFTParams(raw_params; kwargs...))
    end
    return model
end

include("methods/getproperties_SAFT.jl")

end # module
