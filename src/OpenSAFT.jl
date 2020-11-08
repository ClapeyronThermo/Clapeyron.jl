module OpenSAFT

using NLopt, NLsolve, DiffResults, ForwardDiff, LinearAlgebra
include("solvers/Solvers.jl")
using .Solvers
import Unitful
using Combinatorics

include("constants.jl")
include("utils/extractdatabase.jl")
include("models/param_structs.jl")
include("models/model_structs.jl")
include("models/eos/combining_rules.jl")
include("models/import_params.jl")
include("utils/misc.jl")
include("models/eos/SAFT/ideal.jl")

include("models/eos/SAFT/PCSAFT.jl")
include("models/eos/SAFT/sPCSAFT.jl")
include("models/eos/SAFT/SAFTVRMie.jl")
include("models/eos/SAFT/ogSAFT.jl")

include("models/eos/eos.jl")

include("models/system.jl")

export system

include("methods/getproperties_SAFT.jl")

export get_volume, get_sat_pure, get_crit_pure, get_enthalpy_vap, get_pressure, get_entropy, get_chemical_potential, get_internal_energy, get_enthalpy, get_Gibbs_free_energy, get_Helmholtz_free_energy, get_isochoric_heat_capacity, get_isobaric_heat_capacity, get_thermal_compressibility, get_isentropic_compressibility, get_speed_of_sound, get_isobaric_expansitivity, get_Joule_Thomson_coefficient

end # module
