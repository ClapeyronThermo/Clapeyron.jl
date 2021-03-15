module OpenSAFT
using StaticArrays
using LinearAlgebra
using NLopt, NLsolve, NLSolvers,Roots
using  DiffResults, ForwardDiff
include("solvers/Solvers.jl")
using .Solvers
using Combinatorics
using NamedArrays
using DataStructures: DefaultDict
include("constants.jl")
include("utils/macros.jl")
include("utils/OpenSAFTParam.jl")
using CSV, Tables
include("utils/database.jl")
include("utils/misc.jl")

include("models/basetools.jl")
include("models/combiningrules.jl")

include("models/eos.jl")
include("utils/visualisation.jl")


include("models/eos/ideal/BasicIdeal.jl")
include("models/eos/ideal/MonomerIdeal.jl")
include("models/eos/ideal/ReidIdeal.jl")
include("models/eos/ideal/WalkerIdeal.jl")

include("models/eos/SAFT/PCSAFT/PCSAFT.jl")
include("models/eos/SAFT/PCSAFT/variants/sPCSAFT.jl")

include("models/eos/SAFT/SAFTVRMie/SAFTVRMie.jl")
include("models/eos/SAFT/SAFTVRMie/variants/SAFTVRQMie.jl")
include("models/eos/SAFT/SAFTgammaMie/SAFTgammaMie.jl")

include("models/eos/SAFT/CKSAFT/CKSAFT.jl")
include("models/eos/SAFT/CKSAFT/variants/sCKSAFT.jl")

include("models/eos/SAFT/BACKSAFT/BACKSAFT.jl")


include("models/eos/cubic/vdW.jl")
include("models/eos/cubic/RK.jl")
include("models/eos/cubic/SRK.jl")
include("models/eos/cubic/PR.jl")
include("models/eos/cubic/CPA.jl")
include("models/eos/cubic/equations.jl")

include("models/eos/EmpiricHelmholtz/IAPWS95.jl")

# include("models/param_structs.jl")
# include("models/ideal_param_structs.jl")
# include("models/model_structs.jl")
# include("models/eos/combining_rules.jl")
# include("models/import_params.jl")
# include("models/import_ideal_params.jl")

# include("models/eos/Ideal/Basic.jl")
# include("models/eos/Ideal/Monomer.jl")
# include("models/eos/Ideal/Walker.jl")
# include("models/eos/Ideal/Reid.jl")

# include("models/eos/SAFT/PCSAFT.jl")
# include("models/eos/SAFT/CKSAFT.jl")
# include("models/eos/SAFT/sCKSAFT.jl")
# include("models/eos/SAFT/BACKSAFT.jl")
# include("models/eos/SAFT/sPCSAFT.jl")
# include("models/eos/SAFT/SAFTVRMie.jl")
# include("models/eos/SAFT/SAFTVRQMie.jl")
# include("models/eos/SAFT/SAFTVRSW.jl")
# include("models/eos/SAFT/ogSAFT.jl")
# include("models/eos/SAFT/LJSAFT.jl")
# include("models/eos/SAFT/softSAFT.jl")
# include("models/eos/SAFT/SAFTgammaMie.jl")

# include("models/eos/Cubic/vdW.jl")
# include("models/eos/Cubic/RK.jl")
# include("models/eos/Cubic/SRK.jl")
# include("models/eos/Cubic/PR.jl")
# include("models/eos/Cubic/CPA.jl")

# include("models/eos/eos.jl")

# include("models/system.jl")
# include("models/system2.jl")

# export system,System
# export eos,ideal

using Unitful

include("methods/getproperties.jl")
# include("methods/getproperties_SAFT_Unitful.jl")
include("methods/initial_guess_properties.jl")

# export get_volume, get_sat_pure, get_bubble_pressure, get_crit_pure, get_enthalpy_vap, get_pressure, get_entropy, get_chemical_potential, get_internal_energy, get_enthalpy, get_Gibbs_free_energy, get_Helmholtz_free_energy, get_isochoric_heat_capacity, get_isobaric_heat_capacity, get_isothermal_compressibility, get_isentropic_compressibility, get_speed_of_sound, get_isobaric_expansivity, get_Joule_Thomson_coefficient, get_second_virial_coeff
# export eos
# export create_z
end # module
