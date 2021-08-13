module Clapeyron
using StaticArrays
using LinearAlgebra
using NLopt, NLSolvers,Roots
using  DiffResults, ForwardDiff
using NaNMath: log, pow
include("solvers/Solvers.jl")
using .Solvers
using Unitful
import LogExpFunctions
include("constants.jl")
include("models/basetools.jl")
include("utils/ParamOptions.jl")
include("utils/ClapeyronParam.jl")
include("utils/macros.jl")
using CSV, Tables
include("utils/database.jl")
include("utils/misc.jl")

include("models/combiningrules.jl")

include("models/eos.jl")
include("utils/visualisation.jl")


include("models/eos/ideal/BasicIdeal.jl")
include("models/eos/ideal/MonomerIdeal.jl")
include("models/eos/ideal/ReidIdeal.jl")
include("models/eos/ideal/WalkerIdeal.jl")
include("models/eos/ideal/JobackIdeal.jl")

include("models/eos/SAFT/PCSAFT/PCSAFT.jl")
include("models/eos/SAFT/PCSAFT/variants/sPCSAFT.jl")

include("models/eos/SAFT/ogSAFT/ogSAFT.jl")
include("models/eos/SAFT/CPA/CPA.jl")
include("models/eos/SAFT/SAFTVRSW/SAFTVRSW.jl")
include("models/eos/SAFT/LJSAFT/LJSAFT.jl")
include("models/eos/SAFT/softSAFT/softSAFT.jl")
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
include("models/eos/cubic/equations.jl")

include("models/eos/EmpiricHelmholtz/IAPWS95.jl")
include("models/eos/EmpiricHelmholtz/PropaneRef.jl")
include("models/eos/EmpiricHelmholtz/GERG2008/GERG2008.jl")

include("models/eos/SPUNG/SPUNG.jl")

include("models/eos/cached/CachedEoS.jl")

include("methods/methods.jl")

end # module
