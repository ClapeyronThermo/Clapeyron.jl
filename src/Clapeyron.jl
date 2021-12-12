module Clapeyron
using StaticArrays
using LinearAlgebra

#for the assoc solver
import PackedVectorsOfVectors
const PackedVofV = PackedVectorsOfVectors.PackedVectorOfVectors

#for non allocating vectors of zeros and ones
using FillArrays: FillArrays

using Roots: Roots
using NLSolvers
import Metaheuristics
using DiffResults, ForwardDiff
using Scratch 
include("solvers/Solvers.jl")
using .Solvers
using .Solvers: log, sqrt
using Unitful
import LogExpFunctions
include("constants.jl")
include("models/basetools.jl") #type hierarchy
include("database/ParamOptions.jl") 
include("utils/vectors.jl")
include("utils/ClapeyronParam.jl")


include("utils/macros.jl")
using CSV, Tables
include("database/database.jl")
include("utils/misc.jl")


include("models/combiningrules.jl")

include("models/eos.jl")
include("utils/visualisation.jl")
include("utils/split_model.jl")
include("database/UserReader.jl")

include("models/eos/ideal/BasicIdeal.jl") #before macros, because its used there
include("models/eos/ideal/MonomerIdeal.jl")
include("models/eos/ideal/ReidIdeal.jl")
include("models/eos/ideal/WalkerIdeal.jl")
include("models/eos/ideal/JobackIdeal.jl")

include("models/eos/SAFT/PCSAFT/PCSAFT.jl")
include("models/eos/SAFT/PCSAFT/variants/sPCSAFT.jl")

include("models/eos/SAFT/ogSAFT/ogSAFT.jl")
include("models/eos/SAFT/CPA/CPA.jl")
include("models/eos/SAFT/CPA/variants/sCPA.jl")
include("models/eos/SAFT/SAFTVRSW/SAFTVRSW.jl")
include("models/eos/SAFT/LJSAFT/LJSAFT.jl")
include("models/eos/SAFT/softSAFT/softSAFT.jl")
include("models/eos/SAFT/SAFTVRMie/SAFTVRMie.jl")
include("models/eos/SAFT/SAFTVRMie/variants/SAFTVRQMie.jl")
include("models/eos/SAFT/SAFTgammaMie/SAFTgammaMie.jl")

include("models/eos/SAFT/CKSAFT/CKSAFT.jl")
include("models/eos/SAFT/CKSAFT/variants/sCKSAFT.jl")
include("models/eos/SAFT/BACKSAFT/BACKSAFT.jl")

include("models/eos/SAFT/equations.jl")

include("models/eos/cubic/equations.jl")
include("models/eos/cubic/vdW.jl")
include("models/eos/cubic/RK/RK.jl")
include("models/eos/cubic/PR/PR.jl")

include("models/eos/Activity/Wilson/Wilson.jl")
include("models/eos/Activity/NRTL/NRTL.jl")
include("models/eos/Activity/UNIQUAC/UNIQUAC.jl")
include("models/eos/Activity/UNIFAC/UNIFAC.jl")
include("models/eos/Activity/UNIFAC/variants/ogUNIFAC.jl")
include("models/eos/Activity/UNIFAC/variants/PSRK.jl")
include("models/eos/Activity/UNIFAC/variants/VTPR.jl")

include("models/eos/Activity/COSMOSAC/utils.jl")
include("models/eos/Activity/COSMOSAC/COSMOSAC02.jl")
include("models/eos/Activity/COSMOSAC/COSMOSAC10.jl")
include("models/eos/Activity/COSMOSAC/COSMOSACdsp.jl")

include("models/eos/cubic/alphas/alphas.jl")
include("models/eos/cubic/mixing/mixing.jl")
include("models/eos/cubic/translation/translation.jl")

include("models/eos/cubic/RK/variants/SRK.jl")
include("models/eos/cubic/RK/variants/PSRK.jl")
include("models/eos/cubic/PR/variants/PR78.jl")
include("models/eos/cubic/PR/variants/VTPR.jl")
include("models/eos/cubic/PR/variants/UMRPR.jl")

include("models/eos/Activity/equations.jl")

include("models/eos/EmpiricHelmholtz/IAPWS95.jl")
include("models/eos/EmpiricHelmholtz/PropaneRef.jl")
include("models/eos/EmpiricHelmholtz/GERG2008/GERG2008.jl")

include("models/eos/LatticeFluid/SanchezLacombe/SanchezLacombe.jl")

include("models/eos/SPUNG/SPUNG.jl")

include("models/eos/cached/CachedEoS.jl")

include("methods/methods.jl")

end # module
