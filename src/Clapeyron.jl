module Clapeyron
using StaticArrays
using LinearAlgebra
using SparseArrays
#for the assoc solver and the sparse packed VofV
import PackedVectorsOfVectors
const PackedVofV = PackedVectorsOfVectors.PackedVectorOfVectors

#for non allocating vectors of zeros and ones
using FillArrays: FillArrays

using Roots: Roots
using NLSolvers
import BlackBoxOptim
using DiffResults, ForwardDiff
using Scratch 

include("solvers/Solvers.jl")
using .Solvers
using .Solvers: log, sqrt
∂Tag = Solvers.∂Tag
using Unitful
import LogExpFunctions

#The Base of Clapeyron: EoSModel and eos(model,V,T,z)
include("base/constants.jl") 
include("base/EoSModel.jl") 
include("base/eosshow.jl")

#EoSParam, ClapeyronParam
include("params/paramvectors.jl")
include("params/ClapeyronParam.jl")
include("params/combiningrules.jl")

using CSV, Tables
#getparams machinery
include("database/ParamOptions.jl") 
include("database/database.jl")
include("database/UserReader.jl")

#macros, used for defining models
include("utils/macros.jl")

#splitting models, useful for methods.
include("utils/split_model.jl")

#Clapeyron methods (AD, property solvers, etc)
include("methods/methods.jl")

#=
the dependency chain is the following:

base --> params -|-> database ----------------|
                 |-> split_model --> methods -|-> models
                 |-> macros ------------------|

=#

#Clapeyron EoS collection
include("models/types.jl") #type hierarchy
include("models/ideal/ideal.jl")
include("models/ideal/BasicIdeal.jl")
include("models/ideal/MonomerIdeal.jl")
include("models/ideal/ReidIdeal.jl")
include("models/ideal/WalkerIdeal.jl")
include("models/ideal/JobackIdeal.jl")

include("models/SAFT/PCSAFT/PCSAFT.jl")
include("models/SAFT/PCSAFT/variants/sPCSAFT.jl")
include("models/SAFT/ogSAFT/ogSAFT.jl")
include("models/SAFT/CPA/CPA.jl")
include("models/SAFT/CPA/variants/sCPA.jl")
include("models/SAFT/SAFTVRSW/SAFTVRSW.jl")
include("models/SAFT/LJSAFT/LJSAFT.jl")
include("models/SAFT/softSAFT/softSAFT.jl")
include("models/SAFT/SAFTVRMie/SAFTVRMie.jl")
include("models/SAFT/SAFTVRMie/variants/SAFTVRQMie.jl")
include("models/SAFT/SAFTgammaMie/SAFTgammaMie.jl")
include("models/SAFT/CKSAFT/CKSAFT.jl")
include("models/SAFT/CKSAFT/variants/sCKSAFT.jl")
include("models/SAFT/BACKSAFT/BACKSAFT.jl")
include("models/SAFT/equations.jl")

include("models/cubic/equations.jl")
include("models/cubic/vdW.jl")
include("models/cubic/RK/RK.jl")
include("models/cubic/PR/PR.jl")

include("models/Activity/Wilson/Wilson.jl")
include("models/Activity/NRTL/NRTL.jl")
include("models/Activity/UNIQUAC/UNIQUAC.jl")
include("models/Activity/UNIFAC/UNIFAC.jl")
include("models/Activity/UNIFAC/variants/ogUNIFAC.jl")
include("models/Activity/UNIFAC/variants/PSRK.jl")
include("models/Activity/UNIFAC/variants/VTPR.jl")
include("models/Activity/equations.jl")

include("models/Activity/COSMOSAC/utils.jl")
include("models/Activity/COSMOSAC/COSMOSAC02.jl")
include("models/Activity/COSMOSAC/COSMOSAC10.jl")
include("models/Activity/COSMOSAC/COSMOSACdsp.jl")

include("models/cubic/alphas/alphas.jl")
include("models/cubic/mixing/mixing.jl")
include("models/cubic/translation/translation.jl")

include("models/cubic/RK/variants/SRK.jl")
include("models/cubic/RK/variants/PSRK.jl")
include("models/cubic/PR/variants/PR78.jl")
include("models/cubic/PR/variants/VTPR.jl")
include("models/cubic/PR/variants/UMRPR.jl")

include("models/EmpiricHelmholtz/IAPWS95/IAPWS95.jl")
include("models/EmpiricHelmholtz/IAPWS95/IAPWS95Ideal.jl")
include("models/EmpiricHelmholtz/PropaneRef.jl")
include("models/EmpiricHelmholtz/LJRef/LJRef.jl")
include("models/EmpiricHelmholtz/LJRef/LJRefIdeal.jl")
include("models/EmpiricHelmholtz/MultiFluid/multifluid.jl")

include("models/LatticeFluid/SanchezLacombe/SanchezLacombe.jl")

include("models/SPUNG/SPUNG.jl")
include("models/UFTheory/UFTheory.jl")

include("models/cached/CachedEoS.jl")

end # module
