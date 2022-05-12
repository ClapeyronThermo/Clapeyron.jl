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
using Unitful
import LogExpFunctions

include("solvers/Solvers.jl")
using .Solvers
using .Solvers: log, sqrt
∂Tag = Solvers.∂Tag

include("utils/fractions.jl")
import .Fractions
using .Fractions: FractionVector

#Gas constant, Boltzmann Constant
include("base/constants.jl") 

#The Base of Clapeyron: EoSModel and eos(model,V,T,z)
include("base/EoSModel.jl")

#show(model<:EoSModel)
include("base/eosshow.jl")

#EoSParam, ClapeyronParam, All Params
include("database/ClapeyronParam.jl")

#Combining Rules for Single and Pair Params.
include("database/params/combiningrules.jl")

using CSV, Tables
#getparams options
include("database/ParamOptions.jl") 
#getparams definition
include("database/database.jl")
#transform Tables.jl tables to Clapeyron csv files
include("database/UserReader.jl")

#macros, used for defining models
include("utils/macros.jl")

#index reduction
include("utils/index_reduction.jl")

#splitting models, useful for methods.
include("utils/split_model.jl")

#Clapeyron methods (AD, property solvers, etc)
include("methods/methods.jl")

#=
the dependency chain is the following:

base --> database(params)  -|-> split_model --> methods -|-> models                     
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

#Basic utility EoS
include("models/utility/SpecialComp.jl")
include("models/utility/EoSVectorParam.jl")
include("models/utility/ZeroResidual.jl")

#softSAFT2016 uses LJRef. softSAFT uses x0_sat_pure with LJ correlations (from LJRef)
include("models/EmpiricHelmholtz/IAPWS95/IAPWS95.jl")
include("models/EmpiricHelmholtz/IAPWS95/IAPWS95Ideal.jl")
include("models/EmpiricHelmholtz/PropaneRef.jl")
include("models/EmpiricHelmholtz/LJRef/LJRef.jl")
include("models/EmpiricHelmholtz/LJRef/LJRefIdeal.jl")
include("models/EmpiricHelmholtz/MultiFluid/multifluid.jl")

include("models/cubic/equations.jl")
include("models/cubic/vdW.jl")
include("models/cubic/RK/RK.jl")
include("models/cubic/PR/PR.jl")

include("models/SAFT/PCSAFT/PCSAFT.jl")
include("models/SAFT/PCSAFT/variants/sPCSAFT.jl")
include("models/SAFT/PCSAFT/variants/PharmaPCSAFT.jl")
include("models/SAFT/ogSAFT/ogSAFT.jl")
include("models/SAFT/CPA/CPA.jl")
include("models/SAFT/CPA/variants/sCPA.jl")
include("models/SAFT/SAFTVRSW/SAFTVRSW.jl")
include("models/SAFT/LJSAFT/LJSAFT.jl")
include("models/SAFT/softSAFT/softSAFT.jl")
include("models/SAFT/softSAFT/variants/softSAFT2016.jl")

include("models/SAFT/SAFTVRMie/SAFTVRMie.jl")
include("models/SAFT/SAFTVRMie/variants/SAFTVRQMie.jl")
include("models/SAFT/SAFTgammaMie/SAFTgammaMie.jl")
include("models/SAFT/CKSAFT/CKSAFT.jl")
include("models/SAFT/CKSAFT/variants/sCKSAFT.jl")
include("models/SAFT/BACKSAFT/BACKSAFT.jl")
include("models/SAFT/equations.jl")
include("models/SAFT/association.jl")


include("models/Activity/Wilson/Wilson.jl")
include("models/Activity/NRTL/NRTL.jl")
include("models/Activity/UNIQUAC/UNIQUAC.jl")
include("models/Activity/UNIFAC/utils.jl")
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
include("models/cubic/PR/variants/QCPR.jl")
include("models/cubic/PR/variants/EPPR78.jl")

include("models/LatticeFluid/SanchezLacombe/SanchezLacombe.jl")

include("models/SPUNG/SPUNG.jl")
include("models/UFTheory/UFTheory.jl")

include("utils/misc.jl")
end # module
