abstract type LatticeFluidModel <: EoSModel end
abstract type AssociationModel <: LatticeFluidModel end
abstract type CPAModel <: AssociationModel end
abstract type SAFTModel <: AssociationModel end
abstract type GCSAFTModel <: SAFTModel end

abstract type CubicModel <: EoSModel end
abstract type MixingRule <:EoSModel end #mixing rules for cubics (and empiric Helmholtz models)
abstract type ABCubicModel <: CubicModel end #cubics that have an exact polynomial form to solve Z roots, this excludes CPA
abstract type ABCCubicModel <: ABCubicModel end #cubics with a "c" variable parameter, like Clausius and patel-teja
abstract type ActivityModel <: EoSModel end
abstract type IdealModel <: EoSModel end
abstract type EmpiricHelmholtzModel <: EoSModel end
abstract type SatPureAproximation <: EoSModel end

export SAFTModel,CubicModel,EmpiricHelmholtzModel
export IdealModel
export AlphaModel

#check models/CompositeModel/CompositeModel.jl
struct CompositeModel{𝔽,𝕊} <: EoSModel
    components::Vector{String}
    fluid::𝔽
    solid::𝕊
end
 