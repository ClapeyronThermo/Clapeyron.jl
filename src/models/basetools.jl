abstract type EoSModel end
abstract type LatticeFluidModel <: EoSModel end
abstract type SAFTModel <: LatticeFluidModel end
abstract type GCSAFTModel <: SAFTModel end
abstract type CubicModel <: EoSModel end
abstract type ABCubicModel <: CubicModel end #cubics that have an exact polynomial form to solve Z roots, this excludes CPA
abstract type IdealModel <: EoSModel end
abstract type EmpiricHelmholtzModel <: EoSModel end
abstract type EoSParam end
abstract type SatPureAproximation <: EoSModel end
abstract type AlphaModel <:EoSModel end
export EoSModel
export SAFTModel,CubicModel,EmpiricHelmholtzModel
export IdealModel
export EoSParam
export AlphaModel

