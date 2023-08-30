abstract type LatticeFluidModel <: EoSModel end
abstract type SAFTModel <: LatticeFluidModel end
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

"""
    CompositeModel(components;
    gas = BasicIdeal,
    liquid = RackettLiquid,
    saturation = LeeKeslerSat,
    gas_userlocations = String[],
    liquid_userlocations = String[],
    saturation_userlocations = String[]

Composite Model. it is not consistent, but it can hold different correlations that
are faster than a volume or saturation pressure iteration.

"""
struct CompositeModel{𝕍,𝕃,𝕊,𝕃𝕍,𝕃𝕊} <: EoSModel
    components::Vector{String}
    gas::𝕍
    liquid::𝕃
    solid::𝕊
    saturation::𝕃𝕍
    melting::𝕃𝕊
end
