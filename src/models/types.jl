abstract type LatticeFluidModel <: EoSModel end
abstract type AssociationModel <: LatticeFluidModel end
abstract type CPAModel <: AssociationModel end
abstract type SAFTModel <: AssociationModel end
abstract type GCSAFTModel <: SAFTModel end
abstract type GibbsBasedModel <: EoSModel end
abstract type CubicModel <: EoSModel end #cubics that have an exact polynomial form to solve Z roots, this excludes CPA
abstract type DeltaCubicModel <: CubicModel end #cubics with a a/((v-Δ1*b)*(v-Δ2*b)) term
abstract type ABCubicModel <: DeltaCubicModel end #cubics where the critical point is just a function of Δ1,Δ2 (enforced via mixing rules)
abstract type ABCCubicModel <: DeltaCubicModel end #cubics with a "c" variable parameter, like Clausius and patel-teja

abstract type MixingRule <:EoSModel end #mixing rules for cubics (and empiric Helmholtz models)
abstract type ActivityModel <: EoSModel end
abstract type IdealModel <: EoSModel end
abstract type EmpiricHelmholtzModel <: EoSModel end
abstract type SatPureAproximation <: EoSModel end
abstract type AlphaModel <:EoSModel end
abstract type ElectrolyteModel <: EoSModel end
abstract type IonModel <: ElectrolyteModel end
abstract type RSPModel <: ElectrolyteModel end

export SAFTModel,CubicModel,EmpiricHelmholtzModel
export IdealModel
export AlphaModel

#check models/CompositeModel/CompositeModel.jl
struct CompositeModel{𝔽,𝕊} <: EoSModel
    components::Vector{String}
    fluid::𝔽
    solid::𝕊
    mapping::Union{Vector{Pair{Vector{Tuple{String,Int64}},Tuple{String,Int64}}},Nothing}
end
