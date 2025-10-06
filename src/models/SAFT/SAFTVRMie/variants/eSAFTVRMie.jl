abstract type eSAFTVRMieModel <: ESElectrolyteModel end

struct eSAFTVRMie{T<:IdealModel,c<:EoSModel,i<:IonModel} <: eSAFTVRMieModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    charge::Vector{Int64}
    idealmodel::T
    neutralmodel::c
    ionmodel::i
    references::Array{String,1}
end

"""
    eSAFTVRMie(solvents::Array{String,1}, 
        ions::Array{String,1}; 
        idealmodel::IdealModel = BasicIdeal,
        neutralmodel::EoSModel = SAFTVRMie15,
        ionmodel::IonModel = DHBorn,
        RSPmodel::RSPModel = ZuoFurst,
        charges = String[], 
        ideal_userlocations = String[],
        neutralmodel_userlocations = String[],
        ionmodel_userlocations = String[],
        RSPmodel_userlocations = String[],
        assoc_options::AssocOptions = AssocOptions(),
        verbose::Bool = false,
        reference_state = nothing)

## Description
This function is used to create an eSAFTVRMie model which is a combination of the SAFTVR-Mie, Debye-Hückel and Born models.

## Input parameters
### SAFT-VR Mie Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Single Parameter (`Float64`) - Segment Diameter `[Å]`
- `epsilon`: Single Parameter (`Float64`) - Reduced dispersion energy `[K]`
- `lambda_a`: Pair Parameter (`Float64`) - Atractive range parameter (no units)
- `lambda_r`: Pair Parameter (`Float64`) - Repulsive range parameter (no units)
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Parameter (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m³]`
### Debye-Hückel Parameters
- `sigma`: Single Parameter (`Float64`) - Diameter of closest approach `[m]`
- `charge`: Single Parameter (`Float64`) - Charge `[-]`
### Born Parameters
- `sigma_born`: Single Parameter (`Float64`) - Born Diameter `[m]`
- `charge`: Single Parameter (`Float64`) - Charge `[-]`

## Input models
- `idealmodel`: Ideal Model
- `neutralmodel`: Neutral EoS Model
- `ionmodel`: Ion Model

## References
1. Selam, M., Economou, I., Castier, M. (2018). A thermodynamic model for strong aqueous electrolytes based on the eSAFT-VR Mie equation of state. Fluid Phase Equilibria, 464, 47-63.
"""
function eSAFTVRMie(solvents,ions; 
    idealmodel = BasicIdeal,
    neutralmodel = SAFTVRMie15,
    ionmodel = DHBorn,
    RSPmodel = ZuoFurst,
    charges = String[], 
    ideal_userlocations = String[],
    neutralmodel_userlocations = String[],
    ionmodel_userlocations = String[],
    RSPmodel_userlocations = String[],
    assoc_options = AssocOptions(),
    reference_state = nothing,
    verbose = false)

    return ESElectrolyte(solvents,ions;
    idealmodel,neutralmodel,ionmodel,RSPmodel,
    charges,ideal_userlocations,neutralmodel_userlocations,ionmodel_userlocations,RSPmodel_userlocations,assoc_options,reference_state,verbose)
end

export eSAFTVRMie
