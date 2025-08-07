"""
    SAFTgammaEMie(solvents::Array{String,1}, 
        ions::Array{String,1}; 
        idealmodel::IdealModel = BasicIdeal,
        neutralmodel::EoSModel = SAFTgammaMie,
        ionmodel::IonModel = GCMSABorn,
        RSPmodel::RSPModel = Schreckenberg,
        userlocations::Vector{String} = [],
        ideal_userlocations::Vector{String} = [],
        assoc_options::AssocOptions = AssocOptions(),
        verbose::Bool = false)

## Description
This function is used to create an SAFT-gammaE Mie model which is a combination of the SAFT-gamma Mie, MSA and Born models.

## Input parameters
### SAFT-VR Mie Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `vst`: Single Parameter (`Float64`) - Number of segments (no units)
- `S`: Single Parameter (`Float64`) - Shape factor for segment (no units)
- `sigma`: Single Parameter (`Float64`) - Segment Diameter `[Å]`
- `epsilon`: Single Parameter (`Float64`) - Reduced dispersion energy `[K]`
- `lambda_a`: Pair Parameter (`Float64`) - Atractive range parameter (no units)
- `lambda_r`: Pair Parameter (`Float64`) - Repulsive range parameter (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume
### MSA Parameters
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
1. Haslam, A.J., González-Pérez, A., Di Lecce, S., Khalit, S.H., Perdomo, F.A., Kournopoulos, S., Kohns, M., Lindeboom, T., Wehbe, M., Febra, S., Jackson, G., Adjiman, C.S. & Galind, A. (2020). Expanding the Applications of the SAFT-γ Mie Group-Contribution Equation of State: Prediction of Thermodynamic Properties and Phase Behavior of Mixtures. Journal of Chemical Engineering Data, 65(12), 5862–5890
"""
function SAFTgammaEMie(solvents,ions; 
    idealmodel = BasicIdeal,
    neutralmodel = SAFTgammaMie,
    ionmodel = GCMSABorn,
    RSPmodel = Schreckenberg,
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

export SAFTgammaEMie