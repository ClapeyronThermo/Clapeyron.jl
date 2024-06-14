"""
    SAFTgammaEMie(solvents::Array{String,1}, 
        ions::Array{String,1}; 
        idealmodel::IdealModel = BasicIdeal,
        neutralmodel::EoSModel = SAFTgammaMie,
        ionmodel::IonModel = GCMSABorn,
        RSPmodel::RSPModel = Schreckenberg,
        userlocations::Vector{String}=[],
        ideal_userlocations::Vector{String}=[],
        assoc_options::AssocOptions = AssocOptions(),
        verbose::Bool=false)

## Description
This function is used to create an SAFT-gammaE Mie model which is a combination of the SAFT-gamma Mie, MSA and Born models.

## Input parameters
### SAFT-VR Mie Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `vst`: Single Parameter (`Float64`) - Number of segments (no units)
- `S`: Single Parameter (`Float64`) - Shape factor for segment (no units)
- `sigma`: Single Parameter (`Float64`) - Segment Diameter [`A°`]
- `epsilon`: Single Parameter (`Float64`) - Reduced dispersion energy  `[K]`
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
    userlocations=String[], 
    ideal_userlocations=String[],
     verbose=false)
    components = deepcopy(ions)
    prepend!(components,solvents)

    params = getparams(format_components(components), ["Electrolytes/properties/charges.csv"]; userlocations=userlocations, verbose=verbose)
    charge = params["charge"].values

    icomponents = 1:length(components)

    neutral_path = [DB_PATH*"/SAFT/SAFTgammaMie"]

    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    init_neutralmodel = neutralmodel(components;userlocations=userlocations,verbose=verbose)
    init_ionmodel = ionmodel(solvents,ions;RSPmodel=RSPmodel,userlocations=append!(userlocations,neutral_path),verbose=verbose)

    references = String[]
    components = format_components(components)
    model = ESElectrolyte(components,icomponents,charge,init_idealmodel,init_neutralmodel,init_ionmodel,references)
    return model
end

export SAFTgammaEMie