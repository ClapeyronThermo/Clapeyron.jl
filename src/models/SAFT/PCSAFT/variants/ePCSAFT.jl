abstract type ePCSAFTModel <: ESElectrolyteModel end

struct ePCSAFT{T<:IdealModel,c<:EoSModel,i<:IonModel} <: ePCSAFTModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    charge::Vector{Int64}
    idealmodel::T
    neutralmodel::c
    ionmodel::i
    references::Array{String,1}
end

"""
    ePCSAFT(solvents::Array{String,1}, 
        ions::Array{String,1}; 
        idealmodel::IdealModel = BasicIdeal,
        neutralmodel::EoSModel = pharmaPCSAFT,
        ionmodel::IonModel = DH,
        RSPmodel::RSPModel = ConstRSP,
        userlocations::Vector{String}=[],
        ideal_userlocations::Vector{String}=[],
        assoc_options::AssocOptions = AssocOptions(),
        verbose::Bool=false)

## Description
This function is used to create an ePCSAFT model which is a combination of the PC-SAFT and Debye-Hückel model. It is based on the ePC-SAFT Revised variant.

## Input parameters
### PC-SAFT Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Single Parameter (`Float64`) - Segment Diameter [`A°`]
- `epsilon`: Single Parameter (`Float64`) - Reduced dispersion energy  `[K]`
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Paramater (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m^3]`
### Debye-Hückel Parameters
- `sigma`: Single Parameter (`Float64`) - Diameter of closest approach `[m]`
- `charge`: Single Parameter (`Float64`) - Charge `[-]`

## Input models
- `idealmodel`: Ideal Model
- `neutralmodel`: Neutral EoS Model
- `ionmodel`: Ion Model

## References
1. Held, C., Reschke, T., Mohammad, S., Luza, A., Sadowski, G. (2014). ePC-SAFT Revised. Chemical Engineering Research and Design, 92(12), 2884-2897.
"""
function ePCSAFT(solvents,ions; 
    idealmodel = BasicIdeal,
    neutralmodel = pharmaPCSAFT,
    ionmodel = DH,
    RSPmodel = ConstRSP,
    userlocations=String[], 
    ideal_userlocations=String[],
    assoc_options = AssocOptions(),
     verbose=false)
    components = deepcopy(ions)
    prepend!(components,solvents)

    params = getparams(format_components(components), ["Electrolytes/properties/charges.csv"]; userlocations=userlocations, verbose=verbose)
    _charge = params["charge"]
    charge = _charge.values

    icomponents = 1:length(components)

    neutral_path = DB_PATH.*["/SAFT/PCSAFT","/SAFT/PCSAFT/ePCSAFT","/SAFT/PCSAFT/pharmaPCSAFT"]

    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    init_neutralmodel = neutralmodel(components;userlocations=append!(userlocations,neutral_path),verbose=verbose,assoc_options=assoc_options)
    init_ionmodel = ionmodel(solvents,ions;RSPmodel=RSPmodel,userlocations=append!(userlocations,neutral_path),verbose=verbose)


    for i in ions
        init_neutralmodel.params.epsilon[i] = 0. #pure ion has ϵi 
        for j in ions
            if sign(_charge[i]) == sign(_charge[j]) #cation-cation and anion-anion interactions are neglected.
                init_neutralmodel.params.epsilon[i,j] = 0.
            end
        end
    end



    references = ["10.1016/j.cherd.2014.05.017"]
    components = format_components(components)
    model = ePCSAFT(components,icomponents,charge,init_idealmodel,init_neutralmodel,init_ionmodel,references)
    return model
end

function a_res(model::ePCSAFTModel, V, T, z)
    data_pcsaft = data(model.neutralmodel,V,T,z)
    data_ion = data(model.ionmodel,V,T,z)
    data_ion = (data_ion[1],data_pcsaft[1])

    return a_res(model.neutralmodel,V,T,z,data_pcsaft)+a_res(model.ionmodel,V,T,z,data_ion)
end

export ePCSAFT