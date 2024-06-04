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

function eSAFTVRMie(solvents,ions; 
    idealmodel = BasicIdeal,
    neutralmodel = SAFTVRMie15,
    ionmodel = DHBorn,
    RSPmodel = ZuoFirst,
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

    neutral_path = DB_PATH.*["/SAFT/SAFTVRMie/MieKernel"]

    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    init_neutralmodel = neutralmodel(components;userlocations=userlocations,verbose=verbose,assoc_options=assoc_options)
    init_ionmodel = ionmodel(solvents,ions;RSPmodel=RSPmodel,userlocations=append!(userlocations,neutral_path),verbose=verbose)

    references = String[]
    components = format_components(components)
    model = eSAFTVRMie(components,icomponents,charge,init_idealmodel,init_neutralmodel,init_ionmodel,references)
    return model
end

function a_res(model::eSAFTVRMieModel, V, T, z)
    data_saft = data(model.neutralmodel,V,T,z)
    data_ion = data(model.ionmodel,V,T,z)

    return a_res(model.neutralmodel,V,T,z,data_saft)+a_res(model.ionmodel,V,T,z,data_ion)
end

export eSAFTVRMie