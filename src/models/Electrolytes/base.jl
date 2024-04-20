abstract type ESElectrolyteModel <: ElectrolyteModel end

struct ESElectrolyte{T<:IdealModel,c<:EoSModel,i<:IonModel} <: ESElectrolyteModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    charge::Vector{Int64}
    idealmodel::T
    neutralmodel::c
    ionmodel::i
    references::Array{String,1}
end

function ESElectrolyte(solvents,ions; 
    idealmodel = BasicIdeal,
    neutralmodel = pharmaPCSAFT,
    ionmodel = DH,
    RSPmodel = ConstRSP,
    userlocations=String[], 
    ideal_userlocations=String[],
     verbose=false)
    components = deepcopy(ions)
    prepend!(components,solvents)

    params = getparams(components, ["Electrolytes/properties/charges.csv"]; userlocations=userlocations, verbose=verbose)
    charge = params["charge"].values

    icomponents = 1:length(components)

    neutral_path = [DB_PATH*"/"*default_locations(neutralmodel)[1]]

    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    init_neutralmodel = neutralmodel(components;userlocations=userlocations,verbose=verbose)
    init_ionmodel = ionmodel(solvents,ions;RSPmodel=RSPmodel,userlocations=append!(userlocations,neutral_path),verbose=verbose)

    references = String[]
    model = ESElectrolyte(components,icomponents,charge,init_idealmodel,init_neutralmodel,init_ionmodel,references)
    return model
end

export ESElectrolyte

function a_res(model::ESElectrolyteModel, V, T, z)
    return a_res(model.neutralmodel,V,T,z)+a_res(model.ionmodel,V,T,z)
end

function lb_volume(model::ESElectrolyteModel,z=SA[1.])
    return lb_volume(model.neutralmodel,z)
end

function x0_volume_liquid(model::ESElectrolyteModel, T,z=SA[1.])
    return x0_volume_liquid(model.neutralmodel, T,z)*1.15
end

function mw(model::ElectrolyteModel)
    return mw(model.neutralmodel)
end

function p_scale(model::ElectrolyteModel,z=SA[1.])
    return p_scale(model.neutralmodel,z)
end

function T_scale(model::ElectrolyteModel,z=SA[1.])
    return T_scale(model.neutralmodel,z)
end