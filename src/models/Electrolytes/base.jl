abstract type ESElectrolyteModel <: ElectrolyteModel end

struct ESElectrolyte{T<:IdealModel,c<:EoSModel,i<:IonModel} <: ESElectrolyteModel
    components::Array{String,1}
    neutral::Array{String,1}
    ions::Array{String,1}
    icomponents::UnitRange{Int}
    ineutral::UnitRange{Int}
    iions::UnitRange{Int}
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

    icomponents = 1:length(components)
    isolvents = 1:length(solvents)
    iions = (length(solvents)+1):(length(solvents)+length(ions))

    neutral_path = [DB_PATH*"/"*default_locations(neutralmodel)[1]]

    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    init_neutralmodel = neutralmodel(components;userlocations=userlocations,verbose=verbose)
    init_ionmodel = ionmodel(solvents,ions;RSPmodel=RSPmodel,userlocations=append!(userlocations,neutral_path),verbose=verbose)

    references = String[]
    model = ESElectrolyte(components,solvents,ions,icomponents,isolvents,iions,init_idealmodel,init_neutralmodel,init_ionmodel,references)
    return model
end

export ESElectrolyte

function a_res(model::ESElectrolyteModel, V, T, z)
    return a_res(model.neutralmodel,V,T,z)+a_res(model.ionmodel,V,T,z)
end