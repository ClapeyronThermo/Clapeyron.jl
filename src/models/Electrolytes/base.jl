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


"""
    ESElectrolyte(solvents::Array{String,1}, 
        ions::Array{String,1}; 
        idealmodel::IdealModel = BasicIdeal,
        neutralmodel::EoSModel = pharmaPCSAFT,
        ionmodel::IonModel = DH,
        RSPmodel::RSPModel = ConstRSP,
        userlocations::Vector{String}=[],
        ideal_userlocations::Vector{String}=[],
        verbose::Bool=false)

## Description
This function provides the necessary framework to create an electrolyte model by combining ideal, neutral and ion models:
```julia
model = ESElectrolyte(["water"],["sodium","chloride"];
            idealmodel = BasicIdeal,
            neutralmodel = pharmaPCSAFT,
            ionmodel = DH,
            RSPmodel = ConstRSP)  
```
Any of the available models in Clapeyron can be combined in the above. Note that neutral (solvent) species and ions are defined separately. Within Clapeyron, we will only support ion-based electrolyte models; as such, any salt-based approach (i.e. where the salt is treated as a separate species) will not be supported. 
"""
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

function lb_volume(model::ESElectrolyteModel,z)
    return lb_volume(model.neutralmodel,z)
end

function lb_volume(model::ESElectrolyteModel,T,z)
    return lb_volume(model.neutralmodel,T,z)
end

function x0_volume_liquid(model::ESElectrolyteModel,T,z)
    return x0_volume_liquid(model.neutralmodel, T, z)*1.15
end

function mw(model::ElectrolyteModel)
    return mw(model.neutralmodel)
end

function p_scale(model::ElectrolyteModel,z)
    return p_scale(model.neutralmodel,z)
end

function T_scale(model::ElectrolyteModel,z)
    return T_scale(model.neutralmodel,z)
end