abstract type ElectrolyteModel <: EoSModel end
abstract type ElectrolyteSAFTModel <: ElectrolyteModel end

struct ElectrolyteSAFT{c<:SAFTModel,i<:IonModel} <: ElectrolyteSAFTModel
    components::Array{String,1}
    solvent::Array{String,1}
    ions::Array{String,1}
    icomponents::UnitRange{Int}
    isolvents::UnitRange{Int}
    iions::UnitRange{Int}
    SAFTmodel::c
    Ionicmodel::i
    absolutetolerance::Float64
    references::Array{String,1}
end

@registermodel ElectrolyteSAFT
export ElectrolyteSAFT

function ElectrolyteSAFT(solvents,ions; SAFTmodel=PCSAFT,
    ionicmodel = DH,
    RSPmodel = ConstW,
    userlocations=String[], 
     verbose=false)
     components = append!(solvents,ions)
    icomponents = 1:length(components)
    isolvents = 1:sum(components .== solvents)
    iions = (sum(components .== solvents)+1):length(components)

    init_SAFTmodel = SAFTmodel(components)
    init_Ionicmodel = ionicmodel(components;RSPmodel=RSPmodel,SAFTlocations=["SAFT/"*string(SAFTmodel)])
    references = String[]
    model = ElectrolyteSAFT(components,solvents,ions,icomponents,isolvents,iions,init_SAFTmodel,init_Ionicmodel,1e-12,references)
    return model
end

function a_res(model::ElectrolyteSAFTModel, V, T, z)
    return a_res(model.SAFTmodel,V,T,z)+a_ion(model.Ionicmodel,V,T,z)
end