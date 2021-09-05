abstract type ElectrolyteModel <: EoSModel end
abstract type ElectrolyteSAFTModel <: ElectrolyteModel end

struct ElectrolyteSAFT{T<:IdealModel,c<:SAFTModel,i<:IonModel} <: ElectrolyteSAFTModel
    components::Array{String,1}
    solvent::Array{String,1}
    ions::Array{String,1}
    icomponents::UnitRange{Int}
    isolvents::UnitRange{Int}
    iions::UnitRange{Int}
    stoic_coeff::Vector{Vector{Int64}}
    idealmodel::T
    puremodel::c
    Ionicmodel::i
    absolutetolerance::Float64
    references::Array{String,1}
end

@registermodel ElectrolyteSAFT
export ElectrolyteSAFT

function ElectrolyteSAFT(solvents,salts; puremodel=PCSAFT,
    idealmodel = BasicIdeal,
    ionicmodel = DH,
    RSPmodel = ConstW,
    userlocations=String[], 
    ideal_userlocations=String[],
     verbose=false)
    ion_groups = GroupParam(salts, ["Electrolytes/salts.csv"]; verbose=verbose)

    ions = ion_groups.flattenedgroups
    components = deepcopy(solvents)
    append!(components,ions)
    stoichiometric_coeff = ion_groups.n_flattenedgroups
    icomponents = 1:length(components)
    isolvents = 1:length(solvents)
    iions = (length(solvents)+1):length(components)

    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    init_SAFTmodel = puremodel(components)
    init_Ionicmodel = ionicmodel(components;RSPmodel=RSPmodel,SAFTlocations=["SAFT/"*string(puremodel)])
    references = String[]
    model = ElectrolyteSAFT(components,solvents,ions,icomponents,isolvents,iions,stoichiometric_coeff,init_idealmodel,init_SAFTmodel,init_Ionicmodel,1e-12,references)
    return model
end

function a_res(model::ElectrolyteSAFTModel, V, T, z)
    return a_res(model.puremodel,V,T,z)+a_ion(model.Ionicmodel,V,T,z)
end