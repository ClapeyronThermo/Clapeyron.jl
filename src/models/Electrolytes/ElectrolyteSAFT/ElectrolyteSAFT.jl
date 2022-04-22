abstract type ElectrolyteSAFTModel <: ElectrolyteModel end

struct ElectrolyteSAFT{T<:IdealModel,c<:SAFTModel,i<:IonModel,b} <: ElectrolyteSAFTModel
    components::Array{String,1}
    solvents::Array{String,1}
    ions::Array{String,1}
    icomponents::UnitRange{Int}
    isolvents::UnitRange{Int}
    iions::UnitRange{Int}
    stoic_coeff::Vector{Vector{Int64}}
    idealmodel::T
    puremodel::c
    ionicmodel::i
    bornmodel::b
    absolutetolerance::Float64
    references::Array{String,1}
end

@registermodel ElectrolyteSAFT
export ElectrolyteSAFT

function ElectrolyteSAFT(solvents,salts; puremodel=PCSAFT,
    idealmodel = BasicIdeal,
    ionicmodel = DH,
    RSPmodel = ConstW,
    bornmodel = nothing,
    userlocations=String[], 
    ideal_userlocations=String[],
     verbose=false)
    ion_groups = GroupParam(salts, ["Electrolytes/properties/salts.csv"]; verbose=verbose)

    ions = ion_groups.flattenedgroups
    components = deepcopy(solvents)
    append!(components,ions)
    stoichiometric_coeff = ion_groups.n_flattenedgroups
    icomponents = 1:length(components)
    isolvents = 1:length(solvents)
    iions = (length(solvents)+1):length(components)

    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    init_SAFTmodel = puremodel(components)
    init_Ionicmodel = ionicmodel(solvents,salts;RSPmodel=RSPmodel,SAFTlocations=["SAFT/"*string(puremodel)])
    if bornmodel !== nothing
        init_bornmodel = bornmodel(solvents,salts;RSPmodel=RSPmodel,SAFTlocations=["SAFT/"*string(puremodel)])
    else
        init_bornmodel = nothing
    end
    references = String[]
    model = ElectrolyteSAFT(components,solvents,ions,icomponents,isolvents,iions,stoichiometric_coeff,init_idealmodel,init_SAFTmodel,init_Ionicmodel,init_bornmodel,1e-12,references)
    return model
end

function a_res(model::ElectrolyteSAFTModel, V, T, z)
    return a_res(model.puremodel,V,T,z)+a_ion(model,V,T,z,model.ionicmodel)+a_born(model,V,T,z,model.bornmodel)
end

