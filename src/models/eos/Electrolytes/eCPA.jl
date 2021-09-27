abstract type eCPAModel <: ElectrolyteModel end

struct eCPA{T<:IdealModel,c<:SAFTModel,i<:IonModel,b} <: eCPAModel
    components::Array{String,1}
    solvents::Array{String,1}
    salts::Array{String,1}
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

@registermodel eCPA
export eCPA

function eCPA(solvents,salts; 
    puremodel=CPA,
    cubicmodel=RK, 
    alpha=CPAAlpha, 
    mixing=HVRule,
    activity=eCPANRTL,
    translation=NoTranslation, 
    idealmodel = BasicIdeal,
    ionicmodel = DH,
    RSPmodel = MM1,
    bornmodel = Born,
    userlocations=String[], 
    ideal_userlocations=String[],
    alpha_userlocations=String[],
    activity_userlocations=String[],
    mixing_userlocations=String[],
    translation_userlocations=String[],
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
    init_SAFTmodel =    puremodel(components;
                        idealmodel=idealmodel, 
                        cubicmodel=cubicmodel, 
                        alpha=alpha, 
                        mixing=mixing,
                        activity=activity,
                        translation=translation, 
                        userlocations=userlocations, 
                        ideal_userlocations=ideal_userlocations, 
                        alpha_userlocations=alpha_userlocations,
                        activity_userlocations=activity_userlocations,
                        mixing_userlocations=mixing_userlocations,
                        translation_userlocations=translation_userlocations)
    init_Ionicmodel = ionicmodel(solvents,salts;RSPmodel=RSPmodel,SAFTlocations=["SAFT/"*string(puremodel)])
    if bornmodel !== nothing
        init_bornmodel = bornmodel(solvents,salts;RSPmodel=RSPmodel,SAFTlocations=["SAFT/"*string(puremodel)])
    else
        init_bornmodel = nothing
    end
    references = String[]
    model = eCPA(components,solvents,salts,ions,icomponents,isolvents,iions,stoichiometric_coeff,init_idealmodel,init_SAFTmodel,init_Ionicmodel,init_bornmodel,1e-12,references)
    return model
end

function a_res(model::eCPAModel, V, T, z)
    if model.bornmodel !== nothing
        return a_res(model.puremodel,V,T,z)+a_ion(model.ionicmodel,V,T,z)+a_born(model.bornmodel,V,T,z)
    else
        return a_res(model.puremodel,V,T,z)+a_ion(model.ionicmodel,V,T,z)
    end
end