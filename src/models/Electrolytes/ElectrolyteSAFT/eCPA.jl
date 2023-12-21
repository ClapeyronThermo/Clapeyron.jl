abstract type eCPAModel <: ElectrolyteModel end

struct eCPA{T<:IdealModel,c,i<:IonModel,b,ϵ<:RSPModel} <: eCPAModel
    components::Array{String,1}
    solvents::Array{String,1}
    salts::Array{String,1}
    ions::Array{String,1}
    isolvents::UnitRange{Int}
    iions::UnitRange{Int}
    stoic_coeff::Vector{Vector{Int64}}
    idealmodel::T
    puremodel::c
    ionicmodel::i
    bornmodel::b
    rspmodel::ϵ
    references::Array{String,1}
end

@registermodel eCPA
export eCPA

function eCPA(solvents,salts; 
    puremodel=sCPA,
    cubicmodel=RK, 
    alpha=sCPAAlpha, 
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
    
    ion_groups = GroupParam(salts, ["Electrolytes/properties/salts.csv"]; verbose=verbose)
    salts = ion_groups.components
    ions = ion_groups.flattenedgroups
    components = deepcopy(solvents)
    append!(components,ions)
    stoichiometric_coeff = ion_groups.n_flattenedgroups
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
    init_Ionicmodel = ionicmodel(solvents,salts;RSPmodel=nothing,SAFTlocations=["SAFT/"*string(puremodel)])
    
    if bornmodel !== nothing
        init_bornmodel = bornmodel(solvents,salts;RSPmodel=nothing,SAFTlocations=["SAFT/"*string(puremodel)])
    end
    init_RSPmodel = RSPmodel(solvents,salts)
    references = String[]
    model = eCPA(components,solvents,salts,ions,isolvents,iions,stoichiometric_coeff,init_idealmodel,init_SAFTmodel,init_Ionicmodel,init_bornmodel,init_RSPmodel,references)
    return model
end

function data(model::eCPAModel, V, T, z)
    data_cpa = data(model.puremodel,V,T,z)
    data_rsp = dielectric_constant(model.rspmodel,V,T,z,model.puremodel)
    return (data_cpa,data_rsp)
end


function a_res(model::eCPAModel, V, T, z,_data = @f(data))
    (data_cpa,data_rsp) = _data
    n,ā,b̄,c̄ = data_cpa
    return a_res(model.puremodel,V,T,z,data_cpa) + a_res(model.ionicmodel,V+c̄*n,T,z,data_rsp) + a_res(model.bornmodel,V+c̄*n,T,z,data_rsp)

end

#TODO: fix for salts
lb_volume(model::eCPAModel,z = SA[1.0]) = lb_volume(model.puremodel,z)
p_scale(model::eCPAModel,z = SA[1.0]) = p_scale(model.puremodel,z)
T_scale(model::eCPAModel,z = SA[1.0]) = T_scale(model.puremodel,z)