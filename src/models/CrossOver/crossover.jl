
abstract type CrossOverModel <: EoSModel end

struct CritParam{T} <: ParametricEoSParam{T}
    Tc0::SingleParam{T}
    Pc0::SingleParam{T}
    Vc0::SingleParam{T}
end

struct CrossOver{M,C,T} <: EoSModel
    components::Vector{String}
    params::CritParam{T} #critical point of the base model.
    basemodel::M
    critmodel::C
end


"""
CrossOver(model::EoSModel,crossover::CrossOverModel;verbose = false)
CrossOver(model::EoSModel;crossover = nothing,crossover_userlocations = String[],verbose = false)

Cross-over model. Performs critical point renormalization to exactly match the critical point of a substance.
"""
function CrossOver(model::EoSModel;crossover = nothing,crossover_userlocations = String[],verbose = false)
    components = model.components
    init_crossover = init_model(crossover,components,userlocations,verbose)
    return CrossOver(model,critmodel;verbose)
end

function CrossOver(model::EoSModel,critmodel::CrossOverModel;verbose = false)
    components = model.components
    Tc0 = SingleParam("Tc0",model.components)
    Pc0 = SingleParam("Pc0",model.components)
    Vc0 = SingleParam("Vc0",model.components)
    params = CritParam(Tc0,Pc0,Vc0)
    crossover_model = CrossOver(components,params,model,critmodel)
    recombine_crossover!(crossover_model,critmodel)
    set_reference_state!(crossover_model,model.basemodel.reference_state;verbose)
    return crossover_model
end

function a_res(model::CrossOver,V,T,z)
    return a_res_crossover(model,V,T,z,model.crossover)
end

function recombine_impl!(model::CrossOver)
    basemodel = model.basemodel
    recombine!(basemodel)
    pures = split_pure_model(basemodel)
    crits = crit_pure.(pures)
    Tc = first.(crits)
    Pc = getindex.(crits,2)
    Vc = last.(crits)
    model.params.Pc0 .= Pc
    model.params.Tc0 .= Tc
    model.params.Vc0 .= Vc
    recombine_crossover!(model,critmodel,pures)
end

recombine_crossover!(model,critmodel) = recombine_crossover!(model,critmodel,split_pure_model(model.basemodel)) 

idealmodel(model::CrossOver) = idealmodel(model.basemodel)
lb_volume(model::CrossOver,T,z) = lb_volume(model.basemodel,T,z)
T_scale(model::CrossOver,z) = T_scale(model.basemodel,z)
p_scale(model::CrossOver,z) = p_scale(model.basemodel,z)
Rgas(model::CrossOver) = Rgas(model.basemodel)
molecular_weight(model::CrossOver,z) = molecular_weight(model.basemodel,z)
Base.eltype(x::CrossOver) = Base.promote_eltype(x.basemodel,x.params)
reference_state(model::CrossOver) = reference_state(model.basemodel)


function Base.show(io::IO,mime::MIME"text/plain",model::CrossOver)
    print(io,"Critical cross-over Model")
    length(model) == 1 && print(io, " with 1 component:")
    length(model) > 1 && print(io, " with ", length(model), " components:")
    print(io,'\n'," Base Model: ",model.basemodel)
    print(io,'\n'," Liquid Model: ",model.critmodel)
end