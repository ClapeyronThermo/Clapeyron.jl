abstract type KiselevModel <: CrossOverModel end

struct Kiselev2000Param{T} <: ParametricEoSParam{T}
    Tc::SingleParam{T}
    Pc::SingleParam{T}
    Vc::SingleParam{T} #TODO: add more parameters
end

@newmodelsimple Kiselev2000 KiselevModel Kiselev2000Param{T}

#TODO: actually define the model
function a_res_crossover(model::CrossOver,V,T,z,critmodel::KiselevModel,basedata)
    return a_res(model.basemodel,V,T,z)
end

function recombine_crossover!(model,critmodel::Kiselev2000,pures)
    return model
end