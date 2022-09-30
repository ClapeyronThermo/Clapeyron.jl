struct EoSVirial{I,T} <: VirialModel
    idealmodel::I
    model::T
end

@registermodel EoSVirial

function EoSVirial(model;idealmodel = BasicIdeal())
    return EoSVirial{typeof(idealmodel),typeof(model)}(idealmodel,model)
end

second_virial_coefficient_impl(model::EoSVirial,T,z=SA[1.0]) =  second_virial_coefficient_impl(model.model,T,z)

export EoSVirial