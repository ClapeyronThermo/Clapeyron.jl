struct EoSVirial2{I,T} <: SecondVirialModel
    idealmodel::I
    model::T
end

@registermodel EoSVirial2

function EoSVirial2(model;idealmodel = idealmodel(model))
    return EoSVirial2{typeof(idealmodel),typeof(model)}(idealmodel,model)
end

second_virial_coefficient_impl(model::EoSVirial2,T,z=SA[1.0]) =  second_virial_coefficient_impl(model.model,T,z)

export EoSVirial2