abstract type VirialModel <: EoSModel end

struct SecondVirial{I,T} <: VirialModel
    idealmodel::I
    model::T
end

@registermodel SecondVirial

function SecondVirial(model;idealmodel = BasicIdeal())
    return SecondVirial{typeof(idealmodel),typeof(model)}(idealmodel,model)
end

function a_res(model::VirialModel,V,T,z)
    B = second_virial_coefficient(model,T,z)
    return B/V #/a
end

lb_volume(model::VirialModel,z = SA[1.0]) = zero(eltype(z))

second_virial_coefficient(model::SecondVirial,T,z=SA[1.0]) =  second_virial_coefficient(model.model,T,z)


function volume_impl(model::VirialModel,p,T,z,phase,threaded)
    B = second_virial_coefficient(model,T,z)
    return volume_virial(B,p,T,z)
    #todo, stability check
end

function pressure(model::VirialModel,V,T,z=SA[1.0])
    B = second_virial_coefficient(model,T,z)
    return R̄*T*(1 + B/V)/V
end

include("TsonopoulosVirial/TsonopoulosVirial.jl")
include("AbbottVirial/AbbottVirial.jl")
