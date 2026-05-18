abstract type SecondVirialModel <: EoSModel end

function a_res(model::SecondVirialModel,V,T,z)
    B = second_virial_coefficient(model,T,z)
    return B/V #/a
end

lb_volume(model::SecondVirialModel,T,z) = zero(eltype(z))

function volume_impl(model::SecondVirialModel,p,T,z,phase,threaded,vol0)
    B = second_virial_coefficient(model,T,z)
    return volume_virial(B,p,T,z)
end

include("EoSVirial/EoSVirial.jl")
include("TsonopoulosVirial/TsonopoulosVirial.jl")
include("AbbottVirial/AbbottVirial.jl")
