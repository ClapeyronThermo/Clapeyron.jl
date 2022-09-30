abstract type VirialModel <: EoSModel end

function a_res(model::VirialModel,V,T,z)
    B = second_virial_coefficient(model,T,z)
    return B/V #/a
end

lb_volume(model::VirialModel,z = SA[1.0]) = zero(eltype(z))

function volume_impl(model::VirialModel,p,T,z,phase,threaded)
    B = second_virial_coefficient(model,T,z)
    return volume_virial(B,p,T,z)
end

function ∂f∂V(model::VirialModel,V,T,z=SA[1.0])
    B = second_virial_coefficient(model,T,z)
    return -R̄*T*(1 + B/V)/V
end

include("EoSVirial/EoSVirial.jl")
include("TsonopoulosVirial/TsonopoulosVirial.jl")
include("AbbottVirial/AbbottVirial.jl")
