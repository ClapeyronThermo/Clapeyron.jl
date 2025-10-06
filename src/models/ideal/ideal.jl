
a_res(model::IdealModel,V,T,z) = zero(Base.promote_eltype(model,V,T,z))

function volume_impl(model::IdealModel,p,T,z,phase,threaded,vol0)
    return sum(z)*R̄*T/p
end

lb_volume(model::IdealModel,z) = zero(eltype(z))

idealmodel(model::IdealModel) = model

@newmodelsingleton ZeroIdeal IdealModel
a_ideal(::ZeroIdeal,V,T,z) = zero(Base.promote_eltype(V,T,z))

#just for completion
function eos_g(model::IdealModel,p,T,z)
    R = Rgas(model)
    RT = R*T
    n = sum(z)
    V = n*RT/p
    return n*RT*(a_ideal(model,V,T,z) + 1)
end
