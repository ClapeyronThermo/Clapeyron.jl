
function eos(model::IdealModel, V, T, z=SA[1.0])
    negative_vt(V,T) && return nan_num(V,T,z)
    return N_A*k_B*sum(z)*T * a_ideal(model,V,T,z)
end

function eos_res(model::IdealModel, V, T, z=SA[1.0])
    return zero(V+T+first(z))
end

function volume_impl(model::IdealModel,p,T,z=SA[1.0],phase=:unknown,threaded=false,vol0 = nothing)
    return sum(z)*R̄*T/p
end

function VT_chemical_potential_res(model::IdealModel, vol, T, z)
    return false .* z
end

lb_volume(model::IdealModel,z=SA[1.0]) = zero(eltype(z))

idealmodel(model::IdealModel) = nothing

