
function eos(model::IdealModel, V, T, z=SA[1.0])
    return N_A*k_B*sum(z)*T * a_ideal(model,V,T,z)
end

for f in (:eos_res,:a_res,:VT_entropy_res,:VT_gibbs_free_energy_res,:VT_helmholtz_free_energy_res)
    @eval begin
        function $f(model::IdealModel, V, T, z=SA[1.])
            return zero(V+T+first(z))
        end
    end
end


function volume_impl(model::IdealModel,p,T,z=SA[1.0],phase=:unknown,threaded=false,vol0 = nothing)
    return sum(z)*R̄*T/p
end

function VT_chemical_potential_res(model::IdealModel, V, T, z)
    return false .* z
end

function VT_entropy(model::IdealModel, V, T, z=SA[1.])
    return -∂f∂T(model,V,T,z)
end

function VT_internal_energy(model::IdealModel, V, T, z=SA[1.])
    dA, A = ∂f(model,V,T,z)
    ∂A∂V, ∂A∂T = dA
    return A - T*∂A∂T
end

function VT_enthalpy(model::IdealModel, V, T, z=SA[1.])
    dA, A = ∂f(model,V,T,z)
    ∂A∂V, ∂A∂T = dA
    return A - V*∂A∂V - T*∂A∂T
end

function VT_gibbs_free_energy(model::IdealModel, V, T, z=SA[1.])
    A = eos(model,V,T,z)
    ∂A∂V = -sum(z)*R̄*T/V
    return A - V*∂A∂V
end

lb_volume(model::IdealModel,z=SA[1.0]) = zero(eltype(z))

idealmodel(model::IdealModel) = nothing

