
function eos_impl(model::IdealModel, V, T, z)
    return Rgas(model)*sum(z)*T *(a_ideal(model,V,T,z) + reference_state_eval(model,V,T,z))
end

for f in (:eos_res,:a_res,:VT_entropy_res,:VT_gibbs_free_energy_res,:VT_helmholtz_free_energy_res)
    @eval begin
        function $f(model::IdealModel, V, T, z=SA[1.])
            return zero(V+T+first(z))
        end
    end
end

function a_ideal(model::IdealModel, V, T, z)
    #a_ideal = ∑xᵢlogρᵢ + ∑xᵢa₀ᵢ(T)
    ∑xᵢlogρᵢ = sum(Base.Fix2(xlogx,1/V),z)/sum(z)
    ∑xᵢa₀ᵢ = a_ideal_T(model,T,z)
    return ∑xᵢa₀ᵢ + ∑xᵢlogρᵢ
end

function volume_impl(model::IdealModel,p,T,z,phase,threaded,vol0)
    return sum(z)*R̄*T/p
end

function VT_chemical_potential_res(model::IdealModel, V, T, z)
    return false .* z
end

function VT_entropy(model::IdealModel, V, T, z=SA[1.])
    return -∂f∂T(model,V,T,z)
end

function VT_internal_energy(model::IdealModel, V, T, z=SA[1.])
    V₀ = oneunit(V)
    A, ∂A∂T = f∂fdT(model,V₀,T,z)
    return A - T*∂A∂T
end

function VT_enthalpy(model::IdealModel, V, T, z=SA[1.])
    #idealmodels don't have a volume dependence for enthalpy
    if hasmethod(a_ideal_T,Tuple{typeof(model),Any,Any})
    
    end
    V₀ = oneunit(V)
    A, ∂A∂T = f∂fdT(model,V₀,T,z)
    return A + sum(z)*Rgas(model)*T - T*∂A∂T
end

function VT_gibbs_free_energy(model::IdealModel, V, T, z=SA[1.])
    return eos(model,V,T,z) + sum(z)*Rgas(model)*T
end

function VT_helmholtz_free_energy(model::IdealModel, V, T, z=SA[1.])
    return eos(model,V,T,z)
end

lb_volume(model::IdealModel,z) = zero(eltype(z))

idealmodel(model::IdealModel) = model
