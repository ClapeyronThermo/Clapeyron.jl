
a_res(model::IdealModel,V,T,z) = zero(Base.promote_eltype(model,V,T,z))

function a_ideal(model::IdealModel, V, T, z)
    #a_ideal = ∑xᵢlogρᵢ(V,x) + ∑xᵢa₀ᵢ(T,x)
    ∑xᵢlogρᵢ = sum(Base.Fix2(xlogx,1/V),z)/sum(z)
    ∑xᵢa₀ᵢ = a_ideal_T(model,T,z)
    return ∑xᵢa₀ᵢ + ∑xᵢlogρᵢ
end

function volume_impl(model::IdealModel,p,T,z,phase,threaded,vol0)
    return sum(z)*R̄*T/p
end

function VT_gibbs_free_energy(model::IdealModel, V, T, z=SA[1.])
    return eos(model,V,T,z) + sum(z)*Rgas(model)*T
end

lb_volume(model::IdealModel,z) = zero(eltype(z))

idealmodel(model::IdealModel) = model
