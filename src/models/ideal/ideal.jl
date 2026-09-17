
a_res(model::IdealModel, V, T, z) = zero(Base.promote_eltype(model, V, T, z))

function volume_impl(model::IdealModel, p, T, z, phase, threaded, vol0)
    return sum(z)*R̄*T/p
end

lb_volume(model::IdealModel, T, z) = zero(eltype(z))

idealmodel(model::IdealModel) = model

@newmodelsingleton ZeroIdeal IdealModel
a_ideal(::ZeroIdeal, V, T, z) = zero(Base.promote_eltype(V, T, z))

#just for completion
function eos_g(model::IdealModel, p, T, z)
    R = Rgas(model)
    RT = R*T
    n = sum(z)
    V = n*RT/p
    return n*RT*(a_ideal(model, V, T, z) + 1)
end

"""
update_mw!(model,mw)

updates the `Mw` field from the ideal model.
Used on GC models that calculate the molecular weight as sums of contributions instead of looking for it in a database.
"""
function update_mw!(model::M, mw) where M
    if hasfield(M, :params)
        params = model.params
        if hasfield(typeof(params), :Mw)
            Mw = params.Mw.values
            Mw .= mw
        end
    end
end
