struct Thol2016LJ{T} <: LJPotentialModel
    model::T
end

Thol2016LJ() = Thol2016LJ(TholLJ_constructor())

a_potential(model::Thol2016LJ,V,T,_) = a_res(model.model,V,T,SA[1.0])

