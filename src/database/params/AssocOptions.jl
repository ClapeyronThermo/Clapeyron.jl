"""
    AssocOptions(;rtol = 1e-12,atol = 1e-12,max_iters = 1000,dampingfactor = 0.5,combining =:nocombining)

Struct containing iteration parameters for the solver of association sites.

the combining option controls the type of combining rule applied to the association strength:
- `sparse_nocombining` (default). Does not perform any combination rules over the association strength, and returns a sparse matrix.
- `dense_nocombining`. Does not perform any combination rules over the association strength, and returns a dense matrix.
- `elliott` combining rule: `Δ(i,j,a,b) = √(Δ(i,i,a,b)*Δ(j,j,a,b))`. Returns a dense matrix.

!!! info "Association Scheme matters"

    Elliott's rule requires that both `Δ(i,i,a,b)` and  `Δ(j,j,a,b)` are non-zero, that means that components that don't self associate will not be combined.
"""
@Base.kwdef struct AssocOptions <: ClapeyronParam
    rtol::Float64 = 1e-12
    atol::Float64 = 1e-12
    max_iters::Int = 1000
    dampingfactor::Float64 = 0.5
    combining::Symbol = :sparse_nocombining
end

is_splittable(::AssocOptions) = false

