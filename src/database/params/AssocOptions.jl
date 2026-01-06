"""
    AssocOptions(;rtol = 1e-12,atol = 1e-12,max_iters = 1000,dampingfactor = 0.5,combining =:nocombining)
Struct containing iteration parameters for the solver of association sites.
the combining option controls the type of combining rule applied to the association strength:
- `nocombining` (default). Does not perform any combination rules.
- `:cr1`: "combining rule - 1":
    ```
    ε[i,j][a,b] = (ε[i,i][a,b] + ε[j,j][a,b])/2
    β[i,j][a,b] = √(β[i,i][a,b] * β[j,j][a,b])
    ```
- `:esd`,`:elliott`: Elliott–Suresh–Donohue combining rule:
    ```
    ε[i,j][a,b] = (ε[i,i][a,b] + ε[j,j][a,b])/2
    β[i,j][a,b] = √(β[i,i][a,b] * β[j,j][a,b]) * (σ[i]*σ[j]/σ[i,j])^3
    ```
- `:esd_runtime`,`:elliott_runtime`: combining rule, performed at runtime:
    ```
    Δ[i,j][a,b] = √(Δ[i,i][a,b] * Δ[j,j][a,b]) 
    ```
- `:dufal`,`mie15`: combining rule used for `SAFTVRMie15`
    ```
    ε[i,j][a,b] = (ε[i,i][a,b] + ε[j,j][a,b])/2
    β[i,j][a,b] = (∛β[i,i][a,b] + ∛β[j,j][a,b])^3 / 8
    ```

!!! info "Association Scheme matters"
    all combining rules implicitly requires that both `Δ(i,i,a,b)` and  `Δ(j,j,a,b)` are non-zero, that means that components that don't self associate will not be combined.
"""
@Base.kwdef struct AssocOptions <: ClapeyronParam
    rtol::Float64 = 1e-12
    atol::Float64 = 1e-12
    max_iters::Int = 1000
    dampingfactor::Float64 = 0.5
    combining::Symbol = :nocombining
end

#allows overloading default assoc_options
default_assoc_options(m::EoSModel) = default_assoc_options(parameterless_type(m))
default_assoc_options(m) = AssocOptions()

Base.show(io::IO,::MIME"text/plain",options::AssocOptions) = show_as_namedtuple(io,options)
Base.show(io::IO,options::AssocOptions) = show_as_namedtuple(io,options)

is_splittable(::AssocOptions) = false

__init_assoc_options_kw(::Nothing) = AssocOptions()
__init_assoc_options_kw(s::Symbol) = AssocOptions(combining = s)
__init_assoc_options_kw(ref::AssocOptions) = ref