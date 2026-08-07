#=
struct FixedComps{𝕋,𝕍,𝕀} <: EoSModel
    model::𝕋
    indices::𝕀
    values::𝕍
    n::Int
end

Base.length(model::FixedComps) = length(model.indices)


function expand_compositions(fracmodel::FixedComps,z)
    fixed_components = model.n
    @assert fixed_components == length(z)
    z2 = index_expansion(z,model.indices)
    k = 0
    @inbounds for i in eachindex(z2)
        if iszero(z2[i])
            k += 1
            z2[i] = fracmodel.values[k] 
        end
    end
    return z2
end

for f in (:eos, :eos_res, :a_res)
    @eval begin
    ($f)(model::FixedComps,V,T,z) = ($f)(model.model,V,T,expand_compositions(model,z))
    end
end =#