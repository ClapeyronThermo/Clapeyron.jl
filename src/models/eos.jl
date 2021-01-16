function eos(model::EoSModel, V, T, z)
    return N_A*k_B*sum(z)*T * (a_ideal(model.idealmodel,V,T,z)+a_res(model,V,T,z))
end

function eos(model::IdealModel, V, T, z)
    return N_A*k_B*sum(z)*T * a_ideal(model,V,T,z)
end

# function eos(model::CubicModel, V, T, z)
#     return N_A*k_B*sum(z)*T * a_tot(model,V,T,z)
# end

"""
The EoS is extensible to other types of equations.

julia > function Eos(model::SRK, p, T, z)
            return some_function(p, T, z)
        end


Then create some get_properties functions for other equations.
"""



"""
    quadratic_mixing_rule(op, x)

returns a quadratic mixing rule, given a function `op(i,j)` and a vector of molar fractions `x`, it does a simplification of:

    ∑(op(i,j)*x[i]*x[j]) ,i = 1:N, j = 1:N

considering `op(i,j) == op(j,i)` this can be simplified to:

    ∑(op(i,j)*x[i]*x[j]) = ∑op(i,i)x[i]*x[i] + 2∑(op(i,j)*x[i]*x[j]),i = 1:N, j = 1:i - 1
"""

function quadratic_mixing_rule(op, x)
    N = length(x)
    @boundscheck checkbounds(x, N)
    @inbounds begin
        res1 = zero(eltype(x))
        for i = 1:N
            res1 += op(i,i) * x[i]^2
            for j = 1:i - 1
                res1 += 2 * x[i] * x[j] * op(i, j)
            end
        end
    end
    return res1
end
