function eos(model::EoSModel, V, T, z)
    return N_A*k_B*sum(z)*T * (a_ideal(idealmodel(model),V,T,z)+a_res(model,V,T,z))
end

idealmodel(model::EoSModel) = model.idealmodel
idealmodel(model::IAPWS95) = IAPWS95Ideal()

function eos(model::IdealModel, V, T, z)
    return N_A*k_B*sum(z)*T * a_ideal(model,V,T,z)
end

function a_res(model::CubicModel, V, T, z)
    return @f(a_tot) + log(V)  # + f(x)
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

