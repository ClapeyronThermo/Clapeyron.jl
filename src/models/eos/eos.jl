function eos(model::SAFT, z, v, T)
    return N_A*k_B*sum(z) * T*(a_ideal(model,z,v,T)+a_res(model,z,v,T))
end
    
"""
The EoS is extensible to other types of equations.

julia > function Eos(model::SRK, p, T, z)
            return some_function(p, T, z)
        end


Then create some get_properties functions for other equations.
"""
