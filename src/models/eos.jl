function eos(model::EoSModel, V, T, z=@SVector [1.0])
    return N_A*k_B*sum(z)*T * (a_ideal(idealmodel(model),V,T,z)+a_res(model,V,T,z))
end

idealmodel(model::EoSModel) = model.idealmodel

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

#=fast shortcut to evaluate cubics, pressure is known
function ∂f∂v(model::ABCubicModel,v,t,z)
    @info "fast shortcut"
    a,b,p = cubic_abp(model,v,t,z)
    return -p
end
=#


function v_rackett(model,T)
    tc = only(paramvals(model.params.Tc))
    pc = only(paramvals(model.params.Pc))
    vc = only(paramvals(model.params.Vc))
    TT = promote_type(typeof(tc),typeof(T))
    zc = pc*vc/(R̄*tc)
    _2_7 = TT(2//7)
    _1 = TT(1.0)
    R̄*tc/pc*zc^(_1 + (_1 - T/tc)^(_2_7))
end