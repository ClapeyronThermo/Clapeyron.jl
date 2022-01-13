function init_model(idealmodel::Type{<:IdealModel},components,userlocations,verbose)
    verbose && @info("""Now creating ideal model:
    $idealmodel""")
    return idealmodel(components;userlocations,verbose)
end

function eos(model::IdealModel, V, T, z=SA[1.0])
    negative_vt(V,T) && return nan_num(V,T,z)
    return N_A*k_B*sum(z)*T * a_ideal(model,V,T,z)
end

function volume(model::IdealModel,p,T,z=SA[1.0];phase=:unknown,threaded=false)
    return sum(z)*R̄*T/p
end