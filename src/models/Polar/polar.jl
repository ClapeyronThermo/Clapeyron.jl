function a_polar(model::EoSModel,V,T,z,_data = @f(data))
    polar_model = polar_model(model)
    _polar_data = @f(polar_data,polar_model,_data)
    return a_polar(model,V,T,z,polar_model,_data,_polar_data)
end

function polar_comps(model, V, T, z)
    μ̄² = model.params.dipole2.values
    Q̄² = model.params.quadrupole2.values
    dipole_comps = findall(!iszero,μ̄²)
    quadrupole_comps = findall(!iszero,Q̄²)
    return dipole_comps, quadrupole_comps
end

include("NoPolar.jl")
#=include("GrossVrabec.jl")
include("GrossVrabecQuatrupolar.jl") ="
