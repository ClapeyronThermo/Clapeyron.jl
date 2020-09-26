using NamedArrays

export create_z

function create_z(model::EoS, z::Array{Float64,1})
    return NamedArray(z, model.components)
end
