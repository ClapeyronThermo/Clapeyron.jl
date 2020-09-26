using NamedArrays

export create_z

function create_z(model::EoS, z::AbstractArray)
    return NamedArray(z, model.components)
end
