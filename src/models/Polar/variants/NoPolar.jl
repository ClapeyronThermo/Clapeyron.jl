@newmodelsingleton NoPolar PolarModel

function polar_model(model::EoSModel)
    return NoPolar()
end

function polar_data(model,V,T,z,polar_model,_data = @f(data))
    return nothing
end

function a_polar(model,V,T,z,polar_model::NoPolar,_data,_polar_data)
    return zero(@f(Base.promote_eltype))
end
