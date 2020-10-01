function a_ideal(model::SAFT, z, v, T)
    x = z/sum(z)
    return sum(x[i]*log(z[i]/v) for i in model.components)-1
end
