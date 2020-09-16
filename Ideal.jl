function a_ideal(model,z,v,T)
    x = z/sum(z[i] for i in model.components)
    return sum(x[i]*log(z[i]/v) for i in model.components)-1
end
