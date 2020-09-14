function a_ideal(model::PcSaftFamily,z,v,T)
    x = z/sum(z[i] for i in model.components)
    return sum(x[i]*log(x[i]/v) for i in model.components)-1
end
