function a_res(model::RKFamily,z,v,T)
    x = z/sum(z[i] for i in model.components)
    n = sum(z)

    a = sum(sum(model.params.a[union(i,j)]*x[i]*x[j] for j in model.components) for i in model.components)
    b = sum(sum(model.params.b[union(i,j)]*x[i]*x[j] for j in model.components) for i in model.components)
    Tc = model.params.Tc
    return -log(1-n*b/v)-a/(R̄*T*b*√(T/Tc))*log(1+n*b/v)
end
