function a_res(model::vdWFamily,z,v,T)
    x = z/sum(z[i] for i in model.components)
    n = sum(z)
    a = sum(sum(model.params.a[union(i,j)]*x[i]*x[j] for j in model.components) for i in model.components)
    b = sum(sum(model.params.b[union(i,j)]*x[i]*x[j] for j in model.components) for i in model.components)
    return -log(1-n*b/v)-a*n/(R̄*T*v)
end
