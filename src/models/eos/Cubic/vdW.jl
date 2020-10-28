function a_tot(model::vdWFamily,z,v,T)
    x = z/sum(z[i] for i in model.components)
    N = sum(z)
    a = sum(sum(model.params.a_vdW[union(i,j)]*x[i]*x[j] for j in model.components) for i in model.components)
    b = sum(sum(model.params.b_vdW[union(i,j)]*x[i]*x[j] for j in model.components) for i in model.components)
    return -log(V-N*b)-a*N/(k_B*T*V)
end
