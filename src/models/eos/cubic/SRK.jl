function a_res(model::SRKFamily,z,v,T)
    x = z/sum(z[i] for i in model.components)
    n = sum(z)

    āᾱ = sum(sum(model.params.a[union(i,j)]*√(α(model,T,i)*α(model,T,j))*x[i]*x[j] for j in model.components) for i in model.components)
    b̄  = sum(sum(model.params.b[union(i,j)]*x[i]*x[j] for j in model.components) for i in model.components)
    return -log(1-n*b̄/v)-āᾱ/(R̄*T*b̄)*log(1+n*b̄/v)
end

function α(model::SRKFamily,T,i)
    Tc = model.params.Tc[i]
    ω  = model.params.acentric_fac[i]
    return (1+(0.480+1.547*ω-0.176*ω^2)*(1-√(T/Tc)))^2
end
