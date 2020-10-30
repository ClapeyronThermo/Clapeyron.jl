function a_tot(model::PRFamily,z,v,T)
    x = z/sum(z[i] for i in model.components)
    n = sum(z)
    āᾱ = sum(sum(model.params.a[union(i,j)]*√(α(model,T,i)*α(model,T,j))*x[i]*x[j] for j in model.components) for i in model.components)
    b̄  = sum(sum(model.params.b[union(i,j)]*x[i]*x[j] for j in model.components) for i in model.components)
    return -log(v-n*b̄)+āᾱ/(R̄*T*b̄*2^(3/2))*log((2*v-2^(3/2)*b̄*n+2*b̄*n)/(2*v+2^(3/2)*b̄*n+2*b̄*n))
end

function α(model::PRFamily,T,i)
    Tc = model.params.Tc[i]
    ω  = model.params.acentric_fac[i]
    return (1+(0.37464+1.54226*ω-0.26992*ω^2)*(1-√(T/Tc)))^2
end
