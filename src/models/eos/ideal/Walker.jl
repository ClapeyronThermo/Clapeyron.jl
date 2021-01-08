function a_ideal(model::Walker, z, v, T)
    x = z/sum(z)
    Nrot = model.params.Nrot
    θ_vib = model.params.theta_V
    g_vib = model.params.deg_V
    return sum(x[i]*(log(z[i]*N_A/v*Λ(model, z, v, T, i)^3)-Nrot[i]/2*log(T)+sum(g_vib[i][v]*(θ_vib[i][v]/2/T+log(1-exp(-θ_vib[i][v]/T))) for v in 1:4)) for i in model.components)-1
end

function Λ(model::Walker, z, v, T, i)
    Mw = model.params.Mw[i]
    return h/√(k_B*T*Mw/N_A)
end
