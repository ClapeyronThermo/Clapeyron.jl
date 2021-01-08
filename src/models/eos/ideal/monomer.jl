function a_ideal(model::Monomer, z, v, T)
    x = z/sum(z)
    return sum(x[i]*log(z[i]*N_A/v*Λ(model, z, v, T, i)^3) for i in model.components)-1
end

function Λ(model::Monomer, z, v, T, i)
    Mw = model.params.Mw[i]
    return h/√(k_B*T*Mw/N_A)
end
