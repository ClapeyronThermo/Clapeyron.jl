function correct_composition_derivative(model,V,T)
    z = SA[1.0]
    μ = only(Clapeyron.VT_chemical_potential(model,V,T,z))
    g = VT_gibbs_free_energy(model,V,T,z)
    @show μ
    @show g
    return μ,g
end
