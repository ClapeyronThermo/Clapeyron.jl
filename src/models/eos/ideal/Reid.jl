function a_ideal(model::Reid, z, v, T)
    x = z/sum(z)
    poly_coef = model.params.poly_coef
    return sum(x[i]*(log(z[i]/v) + 1/(R̄*T)*(sum(poly_coef[i][k]/k*(T^k-298^k) for k in 1:4)) -
                                   1/R̄*((poly_coef[i][1]-R̄)*log(T/298)+sum(poly_coef[i][k]/(k-1)*(T^(k-1)-298^(k-1)) for k in 2:4))) for i in model.components)-1
end
