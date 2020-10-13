function a_res(model::SAFTgammaMieFamily, z, V, T)
    return a_mono(model,z,V,T) + a_chain(model,z,V,T) + a_assoc(model,z,V,T)
end

function a_mono(model::SAFTgammaMieFamily, z, V, T)
    return a_HS(model,z,V,T) + a1(model,z,V,T) + a2(model,z,V,T) + a3(model,z,V,T)
end

function a_chain(model::SAFTgammaMieFamily, z, V, T)
    return 0
end

function a_assoc(model::SAFTgammaMieFamily, z, V, T)
    return 0
end

function a_hs(model::SAFTgammaMieFamily, z, V, T)
    v = model.groups
    vst = model.params.segment
    S = model.params.shapefactor
    return sum(z[i] * sum(v[i][k]*vst[k]*S[k] for k in @groups) for i in @components) * @f(as_hs)
end

function as_hs(model::SAFTgammaMieFamily, z, V, T)
    ζ0   = ζn(model, z,Vol,Temp, 0)
    ζ1   = ζn(model, z,Vol,Temp, 1)
    ζ2   = ζn(model, z,Vol,Temp, 2)
    ζ3   = ζn(model, z,Vol,Temp, 3)
    return 6*Vol/π/@f(ρ_s)*(3ζ1*ζ2/(1-ζ3) + ζ2^3/(ζ3*(1-ζ3)^2) + (ζ2^3/ζ3^2-ζ0)*log(1-ζ3))
end
    


function d(model::SAFTgammaMieFamily, z, V, T)
    ϵ           = model.params.epsilon[union(i,i)]
    σ           = model.params.sigma[union(i,i)]
    λR          = model.params.lambdaR[union(i,i)]
    λA          = model.params.lambdaA[union(i,i)]
    u           = [0.26356031971814109102031,1.41340305910651679221800,3.59642577104072208122300,7.08581000585883755692200,12.6408008442757826594300]
    w           = [0.5217556105828086524759,0.3986668110831759274500,7.5942449681707595390e-2,3.6117586799220484545e-3,2.3369972385776227891e-5]
    θ           = CMie(model,z,V,T,i,i)*ϵ/T
    return σ*(1-sum(w[j]*(θ./(θ+u[j]))^(1/λR)*(exp(θ*(1/(θ./(θ+u[j]))^(λA/λR)-1))/(u[j]+θ)/λR) for j in 1:5))
end

function CMie(model::SAFTVRMieFamily, z, V, T, i, j)
    λR          = model.params.lambdaR
    λA          = model.params.lambdaA
    return (λR[union(i,j)]/(λR[union(i,j)]-λA[union(i,j)]))*(λR[union(i,j)]/λA[union(i,j)])^(λA[union(i,j)]/(λR[union(i,j)]-λA[union(i,j)]))
end
