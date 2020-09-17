N_A = 6.02214086e23
k_B = 1.38064852e-23
R   = N_A*k_B
function a_res(model::ogSAFTFamily,z,v,T)
    return a_seg(model,z,v,T) + a_chain(model,z,v,T)
end

function a_seg(model::ogSAFTFamily,z,v,T)
    m = model.parameters.segment
    x = z/sum(z[i] for i in model.components)
    m̄ = sum(x[i]*m[i] for i in model.components)

    return m̄*(a_hs(model,z,v,T)+a_disp(model,z,v,T))
end

function a_chain(model::ogSAFTFamily,z,v,T)
    x = z/sum(z[i] for i in model.components)
    m = model.parameters.segment
    return sum(x[i]*(1-m[i])*log(g_hsij(model,z,v,T, i, i)) for i in model.components)
end

function d(model::ogSAFTFamily,z,v,T, component)
    ϵ = model.parameters.epsilon[component,component]
    σ = model.parameters.sigma[component,component]
    m = model.parameters.segment[component]
    fm = 0.0010477+0.025337*(m-1)/m
    f = (1+0.2977T/ϵ)/(1+0.33163T/ϵ+fm*(T/ϵ)^2)
    return σ * f
end

function dx(model::ogSAFTFamily,z,v,T)
    x = z/sum(z[i] for i in model.components)
    m = model.parameters.segment
    σ = model.parameters.sigma
    ϵ = model.parameters.epsilon

    mx = sum(x[i]*m[i] for i in model.components)
    σx = (sum(x[i]*x[j]*m[i]*m[j]*σ[i,j]^3 for i in model.components for j in model.components)/mx^2)^(1/3)
    ϵx = (sum(x[i]*x[j]*m[i]*m[j]*σ[i,j]^3*ϵ[i,j] for i in model.components for j in model.components)/mx^2)/σx^3

    fm = 0.0010477+0.025337*(mx-1)/mx
    f = (1+0.2977T/ϵx)/(1+0.33163T/ϵx+fm*(T/ϵx)^2)
    return σx * f
end

function ζn(model::ogSAFTFamily,z,v,T, n)
    x = z/sum(z[i] for i in model.components)
    m = model.parameters.segment
    return N_A*sum(z[i] for i in model.components)*π/6/v * sum(x[i]*m[i]*d(model,z,v,T, i)^n for i in model.components)
end

function η(model::ogSAFTFamily,z,v,T)
    m = model.parameters.segment
    x = z/sum(z[i] for i in model.components)
    m̄ = sum(x[i]*m[i] for i in model.components)
    return N_A*sum(z[i] for i in model.components)*π/6/v*dx(model,z,v,T)^3*m̄
end

function g_hsij(model::ogSAFTFamily,z,v,T, i, j)
    di = d(model,z,v,T, i)
    dj = d(model,z,v,T, j)
    ζ2 = ζn(model,z,v,T, 2)
    ζ3 = ζn(model,z,v,T, 3)
    return 1/(1-ζ3) + di*dj/(di+dj)*3ζ2/(1-ζ3)^2 + (di*dj/(di+dj))^2*2ζ2^2/(1-ζ3)^3
end

function a_hs(model::ogSAFTFamily,z,v,T)
    ηx = η(model,z,v,T)
    return (4ηx-3ηx^2)/(1-ηx)^2
end

function a_disp(model::ogSAFTFamily,z,v,T)
    m = model.parameters.segment
    σ = model.parameters.sigma
    ϵ = model.parameters.epsilon
    x = z/sum(z[i] for i in model.components)
    ϵx = sum(x[i]*x[j]*m[i]*m[j]*σ[i,j]^3*ϵ[i,j] for i in model.components for j in model.components)/sum(x[i]*x[j]*m[i]*m[j]*σ[i,j]^3 for i in model.components for j in model.components)
    ηx = η(model,z,v,T)
    ρR = (6/sqrt(2)/π)*ηx
    TR = T/ϵx
    a10 = ρR*(-8.5959-4.5424ρR-2.1268ρR^2+10.285ρR^3)
    a20 = ρR*(-1.9075+9.9724ρR-22.216ρR^2+15.904ρR^3)
    return 1/TR*(a10+a20/TR)
end
