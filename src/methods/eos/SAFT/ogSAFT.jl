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
    ϵ = model.parameters.epsilon[union(component,component)]
    σ = model.parameters.sigma[union(component,component)]
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
    σx = (sum(x[i]*x[j]*m[i]*m[j]*σ[union(i,j)]^3 for i in model.components for j in model.components)/mx^2)^(1/3)
    ϵx = (sum(x[i]*x[j]*m[i]*m[j]*σ[union(i,j)]^3*ϵ[union(i,j)] for i in model.components for j in model.components)/mx^2)/σx^3

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
    ϵx = sum(x[i]*x[j]*m[i]*m[j]*σ[union(i,j)]^3*ϵ[union(i,j)] for i in model.components for j in model.components)/sum(x[i]*x[j]*m[i]*m[j]*σ[union(i,j)]^3 for i in model.components for j in model.components)
    ηx = η(model,z,v,T)
    ρR = (6/sqrt(2)/π)*ηx
    TR = T/ϵx
    a_seg1 = ρR*(-8.5959-4.5424ρR-2.1268ρR^2+10.285ρR^3)
    a_seg2 = ρR*(-1.9075+9.9724ρR-22.216ρR^2+15.904ρR^3)
    return 1/TR*(a_seg1+a_seg2/TR)
end

## This is an attempt to make Twu et al.'s segment term; does not work yet
# function a_seg(model::ogSAFTFamily,z,v,T)
#     Bo = [1.31024,-3.80636,-2.37238,-0.798872,0.198761,1.47014,-0.786367,2.19465,5.75429,6.7822,-9.94904,-15.6162,86.643,18.527,9.04755,8.68282]
#     Ba = [3.79621,-6.14518,-1.84061,-2.77584,-0.420751,-5.66128,19.2144,-3.33443,33.0305,-5.90766,9.55619,-197.883,-61.2535,77.1802,-6.57983,0.0]
#     ω = 0.199
#     A = []
#     for i in 1:16
#         append!(A,Bo[i]+ω*Ba[i])
#     end
#     m = model.parameters.segment
#     σ = model.parameters.sigma
#     ϵ = model.parameters.epsilon
#     x = z/sum(z[i] for i in model.components)
#     mx = sum(x[i]*m[i] for i in model.components)
#     σx = (sum(x[i]*x[j]*m[i]*m[j]*σ[union(i,j)]^3 for i in model.components for j in model.components)/mx^2)^(1/3)
#     ϵx = (sum(x[i]*x[j]*m[i]*m[j]*σ[union(i,j)]^3*ϵ[union(i,j)] for i in model.components for j in model.components)/mx^2)/σx^3
#
#     ρ  = N_A*sum(z[i] for i in model.components)/v
#     ρR = ρ*σx^3
#     TR = T/ϵx
#
#     u_res = (A[2]/TR+2A[3]/TR^2+3A[4]/TR^3+5A[5]/TR^5)*ρR+1/2*A[7]/TR*ρR^2+
#             1/(2*A[16])*(3A[9]/TR^3+4A[10]/TR^4+5A[11]/TR^5)*(1-exp(-A[16]*ρR^2))+
#             1/(2*A[16]^2)*(3A[12]/TR^3+4A[13]/TR^4+5A[14]/TR^5)*(1-(1+A[16]*ρR^2)*exp(-A[16]*ρR^2))+
#             1/5*A[15]/TR*ρR^5
#     s_res = -log(ρ*R*T)-(A[1]-A[3]/TR^2-2A[4]/TR^3-4A[5]/TR^5)*ρR-1/2*A[6]*ρR^2-1/3*A[8]*ρR^3+
#             1/(2*A[16])*(2A[9]/TR^3+3A[10]/TR^4+4A[11]/TR^5)*(1-exp(-A[16]*ρR^2))+
#             1/(2*A[16]^2)*(2A[12]/TR^3+3A[13]/TR^4+4A[14]/TR^5)*(1-(1+A[16]*ρR^2)*exp(-A[16]*ρR^2))
#     a_res = u_res-s_res
#     println(A)
#     println(u_res)
#     println(s_res)
#     return mx*(a_res)
# end
