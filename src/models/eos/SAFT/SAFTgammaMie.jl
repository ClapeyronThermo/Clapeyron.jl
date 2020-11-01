function Âres(model::SAFTgammaMieFamily, z, V, T)
    return Âmono(model,z,V,T) + Âchain(model,z,V,T) + Âassoc(model,z,V,T)
end

function Âmono(model::SAFTgammaMieFamily, z, V, T)
    return ÂHS(model,z,V,T) + Â_1(model,z,V,T) + Â_2(model,z,V,T) + Â_3(model,z,V,T)
end

function Âchain(model::SAFTgammaMieFamily, z, V, T)
    x = z/sum(z)
    v = model.group_multiplicities
    vst = model.params.segment
    S = model.params.shapefactor
    return -sum(x[i] * sum((v[i][k]*vst[k]*S[k]-1) for k in @groups(i)) * log(@f(g_Mie,i)) for i in @comps)
end

function Âassoc(model::SAFTgammaMieFamily, z, V, T)
    return 0
end

function ÂHS(model::SAFTgammaMieFamily, z, V, T)
    ζ_0   = @f(ζ,0)
    ζ_1   = @f(ζ,1)
    ζ_2   = @f(ζ,2)
    ζ_3   = @f(ζ,3)
    ρ = sum(z)*N_A/V
    return 6/π/ρ * (3ζ_1*ζ_2/(1-ζ_3) + ζ_2^3/(ζ_3*(1-ζ_3)^2) + (ζ_2^3/ζ_3^2-ζ_0)*log(1-ζ_3))
end

#= for n in 1:3 =#
#=     @eval   function $(Symbol(:a,Symbol(n)))(model::SAFTgammaMieFamily, z, V, T) =#
#=                 x = z/sum(z) =#
#=                 v = model.group_multiplicities =#
#=                 vst = model.params.segment =#
#=                 S = model.params.shapefactor =#

#=                 return 1/(kB*T)^n * sum(x[i] * sum(v[i][k]*vst[k]*S[k] for k in @groups(i)) for i in @comps) * @eval $(@f(Symbol(:â,Symbol(n)))) =#
#=             end =#

#=     @eval   function $(Symbol(:â,Symbol(n)))(model::SAFTgammaMieFamily, z, V, T) =#
#=                 return sum(@f(x_S,k)*@f(x_S,l) * @eval $(@f(Symbol(:â,Symbol(n)),k,l)) for k in @groups for l in @groups) =#
#=             end =#
#= end =#
#
#= function â(model::SAFTgammaMieFamily, z, V, T) =#

function Â_1(model::SAFTgammaMieFamily, z, V, T)
    x = z/sum(z)
    v = model.group_multiplicities
    vst = model.params.segment
    S = model.params.shapefactor

    return 1/T * sum(x[i]*sum(v[i][k]*vst[k]*S[k] for k in @groups(i)) for i in @comps) * @f(a_1)
end
function Â_2(model::SAFTgammaMieFamily, z, V, T)
    x = z/sum(z)
    v = model.group_multiplicities
    vst = model.params.segment
    S = model.params.shapefactor

    return 1/T^2 * sum(x[i]*sum(v[i][k]*vst[k]*S[k] for k in @groups(i)) for i in @comps) * @f(a_2)
end
function Â_3(model::SAFTgammaMieFamily, z, V, T)
    x = z/sum(z)
    v = model.group_multiplicities
    vst = model.params.segment
    S = model.params.shapefactor

    return 1/T^3 * sum(x[i]*sum(v[i][k]*vst[k]*S[k] for k in @groups(i)) for i in @comps) * @f(a_3)
end

function a_1(model::SAFTgammaMieFamily, z, V, T)
    return sum(@f(x_S,k)*@f(x_S,l)*@f(a_1,k,l) for k in @groups for l in @groups)
end
function a_2(model::SAFTgammaMieFamily, z, V, T)
    return sum(@f(x_S,k)*@f(x_S,l)*@f(a_2,k,l) for k in @groups for l in @groups)
end
function a_3(model::SAFTgammaMieFamily, z, V, T)
    return sum(@f(x_S,k)*@f(x_S,l)*@f(a_3,k,l) for k in @groups for l in @groups)
end

function a_1(model::SAFTgammaMieFamily, z, V, T, k, l)
    σ = model.params.sigma[union(k,l)]
    λa = model.params.lambda_a[union(k,l)]
    λr = model.params.lambda_r[union(k,l)]

    x_0 = σ/@f(d,k,l)
    return @f(C,λa,λr) * (x_0^λa*(@f(aS_1,k,l,λa)+@f(B,k,l,λa)) - x_0^λr*(@f(aS_1,k,l,λr)+@f(B,k,l,λr)))
end
function a_2(model::SAFTgammaMieFamily, z, V, T, k, l)
    σ = model.params.sigma[union(k,l)]
    ϵ = model.params.epsilon[union(k,l)]
    λa = model.params.lambda_a[union(k,l)]
    λr = model.params.lambda_r[union(k,l)]

    x_0 = σ/@f(d,k,l)
    return 1/2*@f(KHS)*(1+@f(𝜒,k,l))*ϵ*@f(C,λa,λr)^2 * (
           x_0^(2λa)*(@f(aS_1,k,l,2λa) + @f(B,k,l,2λa))
         - 2x_0^(λa+λr)*(@f(aS_1,k,l,λa+λr) + @f(B,k,l,λa+λr))
         + x_0^(2λr)*(@f(aS_1,k,l,2λr) + @f(B,k,l,2λr)) )
end
function a_3(model::SAFTgammaMieFamily, z, V, T, k, l)
    ϵ = model.params.epsilon[union(k,l)]

    ζst_X_ = @f(ζst_X)
    return -ϵ^3*@f(f,k,l,4)*ζst_X_ * exp(@f(f,k,l,5)*ζst_X_ + @f(f,k,l,6)*ζst_X_^2)
end

function B(model::SAFTgammaMieFamily, z, V, T, k, l, λ)
    ϵ = model.params.epsilon[union(k,l)]
    σ = model.params.sigma[union(k,l)]
    x_0 = σ/@f(d,k,l)
    ζ_X_ = @f(ζ_X)
    I = (1-x_0^(3-λ))/(λ-3)
    J = (1-x_0^(4-λ)*(λ-3)+x_0^(3-λ)*(λ-4))/((λ-3)*(λ-4))
    return 2π*@f(ρ_S)*@f(d,k,l)^3*ϵ * ((1-ζ_X_/2)/(1-ζ_X_)^3*I-9ζ_X_*(1+ζ_X_)/(2(1-ζ_X_)^3)*J)
end

function aS_1(model::SAFTgammaMieFamily, z, V, T, k, l, λ)
    ϵ = model.params.epsilon[union(k,l)]
    ζeff_ = @f(ζeff, λ)
    return -2π*@f(ρ_S) * ϵ*@f(d,k,l)^3/(λ-3) * (1-ζeff_/2)/(1-ζeff_)^3
end

function ζ_X(model::SAFTgammaMieFamily, z, V, T)
    return π*@f(ρ_S)/6 * sum(@f(x_S,k)*@f(x_S,l)*@f(d,k,l)^3 for k in @groups for l in @groups)
end

function ζst_X(model::SAFTgammaMieFamily, z, V, T)
    σ = model.params.sigma
    return π*@f(ρ_S)/6 * sum(@f(x_S,k)*@f(x_S,l)*σ[union(k,l)]^3 for k in @groups for l in @groups)
end

function ζeff(model::SAFTgammaMieFamily, z, V, T, λ)
    A = [ 0.81096    1.7888   -37.578   92.284;
          1.02050  -19.341    151.26  -463.50 ;
         -1.90570   22.845   -228.14   973.92 ;
          1.08850   -6.1962   106.98  -677.64  ]
    ζ_X_ = @f(ζ_X)
    return A * [1; 1/λ; 1/λ^2; 1/λ^3] ⋅ [ζ_X_; ζ_X_^2; ζ_X_^3; ζ_X_^4]
end

function KHS(model::SAFTgammaMieFamily, z, V, T)
    ζ_X_ = @f(ζ_X)
    return (1-ζ_X_)^4/(1+4ζ_X_+4ζ_X_^2-4ζ_X_^3+ζ_X_^4)
end

function 𝜒(model::SAFTgammaMieFamily, z, V, T, k, l)
    ζst_X_ = @f(ζst_X)
    return @f(f,k,l,1)*ζst_X_ + @f(f,k,l,2)*ζst_X_^5 + @f(f,k,l,3)*ζst_X_^8
end

function f(model::SAFTgammaMieFamily, z, V, T, k, l, m)
    ϕ = [[7.5365557, -359.440,  1550.9, -1.199320, -1911.2800,  9236.9  ],
         [-37.604630,  1825.60, -5070.1,  9.063632,  21390.175, -129430 ],
         [71.745953, -3168.00,  6534.6, -17.94820, -51320.700,  357230  ],
         [-46.835520,  1884.20, -3288.7,  11.34027,  37064.540, -315530 ],
         [-2.4679820,- 0.82376, -2.7171,  20.52142,  1103.7420,  1390.2 ],
         [-0.5027200, -3.19350,  2.0883, -56.63770, -3264.6100, -4518.2 ],
         [8.0956883,  3.70900,  0.0000,  40.53683,  2556.1810,  4241.6  ]];
    λa = model.params.lambda_a[union(k,l)]
    λr = model.params.lambda_r[union(k,l)]

    α = @f(C,λa,λr)*(1/(λa-3)-1/(λr-3))
    return sum(ϕ[i+1][m]*α^i for i in 0:3)/(1+sum(ϕ[i+1][m]*α^(i-3) for i in 4:6))
end

function ζ(model::SAFTgammaMieFamily, z, V, T, m)
    return π/6*@f(ρ_S)*sum(@f(x_S,k)*@f(d,k)^m for k in @groups)
end

function ρ_S(model::SAFTgammaMieFamily, z, V, T)
    x = z/sum(z)
    v = model.group_multiplicities
    vst = model.params.segment
    S = model.params.shapefactor

    ρ = sum(z)*N_A/V
    return ρ * sum(x[i] * sum(v[i][k]*vst[k]*S[k] for k in @groups(i)) for i in @comps)
end

function x_S(model::SAFTgammaMieFamily, z, V, T, k)
    x = z/sum(z)
    v = model.group_multiplicities
    vst = model.params.segment
    S = model.params.shapefactor

    return sum(x[i]*v[i][k]*vst[k]*S[k] for i in @comps) / sum(x[i] * sum(v[i][l]*vst[l]*S[l] for l in @groups(i)) for i in @comps)
end

function d(model::SAFTgammaMieFamily, z, V, T, k)
    ϵ = model.params.epsilon[k]
    σ = model.params.sigma[k]
    λa = model.params.lambda_a[k]
    λr = model.params.lambda_r[k]
    u = [0.26356031971814109102031,1.41340305910651679221800,3.59642577104072208122300,7.08581000585883755692200,12.6408008442757826594300]
    w = [0.5217556105828086524759,0.3986668110831759274500,7.5942449681707595390e-2,3.6117586799220484545e-3,2.3369972385776227891e-5]
    θ = @f(C,λa,λr)*ϵ/T
    return σ*(1-sum(w[j]*(θ./(θ+u[j]))^(1/λr)*(exp(θ*(1/(θ./(θ+u[j]))^(λa/λr)-1))/(u[j]+θ)/λr) for j in 1:5))
end

function d(model::SAFTgammaMieFamily, z, V, T, k, l)
    # Lorentz mixing rule
    if k == l
        return @f(d, k)
    else
        return (@f(d, k) + @f(d, l))/2
    end
end

function C(model::SAFTgammaMieFamily, z, V, T, λa, λr)
    return (λr/(λr-λa)) * (λr/λa)^(λa/(λr-λa))
end

function ẑ(model::SAFTgammaMieFamily, z, V, T, i, k)
    v = model.group_multiplicities
    vst = model.params.segment
    S = model.params.shapefactor
    return v[i][k]*vst[k]*S[k] / sum(v[i][l]*vst[l]*S[l] for l in @groups(i))
end

function σ̄(model::SAFTgammaMieFamily, z, V, T, i)
    σ = model.params.sigma
    return cbrt(sum(sum(@f(ẑ,i,k)*@f(ẑ,i,l)*σ[union(k,l)]^3 for l in @groups) for k in @groups))
end

function d̄(model::SAFTgammaMieFamily, z, V, T, i)
    return cbrt(sum(sum(@f(ẑ,i,k)*@f(ẑ,i,l)*@f(d,k,l)^3 for l in @groups) for k in @groups))
end

function ϵ̄(model::SAFTgammaMieFamily, z, V, T, i)
    ϵ = model.params.epsilon
    return sum(sum(@f(ẑ,i,k)*@f(ẑ,i,l)*ϵ[union(k,l)] for l in @groups) for k in @groups)
end

function λ̄a(model::SAFTgammaMieFamily, z, V, T, i)
    λa = model.params.lambda_a
    return sum(sum(@f(ẑ,i,k)*@f(ẑ,i,l)*λa[union(k,l)] for l in @groups) for k in @groups)
end

function λ̄r(model::SAFTgammaMieFamily, z, V, T, i)
    λr = model.params.lambda_r
    return sum(sum(@f(ẑ,i,k)*@f(ẑ,i,l)*λr[union(k,l)] for l in @groups) for k in @groups)
end

function g_Mie(model::SAFTgammaMieFamily, z, V, T, i)
    ϵ̄_=@f(ϵ̄,i)
    return @f(g_HS,i)*exp(ϵ̄_/T*@f(g_1,i)/@f(g_HS,i)+(ϵ̄_/T)^2*@f(g_2,i)/@f(g_HS,i))
end

function g_HS(model::SAFTgammaMieFamily, z, V, T, i)
    x̄_0 = @f(σ̄,i)/@f(d̄,i)
    ζ_X_ = @f(ζ_X)
    k_0 = -log(1-ζ_X_) + (42ζ_X_-39ζ_X_^2+9ζ_X_^3-2ζ_X_^4)/(6*(1-ζ_X_)^3)
    k_1 = (ζ_X_^4+6ζ_X_^2-12ζ_X_)/(2*(1-ζ_X_)^3)
    k_2 = -3ζ_X_^2/(8(1-ζ_X_)^2)
    k_3 = (-ζ_X_^4+3ζ_X_^2+3ζ_X_)/(6(1-ζ_X_)^3)
    return exp(k_0+k_1*x̄_0+k_2*x̄_0^2+k_3*x̄_0^3)
end

function g_1(model::SAFTgammaMieFamily, z, V, T, i)
    λ̄a_ = @f(λ̄a,i)
    λ̄r_ = @f(λ̄r,i)
    σ̄_ = @f(σ̄,i)
    ϵ̄_ = @f(ϵ̄,i)
    d̄_ = @f(d̄,i)
    ρ_S_ = @f(ρ_S)
    x̄_0 = σ̄_/@f(d̄,i)
    C̄_ = @f(C,λ̄a_,λ̄r_)
    return 1/(2π*ϵ̄_*d̄_^3)*(3*@f(∂ā_1∂ρ_S,i) - C̄_*λ̄a_*x̄_0^λ̄a_*(@f(āS_1,i,λ̄a_)+@f(B̄,i,λ̄a_)/ρ_S_) + C̄_*λ̄r_*x̄_0^λ̄r_*(@f(āS_1,i,λ̄r_)+@f(B̄,i,λ̄r_)/ρ_S_))
end

function g_2(model::SAFTgammaMieFamily, z, V, T, i)
    return (1+@f(γ_c,i))*@f(gMCA_2,i)
end

function B̄(model::SAFTgammaMieFamily, z, V, T, i, λ̄)
    ϵ̄_ = @f(ϵ̄,i)
    σ̄_ = @f(σ̄,i)
    d̄_ = @f(d̄,i)
    x̄_0 = σ̄_/d̄_
    ζ_X_ = @f(ζ_X)
    I = (1-x̄_0^(3-λ̄))/(λ̄-3)
    J = (1-x̄_0^(4-λ̄)*(λ̄-3)+x̄_0^(3-λ̄)*(λ̄-4))/((λ̄-3)*(λ̄-4))
    return 2π*@f(ρ_S)*d̄_^3*ϵ̄_ * ((1-ζ_X_/2)/(1-ζ_X_)^3*I-9ζ_X_*(1+ζ_X_)/(2(1-ζ_X_)^3)*J)
end

function ∂B∂ρ_S(model::SAFTgammaMieFamily, z, V, T, i, λ̄)
    ϵ̄_ = @f(ϵ̄,i)
    σ̄_ = @f(σ̄,i)
    d̄_ = @f(d̄,i)
    x̄_0 = σ̄_/d̄_
    ζ_X_ = @f(ζ_X)
    I = (1-x̄_0^(3-λ̄))/(λ̄-3)
    J = (1-x̄_0^(4-λ̄)*(λ̄-3)+x̄_0^(3-λ̄)*(λ̄-4))/((λ̄-3)*(λ̄-4))
    return @f(B̄,i,λ̄) + 2π*d̄_^3*ϵ̄_ * ζ_X_*((3*(1-ζ_X_/2)*(1-ζ_X_)^2
            -0.5*(1-ζ_X_)^3)*I/(1-ζ_X_)^6-9*J*((1+2*ζ_X_)*(1-ζ_X_)^3+ζ_X_*(1+ζ_X_)*3*(1-ζ_X_)^2)/(2*(1-ζ_X_)^6))
end

function ā_1(model::SAFTgammaMieFamily, z, V, T, i)
    λ̄a_ = @f(λ̄a,i)
    λ̄r_ = @f(λ̄r,i)
    x̄_0 = @f(σ̄,i)/@f(d̄,i)
    return @f(C,λ̄a_,λ̄r_)*(x̄_0^λ̄a_*(@f(āS_1,λ̄a_)+@f(B̄,λ̄a_))-x̄_0^λ̄r_*(@f(āS_1,λ̄r_)+@f(B̄,λ̄r_)))
end

function ∂ā_1∂ρ_S(model::SAFTgammaMieFamily, z, V, T, i)
    λ̄a_ = @f(λ̄a,i)
    λ̄r_ = @f(λ̄r,i)
    x̄_0 = @f(σ̄,i)/@f(d̄,i)
    return @f(C,λ̄a_,λ̄r_)*(x̄_0^λ̄a_*(@f(∂āS_1∂ρ_S,i,λ̄a_)+@f(∂B∂ρ_S,i,λ̄a_))-x̄_0^λ̄r_*(@f(∂āS_1∂ρ_S,i,λ̄r_)+@f(∂B∂ρ_S,i,λ̄r_)))
end

function āS_1(model::SAFTgammaMieFamily, z, V, T, i, λ̄)
    ϵ̄_ = @f(ϵ̄,i)
    σ̄_ = @f(σ̄,i)
    d̄_ = @f(d̄,i)
    ζ̄eff_ = @f(ζeff, λ̄)
    return -2π*@f(ρ_S) * ϵ̄_*d̄_^3/(λ̄-3) * (1-ζ̄eff_/2)/(1-ζ̄eff_)^3
end

function ∂āS_1∂ρ_S(model::SAFTgammaMieFamily, z, V, T, i, λ̄)
    ϵ̄_ = @f(ϵ̄,i)
    σ̄_ = @f(σ̄,i)
    d̄_ = @f(d̄,i)
    ζ̄eff_ = @f(ζeff, λ̄)
    A = [ 0.81096    1.7888   -37.578   92.284;
          1.02050  -19.341    151.26  -463.50 ;
         -1.90570   22.845   -228.14   973.92 ;
          1.08850   -6.1962   106.98  -677.64  ]
    ζ_X_ = @f(ζ_X)
    ∂ζ̄eff∂ρ_S = A * [1; 1/λ̄; 1/λ̄^2; 1/λ̄^3] ⋅ [1; 2ζ_X_; 3ζ_X_^2; 4ζ_X_^3]
    return @f(āS_1,i,λ̄) - 2π*(ϵ̄_*d̄_^3)/(λ̄-3) * ((3*(1-ζ̄eff_/2)*(1-ζ̄eff_)^2-1/2*(1-ζ̄eff_)^3)/(1-ζ̄eff_)^6 * ∂ζ̄eff∂ρ_S*ζ_X_)
end

function γ_c(model::SAFTgammaMieFamily, z, V, T, i)
    λ̄a_ = @f(λ̄a,i)
    λ̄r_ = @f(λ̄r,i)
    ϵ̄_ = @f(ϵ̄,i)
    ᾱ = @f(C,λ̄a_,λ̄r_)*(1/(λ̄a_-3)-1/(λ̄r_-3))
    θ = exp(ϵ̄_/T)-1
    ζst_X_ = @f(ζst_X)
    return 10 * (tanh(10*(0.57-ᾱ))+1) * ζst_X_*θ*exp(-6.7*ζst_X_-8ζst_X_)
end

function gMCA_2(model::SAFTgammaMieFamily, z, V, T, i)
    ϵ̄_ = @f(ϵ̄,i)
    σ̄_ = @f(σ̄,i)
    d̄_ = @f(d̄,i)
    x̄_0 = σ̄_/d̄_
    λ̄a_ = @f(λ̄a,i)
    λ̄r_ = @f(λ̄r,i)
    ρ_S_ = @f(ρ_S)
    C̄_ = @f(C,λ̄a_,λ̄r_)
    KHS_ = @f(KHS)
    return 1/(2π*ϵ̄_*d̄_^3)*(3*@f(∂ā_2∂ρ_S,i)-
        ϵ̄_*KHS_*C̄_^2*λ̄r_*x̄_0^(2λ̄r_)*(@f(āS_1,i,2λ̄r_)+@f(B̄,i,2λ̄r_))/ρ_S_+
        ϵ̄_*KHS_*C̄_^2*(λ̄a_+λ̄r_)*x̄_0^(λ̄a_+λ̄r_)*(@f(āS_1,i,λ̄a_+λ̄r_)+@f(B̄,i,λ̄a_+λ̄r_))-
        ϵ̄_*KHS_*C̄_^2*λ̄a_*x̄_0^(2λ̄a_)*(@f(āS_1,i,2λ̄a_)+@f(B̄,i,2λ̄a_))/ρ_S_)
end

function ∂ā_2∂ρ_S(model::SAFTgammaMieFamily, z, V, T, i)
    ϵ̄_ = @f(ϵ̄,i)
    σ̄_ = @f(σ̄,i)
    d̄_ = @f(d̄,i)
    x̄_0 = σ̄_/d̄_
    ζ_X_ = @f(ζ_X)
    ∂KHS∂ρ_S = -(4*(1-ζ_X_)^3*(1+4ζ_X_+4ζ_X_^2-4ζ_X_^3+ζ_X_^4) + (1-ζ_X_)^4*(4+8ζ_X_-12ζ_X_^2+4ζ_X_^3))/(1+4ζ_X_+4ζ_X_^2-4ζ_X_^3+ζ_X_^4)^2*ζ_X_/@f(ρ_S)
    λ̄a_ = @f(λ̄a,i)
    λ̄r_ = @f(λ̄r,i)
    return 1/2*ϵ̄_*@f(C,λ̄a_,λ̄r_)*(∂KHS∂ρ_S*(x̄_0^(2λ̄a_)*(@f(āS_1,i,2λ̄a_)+@f(B̄,i,2λ̄a_))
        -2x̄_0^(λ̄a_+λ̄r_)*(@f(āS_1,i,λ̄a_+λ̄r_)+@f(B̄,i,λ̄a_+λ̄r_))
        +x̄_0^(2λ̄r_)*(@f(āS_1,i,2λ̄r_)+@f(B̄,i,2λ̄r_)))
        +@f(KHS)*(x̄_0^(2λ̄a_)*(@f(∂āS_1∂ρ_S,i,2λ̄a_)+@f(∂B∂ρ_S,i,2λ̄a_))
        -2x̄_0^(λ̄a_+λ̄r_)*(@f(∂āS_1∂ρ_S,i,λ̄a_+λ̄r_)+@f(∂B∂ρ_S,i,λ̄a_+λ̄r_))
        +x̄_0^(2λ̄r_)*(@f(∂āS_1∂ρ_S,i,2λ̄r_)+@f(∂B∂ρ_S,i,2λ̄r_))))
end
