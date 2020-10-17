function a_res(model::SAFTgammaMieFamily, z, V, T)
    return a_mono(model,z,V,T) + a_chain(model,z,V,T) + a_assoc(model,z,V,T)
end

function a_mono(model::SAFTgammaMieFamily, z, V, T)
    return a_hs(model,z,V,T) + a1(model,z,V,T) + a2(model,z,V,T) + a3(model,z,V,T)
end

function a_chain(model::SAFTgammaMieFamily, z, V, T)
    return 0
end

function a_assoc(model::SAFTgammaMieFamily, z, V, T)
    return 0
end

function a_hs(model::SAFTgammaMieFamily, z, V, T)
    ζ0   = @f(ζ,0)
    ζ1   = @f(ζ,1)
    ζ2   = @f(ζ,2)
    ζ3   = @f(ζ,3)
    ρ = sum(z)*N_A/V
    return 6/π/ρ * (3ζ1*ζ2/(1-ζ3) + ζ2^3/(ζ3*(1-ζ3)^2) + (ζ2^3/ζ3^2-ζ0)*log(1-ζ3))
end

#= for n in 1:3 =#
#=     @eval   function $(Symbol(:a,Symbol(n)))(model::SAFTgammaMieFamily, z, V, T) =#
#=                 x = z/sum(z) =#
#=                 v = model.group_multiplicities =#
#=                 v_st = model.params.segment =#
#=                 S = model.params.shapefactor =#

#=                 return 1/(kB*T)^n * sum(x[i] * sum(v[i][k]*v_st[k]*S[k] for k in @groups(i)) for i in @comps) * @eval $(@f(Symbol(:â,Symbol(n)))) =#
#=             end =#

#=     @eval   function $(Symbol(:â,Symbol(n)))(model::SAFTgammaMieFamily, z, V, T) =#
#=                 return sum(@f(xS,k)*@f(xS,l) * @eval $(@f(Symbol(:â,Symbol(n)),k,l)) for k in @groups for l in @groups) =#
#=             end =#
#= end =#
#
#= function â(model::SAFTgammaMieFamily, z, V, T) =#

function a1(model::SAFTgammaMieFamily, z, V, T)
    x = z/sum(z)
    v = model.group_multiplicities
    v_st = model.params.segment
    S = model.params.shapefactor

    return 1/T * sum(x[i]*sum(v[i][k]*v_st[k]*S[k] for k in @groups(i)) for i in @comps) * @f(â1)
end
function a2(model::SAFTgammaMieFamily, z, V, T)
    x = z/sum(z)
    v = model.group_multiplicities
    v_st = model.params.segment
    S = model.params.shapefactor

    return 1/T^2 * sum(x[i]*sum(v[i][k]*v_st[k]*S[k] for k in @groups(i)) for i in @comps) * @f(â2)
end
function a3(model::SAFTgammaMieFamily, z, V, T)
    x = z/sum(z)
    v = model.group_multiplicities
    v_st = model.params.segment
    S = model.params.shapefactor

    return 1/T^3 * sum(x[i]*sum(v[i][k]*v_st[k]*S[k] for k in @groups(i)) for i in @comps) * @f(â3)
end

function â1(model::SAFTgammaMieFamily, z, V, T)
    return sum(@f(xS,k)*@f(xS,l)*@f(â1,k,l) for k in @groups for l in @groups)
end
function â2(model::SAFTgammaMieFamily, z, V, T)
    return sum(@f(xS,k)*@f(xS,l)*@f(â2,k,l) for k in @groups for l in @groups)
end
function â3(model::SAFTgammaMieFamily, z, V, T)
    return sum(@f(xS,k)*@f(xS,l)*@f(â3,k,l) for k in @groups for l in @groups)
end

function â1(model::SAFTgammaMieFamily, z, V, T, k, l)
    σ = model.params.sigma[union(k,l)]
    λ_a = model.params.lambda_r[union(k,l)]
    λ_r = model.params.lambda_a[union(k,l)]

    x0 = σ/@f(d,k,l)
    return @f(C,k,l) * (x0^λ_a*(@f(â1_S,k,l,λ_a)+@f(B,k,l,λ_a)) - x0^λ_r*(@f(â1_S,k,l,λ_r)+@f(B,k,l,λ_r)))
end
function â2(model::SAFTgammaMieFamily, z, V, T, k, l)
    σ = model.params.sigma[union(k,l)]
    ϵ = model.params.epsilon[union(k,l)]
    λ_a = model.params.lambda_r[union(k,l)]
    λ_r = model.params.lambda_a[union(k,l)]

    x0 = σ/@f(d,k,l)
    return 1/2*@f(K_hs)*(1+@f(𝜒,k,l))*@f(C,k,l)^2 * (
           x0^(2λ_a)*(@f(â1_S,k,l,2λ_a) + @f(B,k,l,2λ_a))
         - 2x0^(λ_a+λ_r)*(@f(â1_S,k,l,λ_a+λ_r) + @f(B,k,l,λ_a+λ_r))
         + x0^(2λ_r)*(@f(â1_S,k,l,2λ_r) + @f(B,k,l,2λ_r)) )
end
function â3(model::SAFTgammaMieFamily, z, V, T, k, l)
    ϵ = model.params.epsilon[union(k,l)]

    ζX_st_ = @f(ζX_st)
    return -ϵ^3*@f(f,k,l,4)*ζX_st_ * exp(@f(f,k,l,5)*ζX_st_ + @f(f,k,l,6)*ζX_st_^2)
end

function B(model::SAFTgammaMieFamily, z, V, T, k, l, λ)
    ϵ = model.params.epsilon[union(k,l)]
    ζX_ = @f(ζX)
    return 2π*@f(ρS)*@f(d,k,l)^3*ϵ * ((1-ζX_/2)/(1-ζX_)^3*@f(I,k,l,λ)-9ζX_*(1+ζX_)/(2(1-ζX_)^3)*@f(J,k,l,λ))
end

function ζX(model::SAFTgammaMieFamily, z, V, T)
    return π*@f(ρS)/6 * sum(@f(xS,k)*@f(xS,l)*@f(d,k,l)^3 for k in @groups for l in @groups)
end

function ζX_st(model::SAFTgammaMieFamily, z, V, T)
    σ = model.params.sigma
    return π*@f(ρS)/6 * sum(@f(xS,k)*@f(xS,l)*σ[union(k,l)]^3 for k in @groups for l in @groups)
end

function I(model::SAFTgammaMieFamily, z, V, T, k, l, λ)
    σ = model.params.sigma[union(k,l)]
    x0 = σ/@f(d,k,l)
    return (1-x0^(3-λ))/(λ-3)
end

function J(model::SAFTgammaMieFamily, z, V, T, k, l, λ)
    σ = model.params.sigma[union(k,l)]
    x0 = σ/@f(d,k,l)
    return (1-x0^(4-λ)*(λ-3)+x0^(3-λ)*(λ-4))/((λ-3)*(λ-4))
end

function â1_S(model::SAFTgammaMieFamily, z, V, T, k, l, λ)
    ϵ = model.params.epsilon[union(k,l)]
    ζ_eff_ = @f(ζ_eff, λ)
    return -2π*@f(ρS) * ϵ*@f(d,k,l)^3/(λ-3) * (1-ζ_eff_/2)/(1-ζ_eff_)^3
end

function ζ_eff(model::SAFTgammaMieFamily, z, V, T, λ)
    A = [ 0.81096    1.7888   -37.578   92.284;
          1.02050  -19.341    151.26  -463.50 ;
         -1.90570   22.845   -228.14   973.92 ;
          1.08850   -6.1962   106.98  -677.64  ]
    ζX_ = @f(ζX)
    return A * [1; 1/λ; 1/λ^2; 1/λ^3] ⋅ [ζX_; ζX_^2; ζX_^3; ζX_^4]
end

function K_hs(model::SAFTgammaMieFamily, z, V, T)
    ζX_ = @f(ζX)
    return (1-ζX_)^4/(1+4ζX_+4ζX_^2-4ζX_^3+ζX_^4)
end

function 𝜒(model::SAFTgammaMieFamily, z, V, T, k, l)
    ζX_st_ = @f(ζX_st)
    return @f(f,k,l,1)*ζX_st_ + @f(f,k,l,2)*ζX_st_^5 + @f(f,k,l,3)*ζX_st_^8
end

function f(model::SAFTgammaMieFamily, z, V, T, k, l, m)
    ϕ = [[7.5365557, -359.440,  1550.9, -1.199320, -1911.2800,  9236.9,  10.0],
         [-37.604630,  1825.60, -5070.1,  9.063632,  21390.175, -129430,  10.0],
         [71.745953, -3168.00,  6534.6, -17.94820, -51320.700,  357230,  0.57],
         [-46.835520,  1884.20, -3288.7,  11.34027,  37064.540, -315530, -6.70],
         [-2.4679820,- 0.82376, -2.7171,  20.52142,  1103.7420,  1390.2, -8.00],
         [-0.5027200, -3.19350,  2.0883, -56.63770, -3264.6100, -4518.2,   NaN],
         [8.0956883,  3.70900,  0.0000,  40.53683,  2556.1810,  4241.6,   NaN]];
    λ_a = model.params.lambda_r[union(k,l)]
    λ_r = model.params.lambda_a[union(k,l)]

    α = @f(C,k,l)*(1/(λ_a-3)-1/(λ_r-3))
    return sum(ϕ[i+1][m]*α^i for i in 0:3)/(1+sum(ϕ[i+1][m]*α^(i-3) for i in 4:6))
end

function ζ(model::SAFTgammaMieFamily, z, V, T, m)
    return π/6*@f(ρS)*sum(@f(xS,k)*@f(d,k)^m for k in @groups)
end

function ρS(model::SAFTgammaMieFamily, z, V, T)
    x = z/sum(z)
    v = model.group_multiplicities
    v_st = model.params.segment
    S = model.params.shapefactor

    ρ = sum(z)*N_A/V
    return ρ * sum(x[i] * sum(v[i][k]*v_st[k]*S[k] for k in @groups(i)) for i in @comps)
end

function xS(model::SAFTgammaMieFamily, z, V, T, k)
    x = z/sum(z)
    v = model.group_multiplicities
    v_st = model.params.segment
    S = model.params.shapefactor

    return sum(x[i]*v[i][k]*v_st[k]*S[k] for i in @comps) / sum(x[i] * sum(v[i][l]*v_st[l]*S[l] for l in @groups(i)) for i in @comps)
end

function d(model::SAFTgammaMieFamily, z, V, T, k)
    ϵ = model.params.epsilon[k]
    σ = model.params.sigma[k]
    λ_r = model.params.lambda_r[k]
    λ_a = model.params.lambda_a[k]
    u = [0.26356031971814109102031,1.41340305910651679221800,3.59642577104072208122300,7.08581000585883755692200,12.6408008442757826594300]
    w = [0.5217556105828086524759,0.3986668110831759274500,7.5942449681707595390e-2,3.6117586799220484545e-3,2.3369972385776227891e-5]
    θ = @f(C,k,k)*ϵ/T
    return σ*(1-sum(w[j]*(θ./(θ+u[j]))^(1/λ_r)*(exp(θ*(1/(θ./(θ+u[j]))^(λ_a/λ_r)-1))/(u[j]+θ)/λ_r) for j in 1:5))
end

function d(model::SAFTgammaMieFamily, z, V, T, k, l)
    # Lorentz mixing rule
    if k == l
        return @f(d, k)
    else
        return (@f(d, k) + @f(d, l))/2
    end
end

function C(model::SAFTgammaMieFamily, z, V, T, k, l)
    λ_r = model.params.lambda_r[union(k,l)]
    λ_a = model.params.lambda_a[union(k,l)]
    return (λ_r/(λ_r-λ_a)) * (λ_r/λ_a)^(λ_a/(λ_r-λ_a))
end
