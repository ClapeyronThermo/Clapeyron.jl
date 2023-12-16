#TODO: find a better name

struct SAFTVRSMieParam <: EoSParam
    Mw::SingleParam{Float64}
    segment::SingleParam{Float64}
    sigma::PairParam{Float64}
    lambda_a::PairParam{Float64}
    lambda_r::PairParam{Float64}
    epsilon::PairParam{Float64}
end

abstract type SAFTVRSMieModel <: SAFTModel end
@newmodel SAFTVRSMie SAFTVRSMieModel SAFTVRSMieParam
default_references(::Type{SAFTVRSMie}) = ["TODO"]
default_locations(::Type{SAFTVRSMie}) = ["SAFT/SAFTVRMie", "properties/molarmass.csv"]

export SAFTVRSMie

function transform_params(::Type{SAFTVRSMie},params)
    sigma = params["sigma"]
    sigma.values .*= 1E-10
    sigma = sigma_LorentzBerthelot(sigma)
    epsilon = epsilon_HudsenMcCoubrey(params["epsilon"], sigma)
    lambda_a = lambda_LorentzBerthelot(params["lambda_a"])
    lambda_r = lambda_LorentzBerthelot(params["lambda_r"])
    params["sigma"] = sigma
    params["epsilon"] = epsilon
    params["lambda_a"] = lambda_a
    params["lambda_r"] = lambda_r
    return params
end

function a_res(model::SAFTVRSMieModel,V,T,z)
    _data = @f(data)
    return @f(a_mono,_data) + @f(a_chain,_data)
end

function ζ3(model::SAFTVRSMieModel, V, T, z,_d=@f(d),m̄ = dot(z,model.params.segment.values))
    m = model.params.segment
    ζ3 = zero(V+T+first(z)+one(eltype(model)))
    for i ∈ 1:length(z)
        di =_d[i]
        xS = z[i]*m[i]/m̄
        ζ3 += xS*di*di*di
    end
    c = π/6*N_A*m̄/V
    ζ3 = c*ζ3
    return ζ3
end

function data(model::SAFTVRSMieModel,V,T,z)
    ∑z = sum(z)
    m̄ = dot(model.params.segment.values,z)
    _d = d(model,V,T,z)
    η = ζ3(model,V,T,z,_d,m̄)
    ρS = m̄*N_A/V #is this ok?
    ηc = 0.740480489693061
    d̄ = dot(z,_d)/sum(z)
    ρ0 = 6*ηc/(π*d̄*d̄*d̄)
    γ = 4*(1 - ρS/ρ0)
    α = ρ0/ρS - 1
    Z = Zhs_hall(γ,α)
    J̄ = J(model,V,T,z,_d,η,Z)
    return m̄,_d,η,ρS,ρ0,Z,J̄
end

function a_mono(model::SAFTVRSMieModel,V,T,z,_data = @f(data))
    m̄,d,η,ρS,ρ0,Z,J̄ = _data
    ahs = @f(a_hs,_data)
    a1 = @f(a_1,_data)
    return m̄*(ahs + a1)/sum(z)
end

function a_hs(model::SAFTVRSMieModel,V,T,z,_data = @f(data))
    m̄,d,η,ρS,ρ0,Z,J̄ = _data
    ρ0ρs = ρ0/ρS
    γ = 4*(1 - ρS/ρ0)
    α = ρ0ρs - 1
    S₀ = -0.24 #± 0.04
    res = -S₀ + 1 + log(ρ0ρs) - 3*log(2*α/3)
    #=
    #Zhs(γ) - Zref = evalpoly(γ,(-0.442304,0.1253077,0.1762393,-1.053308,2.818621,-2.921934,1.118413))
    #integral(0,γ,evalpoly(γ,pol)/(γ - 4)) = integral(polynomial part) + integral(fractional part)

    #evalpoly(γ,pol)/(γ - 4) =
    1.11841*x^5 + 1.55172*x^4 + 9.02549*x^3 + 35.0487*x^2 + 140.371*x + 561.609 + 2245.99/(x - 4)
    (wolfram alpha)
    p_int = (561.609,140.371,35.0487,9.02549,1.55172,1.11841)
    = evalpoly(x,p_int) +  2245.99/(x - 4)
    =#
    p_int = (2.557696-3,0.1253077,0.1762393,-1.053308,2.818621,-2.921934,1.11841)

    ∫Z = zero(V+T+first(z))
    ∫Z += (log(1-γ/4))*p_int[1]
    ∫Z += (γ+4*log(1-γ/4))*p_int[2]
    ∫Z += (γ^2/2+4*γ+16*log(1-γ/4))*p_int[3]
    ∫Z += (γ^3/3+2*γ^2+16*γ+64*log(1-γ/4))*p_int[4]
    ∫Z += (γ^4/4+4*γ^3/3+8*γ^2+64*γ+256*log(1-γ/4))*p_int[5]
    ∫Z += (γ^5/5+γ^4+16*γ^3/3+32*γ^2+256*γ+1024*log(1-γ/4))*p_int[6]
    ∫Z += (γ^6/6+4*γ^5/5+4*γ^4+64*γ^3/3+128*γ^2+1024*γ+4096*log(1-γ/4))*p_int[7]
    return res + ∫Z
end

function Zhs_hall(γ,α) #eq 5
    return 3/α + evalpoly(γ,(2.557696,0.1253077,0.1762393,-1.053308,2.818621,-2.921934,1.118413))
end

function g_hs(model::SAFTVRSMieModel,η,d,r,J,r₁d = r1_d(η),k₁ = K1(model,η),k₂ = K2(model,η),k = K(model,η)) #eq 11
    g_hs₁ = g_hs_1(model,η,d,r,J,r₁d,k₁,k₂)
    g_hsᵢ = g_hs_i(model,η,d,r,k)
    return g_hs₁ + g_hsᵢ
end

function g_hs_fdf(model::SAFTVRSMieModel,V,T,z,d,r,i::Int) #eq 11, used in the context of the evaluation of d
    mᵢ = model.params.segment.values[i]
    η = (π/6*N_A*mᵢ*z[i]/V)*d*d*d
    ηc = 0.740480489693061
    ρs = z[i]*mᵢ*N_A/V #is this ok?
    ρ0 = 6*ηc/(π*d*d*d)
    γ = 4*(1 - ρs/ρ0)
    α = ρ0/ρs - 1
    Z = Zhs_hall(γ,α)
    r₁d = r1_d(η)
    k₁ = K1(model,η)
    k₂ = K2(model,η)
    k = K(model,η)
    J = g_hs_Ji(model,d,η,Z,k₁,k₂,r₁d,k)
    g_hs₁,∂g_hs₁ = g_hs_1_fdf(model,η,d,r,J,r₁d,k₁,k₂)
    g_hsᵢ,∂g_hsᵢ = g_hs_i_fdf(model,η,d,r,k)

    # @show Z
    # @show k₁,k₂,k
    # @show J
    # @show g_hs₁,∂g_hs₁
    # @show g_hsᵢ,∂g_hsᵢ
    return g_hs₁ + g_hsᵢ, ∂g_hs₁ + ∂g_hsᵢ
end

#eq 12
function g_hs_1(model::SAFTVRSMieModel,η,d,r,J,r1d = r1_d(η),k₁ = K1(model,η),k₂ = K2(model,η)) 
    rd = r/d
    _k12 = -(k₁*(rd-r1d))^2
    _k24 = -(k₂*(rd-r1d))^4
    (J*d/r)*exp(_k12 + _k24)
end

#eq 12 + derivative, used in calculation of d
function g_hs_1_fdf(model::SAFTVRSMieModel,η,d,r,J,r1d = r1_d(η),k₁ = K1(model,η),k₂ = K2(model,η))
    rd = r/d
    dr = d/r
    δr = (rd-r1d)
    _k12 = -(k₁*δr)^2
    _k24 = -(k₂*δr)^4
    _exp = exp(_k12+_k24)
    g = (J*d/r)*_exp
    ∂g = -g*dr*(dr+4*k₂^4*δr^3+2*k₁^2*δr)
    return g,∂g
end

function g_hs_i(model::SAFTVRSMieModel,η,d,r,k = K(model::SAFTVRSMieModel,η)) #eq 13
    n = SAFTVRSMieConsts.r_fcc::Vector{Int}
    g = zero(η+d+r+k)
    rd = r/d
    dr = d/r
    v₀ = π*d^3/(6*η)
    d₀ = cbrt(sqrt(2)*v₀)
    for rᵢd₀2 in 2:length(n) #(rᵢ/d₀)^2
        nᵢ = n[rᵢd₀2]
        nᵢ != 0 && begin
            rᵢ = sqrt(rᵢd₀2)*d₀
            δr = rd - rᵢ/d
            _expᵢ = exp(-(k*δr)^2)
            C = (d/rᵢ)*k*nᵢ/(sqrt(π)*24*η)
            g += dr*C*_expᵢ
        end
    end
    return g
end

function g_hs_i_fdf(model::SAFTVRSMieModel,η,d,r,k = K(model::SAFTVRSMieModel,η)) #eq 13
    n = SAFTVRSMieConsts.r_fcc::Vector{Int}
    g = zero(η+d+r+k)
    ∂g = zero(g)
    rd = r/d
    dr = d/r
    v₀ = π*d^3/(6*η)
    d₀ = cbrt(sqrt(2)*v₀)
    for rᵢd₀2 in 2:length(n) #(rᵢ/d₀)^2
        nᵢ = n[rᵢd₀2]
        nᵢ != 0 && begin
            rᵢ = sqrt(rᵢd₀2)*d₀
            δr = rd - rᵢ/d
            _expᵢ = exp(-(k*δr)^2)
            C = (d/rᵢ)*k*nᵢ/(sqrt(π)*24*η) #independent of d/r
            g += dr*C*_expᵢ
            ∂g += -g*(dr + 2*k*k*δr)*dr
        end
    end
    return g,∂g
end

function r1_d(η) #SA, eq 1
    ηc = 0.740480489693061
    η✷ = ηc - η
    num = evalpoly(η✷,(1.0,-8.0521,18.003))
    denom = evalpoly(η✷,(1.0,-8.2973,20.546,-13.828,103.95,-582.74,1245.7))
    return num/denom
end

function J(model::SAFTVRSMieModel,V,T,z,_d = @f(d),η = @f(ζ3,_d),Z = @f(Zhs_hall))
    #=
    g_hs(1,η) = (z-1)/(4*η)
    g_hs_1(1,η) + g_hs_i(1,η) = (z-1)/(4*η)
    g_hs_1(1,η) = (z-1)/(4*η) - g_hs_i(1,η)
    J*ghs_1_divJ(1,η) =  (z-1)/(4*η) - g_hs_i(1,η)
    J = 
    =#
    J̄ = fill(zero(V+T+first(z)+one(eltype(model))),length(model))

    for i in @comps
        r1d = r1_d(η)
        k₁ = K1(model,η)
        k₂ = K2(model,η)
        k = K(model,η)
        J̄[i] = g_hs_Ji(model,_d[i],η,Z,k₁,k₂,r1d,k)
    end
    return J̄
end

#specialization for single component, one less allocation
function J(model::SAFTVRSMieModel,V,T,z::SingleComp,_d = @f(d),η = @f(ζ3,_d),Z = @f(Zhs_hall))
    r1d = r1_d(η)
    k₁ = K1(model,η)
    k₂ = K2(model,η)
    k = K(model,η)
    return SA[g_hs_Ji(model,_d[1],η,Z,k₁,k₂,r1d,k)]
end

function g_hs_Ji(model::SAFTVRSMieModel,di,η,Z,k₁,k₂,r1d,k)
    ghs_at_1 = g_hs_i(model,η,di,di,k)
    _k12 = -(k₁*(1-r1d))^2
    _k24 = -(k₂*(1-r1d))^4
    ghs_1_divJ = exp(_k12+_k24)
    # @show ghs_at_1,ghs_1_divJ

    return ((Z-1)/(4*η) - ghs_at_1)/ghs_1_divJ
end

function K1(model::SAFTVRSMieModel,η)   #eq 15
    ηc = 0.740480489693061
    η✷ = ηc - η
    return 1.5338/η✷ - 0.37687*exp(-989.6*(η - 0.52)^2) + evalpoly(η✷,(-2.5146,-1.3574,-8.5038))
end

function K2(model::SAFTVRSMieModel,η)  #eq 16
    ηc = 0.740480489693061
    η✷ = ηc - η
    return 0.80313/η✷ - 1.208*exp(5.6128*η✷) + 67.808*η✷*η✷ - 67.918*η✷*η✷*η✷
end

function K(model::SAFTVRSMieModel,η)  #eq 17
    ηc = 0.740480489693061
    η✷ = ηc - η
    return 1.9881/η✷ + evalpoly(η✷,(-3.5276,6.9762,-26.205))
end

function a_1(model::SAFTVRSMieModel,V,T,z,_data = @f(data))
    m̄,d,η,ρS,ρ0,Z,J̄ = _data
    m̄inv = 1/m̄
    m = model.params.segment.values
    σ = model.params.sigma.values
    λa = model.params.lambda_a.values
    λr = model.params.lambda_r.values
    ϵ = model.params.epsilon.values
    a₁ = zero(V+T+first(z)+one(eltype(model)))
    for i ∈ @comps
        ρS = z[i]*m[i]*N_A/V
        λ = cbrt(sqrt(2)/ρS)*one(T)
        λ̄ = λ/σ[i]
        
        xsᵢ = z[i]*m[i]*m̄inv
        dᵢᵢ = d[i]/σ[i]
        ϵ̄ = ϵ[i]/T
        r₁dᵢ = r1_d(η)
        k₁ᵢ = K1(model,η)
        k₂ᵢ = K2(model,η)
        kᵢ = K(model,η)
        Cᵢ = Cλ_mie(λa[i], λr[i])
        u(r) = Cᵢ*ϵ̄*(r^-λr[i] -r^-λa[i])
        du(r) = -Cᵢ*ϵ̄*(λr[i]*r^-(λr[i]+1) -λa[i]*r^-(λa[i]+1))


        function W(r)
            if r >= λ̄
                return u(r)
            else
                uλ = u(λ̄)
                duλ = du(λ̄)
                #du/dr at r = λ
                return uλ - duλ*(λ̄ - r)
            end
        end


        ghsWrᵢ = r -> begin 
            _Wr = W(r)
            g_hs(model,η,dᵢᵢ*σ[i],r*σ[i],J̄[i],r₁dᵢ,k₁ᵢ,k₂ᵢ,kᵢ)*r*r*_Wr
        end
        Wrᵢ = r -> begin
            W(r)*r*r
        end
        #we separate the integration between [d,3.3d] and [3.3d,∞]
        #the second part can be solved more efficiently.
        a₁ᵢᵢ = Solvers.integral64(ghsWrᵢ,dᵢᵢ,3.3*dᵢᵢ)*σ[i]^3 #TODO is this grade ok? maybe a 10-point quadrature would suffice
        a₁ᵢᵢ += Solvers.integral64(Wrᵢ,3.3dᵢᵢ,10*dᵢᵢ)*σ[i]^3
        a₁ += a₁ᵢᵢ*xsᵢ*xsᵢ
        for j ∈ 1:(i-1)
            xsⱼ = z[i]*m[i]*m̄inv
            dᵢⱼ = 0.5*(dᵢᵢ+d[j])
            a₁ᵢⱼ = 1 #TODO: what is the definition of g_hs here?
            a₁ += 2*a₁ᵢⱼ*xsᵢ*xsⱼ
        end
    end
    return 2*π*ρS*a₁
end


function a_chain(model::SAFTVRSMieModel,V,T,z,_data = @f(data))
    m̄,d,η,ρS,ρ0,Z,J̄ = _data
    β = 1/T #our epsilon is already divided by kᵦ
    σ = model.params.sigma
    λa = model.params.lambda_a
    λr = model.params.lambda_r
    ϵ = model.params.epsilon
    achain = zero(V+T+first(z)+one(eltype(model)))
    m = model.params.segment.values
    # eq 8, using linear mixing (from SAFTVRMie)
    λ = cbrt(sqrt(2)/ρS)*one(T)
    for i in 1:length(model)
        σᵢ,λaᵢ,λrᵢ,ϵᵢ,dᵢ = σ[i],λa[i],λr[i],ϵ[i],d[i]
        Cᵢ = Cλ_mie(λaᵢ, λrᵢ)
        λ̄ = λ/σᵢ
        ϵ̄ = β*ϵᵢ
        uλ = Cᵢ*ϵ̄*(λ̄^-λrᵢ -λ̄^-λaᵢ)
        duλ = -Cᵢ*ϵ̄*(λrᵢ*λ̄^-(λrᵢ+1) -λaᵢ*λ̄^-(λaᵢ+1))
        
        g_hsᵢ = g_hs(model,η,dᵢ,σᵢ,J̄[i])
        βV₀ᵢ = - uλ + duλ*(λ̄ - 1)
        
        y_hsᵢ = g_hsᵢ
        g_Mieᵢ = y_hsᵢ*exp(-βV₀ᵢ)
        achain -= z[i]*(log(g_Mieᵢ)*(m[i] - 1))
    end
    return achain/sum(z)
end

function d(model::SAFTVRSMieModel,V,T,z)
    n = length(z)

    _d = fill(zero(V+T+first(z)+one(eltype(model))),n)
    for k ∈ 1:n
        _d[k] = d_vrs(model,V,T,z,k)
    end
    return _d
end

function d(model::SAFTVRSMieModel,V,T,z::SingleComp)
    return SA[d_vrs(model::SAFTVRSMieModel,V,T,z,1)]
end

#TODO: this is a sketch
function d_vrs(model::SAFTVRSMieModel,V,T,z,i::Int)
    _0 = zero(T+first(z))
    mᵢ = model.params.segment.values[i]
    ρS = z[i]*mᵢ*N_A/V
    λ = cbrt(sqrt(2)/ρS)*one(T)
    σᵢ = model.params.sigma.values[i,i]
    λ̄ = λ/σᵢ

    λa,λr = model.params.lambda_a[i],model.params.lambda_r[i]
    Cᵢ = Cλ_mie(λa, λr)
    β = 1/T
    ϵ = model.params.epsilon[i]
    ϵ̄ = ϵ/T

    u(r) = Cᵢ*ϵ̄*(r^-λr -r^-λa)
    du(r) = -Cᵢ*ϵ̄*(λr*r^-(λr+1) -λa*r^-(λa+1))

    function V₀(r)
        if r > λ̄
            return zero(r)
        else
            uλ = u(λ̄)
            duλ = du(λ̄)
            #du/dr at r = λ
            return u(r) - uλ + duλ*(λ̄ - r)
        end
    end

    function dV₀(r)
        if r > λ̄
            return zero(r)
        else
            duλ = du(λ̄)
            return du(r) - duλ
        end
    end

    dB_f(r) = -expm1(-V₀(r))
    dB = Solvers.integral64(dB_f,_0,λ̄)
    # dB *= σᵢ
    #TODO: derivate this
    δ_f(r)  = -dV₀(r)*exp(-V₀(r))*(r/dB -1)^2
    δ = Solvers.integral64(δ_f,_0,λ̄)
    #iteration 0: d = dB
    d = dB
    d0 = one(d)
    function f0(dj)
        g_hs,∂g_hs = g_hs_fdf(model,V,T,z,dj*σᵢ,dj*σᵢ,i)
        σ₀ = g_hs
        σ₁ = 2*σ₀ + ∂g_hs
        return dB*(1 + δ*σ₁/(2*σ₀))
    end

    return Solvers.fixpoint(f0,d0)*σᵢ
    #=
    while abs(d - d0) > 1E-12 && k < 100
        d0 = d
        g_hs,∂g_hs = g_hs_fdf(model,V,T,z,d0*σᵢ,d0*σᵢ,i)
        # #=TODO 
        # is d > λ ? if so, y_hs = g_hs and it would simplify this calculation a lot
        # =#
        # V₀ = d > λ ? zero(λ) : 
        # expV₀ = exp(-V₀(d)/T)
        # ∂expV₀ = 1 #TODO
        # y_hs = g_hs*expV0
        # dy_hs = g_hs*∂expV₀ + ∂g_hs*expV₀
        σ₀ = g_hs
        σ₁ = 2*σ₀ + ∂g_hs
        d = dB*(1 + δ*σ₁/(2*σ₀))
        k += 1
    end
    return d*σᵢ =#
end

#SA, Table 1
#note: idx = [14,30,46,56,62] are missing.
#for purposes of wiping out an MVE, those values are replaced with 1
SAFTVRSMieConsts = (;
    r_fcc = [12, 6, 24, 12, 24, 8, 48, 6, 36, 24, 24, 24, 72, 0, 48, 12, 48, 30, 72, 24, 48, 24, 48, 8, 84, 24, 96, 48, 24, 0, 96, 6, 96, 48, 48, 36, 120, 24, 48, 24, 48, 48, 120, 24, 120, 0, 96, 24, 108, 30, 48, 72, 72, 32, 144, 0, 96, 72, 72, 48, 120, 0, 144, 12, 48],
    )