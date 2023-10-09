#TODO: find a better name

struct VRSParam <: EoSParam
    Mw::SingleParam{Float64}
    segment::SingleParam{Float64}
    sigma::PairParam{Float64}
    lambda_a::PairParam{Float64}
    lambda_r::PairParam{Float64}
    epsilon::PairParam{Float64}
end

abstract type VRSModel <: SAFTModel end
@newmodel VRS VRSModel VRSParam
default_references(::Type{VRS}) = ["TODO"]
default_locations(::Type{VRS}) = ["SAFT/SAFTVRMie", "properties/molarmass.csv"]

function transform_params(::Type{VRS},params)
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

function a_res(model::VRSModel,V,T,z)
    _data = @f(data)
    return @f(a_mono,_data) + @f(a_chain,_data)
end

function ζ3(model::VRSModel, V, T, z,_d=@f(d),m̄ = dot(z,model.params.segment.values))
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

function data(model::VRSModel,V,T,z)
    ∑z = sum(z)
    m̄ = dot(model.params.segment.values,z)
    _d = d(model,V,T,z)
    η = ζ3(model,V,T,z,_d,m̄)
    ρS = N_A/V*m̄
    ρ0 = 1 #TODO
    γ = 4*(1 - ρs/ρ0)
    α = ρ0ρs - 1
    Z = Zhs_hall(γ,α)
    J̄ = J(model,V,T,z,_d,η,Z)
    return m̄,d,η,ρS,ρ0,Z,J̄
end

function a_mono(model::VRSModel,V,T,z,_data = @f(data))
    ahs = @f(a_hs,_data)
    a1 = @f(a_1,_data)
    return m̄*(ahs + a1)
end

function a_hs(model::VRSModel,V,T,z,_data = @f(data))
    m̄,d,η,ρS,ρ0,Z,J̄ = _data
    ρ0ρs = ρ0/ρs
    γ = 4*(1 - ρs/ρ0)
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
    p_int = (561.609,140.371,35.0487,9.02549,1.55172,1.11841)
    ∫0 = -3113.603272131663 #log(abs(x-4))*2245.99 at x = 0, the polynomial integral evaluated at 0 == 0
    ∫γ = evalpolyint(γ,p_int) + 2245.99*log(abs(γ-4))
    return res + ∫γ - ∫0
end

function Zhs_hall(γ,α) #eq 5
    return 3/α + evalpoly(γ,(2.557696,0.1253077,0.1762393,-1.053308,2.818621,-2.921934,1.118413))
end

function Zhs_hall(model::VRSModel,V,T,z)
    ρS = N_A/V*m̄
    ρ0 = 1 #TODO
    γ = 4*(1 - ρs/ρ0)
    α = ρ0ρs - 1
    return Zhs_hall(γ,α)
end





function g_hs(model,η,d,r,J) #eq 11
    r > 3.3*d && one(η+d+r+J) #mean-field correction
    g_hs₁ = g_hs_1(model,η,d,r,J)
    g_hsᵢ = g_hs_2_n(model,η,d,r)
    return g_hs₁ + g_hsᵢ
end

function g_hs_1(model::VRSModel,η,d,r,J) #eq 12
    k₁ = K1(model,η)
    k₂ = K2(model,η)
    _k12 = (-k₁*(1-r1d))^2
    _k24 = (-k₂*(1-r1d))^4
    (J*d/r)*exp(_k12+_k24)
end

function g_hs_2_n(model::VRSModel,η,d,r) #eq 13
    k = K(model,η)
    n = VRSConsts.r_fcc
    g = zero(typeof(η))
    rd = r/d
    dr = d/r
    v₀ = π*d^3/(6*η)
    d₀ = cbrt(sqrt(2)*v₀)
    for rᵢd₀2 in 2:length(g) #(rᵢ/d₀)^2
        rᵢ = sqrt(rᵢd₀2)*d₀
        nᵢ = n[rᵢd₀2]
        g += (dr*(d/rᵢ)*K*nᵢ*/(sqrt(π)*24*η))*exp(-(k*(rd - rᵢ/d))^2)
    end
    return g
end

function r1_d(η) #SA, eq 1
    ηc = 0.740480489693061
    η✷ = ηc - η
    num = evalpoly(η✷,(1.0,-8.0521,18.003))
    denom = evalpoly(η✷,(1.0,-8.2973,20.546,-13.828,103.95,-582.74,1245.7))
    return num/denom
end

function J(model::VRSModel,V,T,z,_d = @f(d),η = @f(ζ3,_d),Z = @f(Zhs_hall))
    #=
    g_hs(1,η) = (z-1)/(4*η)
    g_hs_1(1,η) + g_hs_2_n(1,η) = (z-1)/(4*η)
    g_hs_1(1,η) = (z-1)/(4*η) - g_hs_2_n(1,η)
    J*ghs_1_divJ(1,η) =  (z-1)/(4*η) - g_hs_2_n(1,η)
    J = 
    =#
    r1d = r1_d(η)
    k₁ = K1(model,η)
    k₂ = K2(model,η)
    J̄ = zeros(length(model),V+T+first(z))
    for i in @comps
        J̄[i] = g_hs_Ji(model,_d[i],η,Z,k₁,k₂,r1d)
    end
    return J̄
end

#specialization for single component, one less allocation
function J(model::VRSModel,V,T,z::SingleComp,_d = @f(d),η = @f(ζ3,_d),Z = @f(Zhs_hall))
    r1d = r1_d(η)
    k₁ = K1(model,η)
    k₂ = K2(model,η)
    return SA[g_hs_Ji(model,_d[i],η,Z,k₁,k₂,r1d)]
end

function g_hs_Ji(model::VRSModel,di,η,Z,k₁,k₂,r1d)
    ghs_at_1 = g_hs_2_n(model,η,di,di)
    _k12 = (-k₁*(1-r1d))^2
    _k24 = (-k₂*(1-r1d))^4
    ghs_1_divJ = exp(_k12+_k24)
    ((Z-1)/(4*η) - ghs_at_1)/ghs_1_divJ
end

function K1(model::VRSModel,η)   #eq 15
    ηc = 0.740480489693061
    η✷ = ηc - η
    return 1.5338/η✷ - 0.37687*exp(-989.6*(η - 0.52)^2) + evalpoly(η✷,(-2.5146,-1.3574,-8.5038))
end

function K2(model::VRSModel,η)  #eq 16
    ηc = 0.740480489693061
    η✷ = ηc - η
    return 0.80313/η✷ - 1.208*exp(5.6128*η✷) + 67.808*η✷*η✷ - 67.918*η✷*η✷*η✷
end

function K(model::VRSModel,η)  #eq 17
    ηc = 0.740480489693061
    η✷ = ηc - η
    return 1.9881/η✷ + evalpoly(η✷,(-3.5276,6.9762,-26.205))
end

function a_1(model::VRSModel,V,T,z,_data = @f(data))
    m̄,d,η,ρS,ρ0,Z,J̄ = _data
    m̄inv = 1/m̄
    a₁ = zero(V+T+first(z)+one(eltype(model)))
    for i ∈ @comps
        xsᵢ = z[i]*m[i]*m̄inv
        dᵢᵢ = d[i]
        ghsWrᵢ = r -> begin 
            _Wr = 1 #TODO: W(r)
            g_hs(model,η,dᵢᵢ,r,J̄[i])*r*r*_Wr
        end
        Wrᵢ = r -> begin
            one(dᵢᵢ)*r*r #TODO: W(r)
        end
        #we separate the integration between [d,3.3d] and [3.3d,∞]
        #the second part can be solved more efficiently.
        a₁ᵢᵢ = Solvers.integral21(ghsWrᵢ,dᵢᵢ,3.3*dᵢᵢ) #TODO is this grade ok? maybe a 10-point quadrature would suffice
        a₁ᵢᵢ += Solvers.integral21(Wrᵢ,3.3dᵢᵢ,10*3dᵢᵢ)
        a₁ += a₁ᵢᵢ*xsᵢ*xsᵢ
        for j ∈ 1:(i-1)
            xsⱼ = z[i]*m[i]*m̄inv
            dᵢⱼ = 0.5*(dᵢᵢ+d[j])
            a₁ᵢⱼ = 1 #TODO: what is the definition of g_hs here?
            a₁ += 2*a₁ᵢⱼ*xsᵢ*xsⱼ
        end
    end
    return 2*π*ρS*a₁/T
end

function Wri_integral(λa,λr,T)
     #TODO
end
function a_chain(model::VRSModel,V,T,z,_data = @f(data))
    m̄,d,η,ρS,ρ0,Z,J̄ = _data
    β = 1/T #our epsilon is already divided by kᵦ
    σ = model.params.sigma
    λa = model.params.lambda_a
    λr = model.params.lambda_r
    achain = zero(V+T+first(z)+one(eltype(model)))
    m = model.params.segment.values
    # eq 8, using linear mixing (from SAFTVRMie)
    for i in 1:length(model)
        σᵢ,λaᵢ,λrᵢ,dᵢ = σ[i],λa[i],λr[i],d[i]
        g_hsᵢ = g_hs(model,η,dᵢ,σᵢ,J̄[i])
        V_hsᵢ = 1 #TODO
        V₀ᵢ = 1#TODO
        y_hsᵢ = g_hsᵢ*exp(β*V_hsᵢ)
        g_Mieᵢ = y_hsᵢ*exp(β*V₀ᵢ)
        achain -= z[i]*(log(g_Mieᵢ)*(m[i] - 1))
    end
    return achain/sum(z)
end

function d(model::VRSModel,V,T,z)
    #TODO
end

#SA, Table 1
#note: idx = [14,30,46,56,62] are missing.
#for purposes of wiping out an MVE, those values are replaced with 1
VRSConsts = (;
    r_fcc = [12, 6, 24, 12, 24, 8, 48, 6, 36, 24, 24, 24, 72, 1, 48, 12, 48, 30, 72, 24, 48, 24, 48, 8, 84, 24, 96, 48, 24, 1, 96, 6, 96, 48, 48, 36, 120, 24, 48, 24, 48, 48, 120, 24, 120, 1, 96, 24, 108, 30, 48, 72, 72, 32, 144, 1, 96, 72, 72, 48, 120, 1, 144, 12, 48],
    )