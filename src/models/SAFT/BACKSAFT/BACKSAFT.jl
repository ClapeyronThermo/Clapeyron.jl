struct BACKSAFTParam <: EoSParam
    Mw::SingleParam{Float64}
    segment::SingleParam{Float64}
    sigma::PairParam{Float64}
    epsilon::PairParam{Float64}
    c::SingleParam{Float64}
    alpha::SingleParam{Float64}
end

abstract type BACKSAFTModel <: SAFTModel end
#BACKSAFT does not define association.
@newmodel BACKSAFT BACKSAFTModel BACKSAFTParam false
default_references(::Type{BACKSAFT}) = ["10.1016/s0378-3812(01)00521-0","10.1016/s0378-3812(02)00093-6","10.1016/j.fluid.2010.08.019"]
default_locations(::Type{BACKSAFT}) = ["SAFT/BACKSAFT","properties/molarmass.csv"]
function transform_params(::Type{BACKSAFT},params)
    k = get(params,"k",nothing)
    l = get(params,"l",nothing)
    sigma = params["vol"]
    sigma.values .*= 6/N_A/1e6/π
    sigma.values .^= 1/3
    epsilon = params["epsilon"]
    params["sigma"] = sigma_LorentzBerthelot(sigma)
    params["epsilon"] = epsilon_LorentzBerthelot(epsilon, k)
    return params
end

function get_k(model::BACKSAFT)   
    return get_k_geomean(model.params.epsilon)
end

function get_l(model::BACKSAFT)   
    return get_k_mean(model.params.sigma)
end

export BACKSAFT

"""
    BACKSAFTModel <: SAFTModel

    BACKSAFT(components; 
    idealmodel = BasicIdeal,
    userlocations = String[],
    ideal_userlocations = String[],
    reference_state = nothing,
    verbose = false,
    assoc_options = AssocOptions())

## Input parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `vol`: Single Parameter (`Float64`) - Segment Volume `[dm³]`
- `epsilon`: Single Parameter (`Float64`) - Reduced dispersion energy `[K·mol⁻¹]`
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Parameter (no units)
- `c`: Single Parameter (`Float64`) - Adjustable parameter (no units)
- `alpha`: Single Parameter (`Float64`) - Non-spherical deviation (no units)

## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Pair Parameter (`Float64`) - Mixed segment Diameter `[m]`
- `epsilon`: Pair Parameter (`Float64`) - Mixed reduced dispersion energy`[K]`
- `c`: Single Parameter (`Float64`) - Adjustable parameter (no units)
- `alpha`: Single Parameter (`Float64`) - Non-spherical deviation (no units)

## Input models
- `idealmodel`: Ideal Model

## Description

SAFT-BACK equation of state. Also known as SAFT-CP.

## References
1. Chen, J., & Mi, J.-G. (2001). Equation of state extended from SAFT with improved results for non-polar fluids across the critical point. Fluid Phase Equilibria, 186(1–2), 165–184. [doi:10.1016/s0378-3812(01)00521-0](https://doi.org/10.1016/s0378-3812(01)00521-0)
2. Mi, J.-G., Chen, J., Gao, G.-H., & Fei, W.-Y. (2002). Equation of state extended from SAFT with improved results for polar fluids across the critical point. Fluid Phase Equilibria, 201(2), 295–307. [doi:10.1016/s0378-3812(02)00093-6](https://doi.org/10.1016/s0378-3812(02)00093-6)
3. Maghari, A., & Hamzehloo, M. (2011). Second-order thermodynamic derivative properties of binary mixtures of n-alkanes through the SAFT-CP equation of state. Fluid Phase Equilibria, 302(1–2), 195–201. [doi:10.1016/j.fluid.2010.08.019](https://doi.org/10.1016/j.fluid.2010.08.019)
"""
BACKSAFT

recombine_impl!(model::BACKSAFTModel) = recombine_saft!(model)


function lb_volume2(model::BACKSAFTModel,T,z)
    V = 0.0
    m = model.params.segment.values
    dᵢ = @f(d)
    #dᵢ = model.params.sigma.values
    val = 2π/6*N_A*sum(z[i]*m[i]*dᵢ[i,i]^3 for i in 1:length(z))
    return val
end

function lb_volume(model::BACKSAFTModel,T,z)
    α = dot(model.params.alpha.values,z)/sum(z)

    #TODO: remember how i got to this
    poly = (1.0,3α-2,3α*α - 3α +1 , -α*α)
    dpoly = Solvers.polyder(poly)
    function pol(x)
        f = evalpoly(x,poly)
        df = evalpoly(x,dpoly)
        return f,f/df
    end
    prob = Roots.ZeroProblem(pol,1.81)
    k = Roots.solve(prob,Roots.Newton())
    
    m = model.params.segment.values
    σᵢᵢ = model.params.sigma.values
    V = 0.0
    σᵢᵢ = @f(d)
    val = π/6*N_A*sum(z[i]*m[i]*σᵢᵢ[i,i]^3 for i in 1:length(z)) #limit at η -> 1
    return k*val
end

function x0_volume_liquid(model::BACKSAFTModel,p,T,z)
    lbv = lb_volume(model,T,z)
    V0 = 1.01*lbv
    for i in 1:10
        pressure(model,V0,T,z) > p && break
        V0 = 0.5*V0 + 0.5*lbv
    end
    return V0
end

function x0_crit_pure(model::BACKSAFTModel,z)
    T = T_scale(model,z)
    lb_v = lb_volume(model,T,z)/sum(z)
    (2.0, log10(3lb_v))
end

function data(model::BACKSAFTModel,V, T, z)
    _d = @f(d)
    n = sum(z)
    m = model.params.segment.values
    η = ζ(model,V,T,z,3,_d,m)
    m̄ = dot(m,z)/n
    return _d, n, η, m̄
end

function a_res(model::BACKSAFTModel ,V, T, z,_data = @f(data))
    a_hcb_ = @f(a_hcb,_data)
    a_disp_ = @f(a_disp,_data)
    a_chain_ = @f(a_chain,_data)
    #@show a_hcb_
    #@show a_chain_
    #@show a_disp_
    #@show (1.75*(a_chain_/a_hcb_))*a_disp_
    return  a_hcb_ + a_chain_ + (1.75*(a_chain_/a_hcb_) + 1)*a_disp_
end

function a_hcb(model::BACKSAFTModel, V, T, z,_data = @f(data))
    _d, n, η, m̄ = _data
    αi = model.params.alpha.values
    α = dot(αi,z)/n
    α2 = α*α
    x = 1/(1 - η)
    #@show α2*x*x
    #@show - (α2-3α)*x
    #@show - (1-α2)*log1p(-η)
    #@show - 3α
    return m̄*(α2*x*x - (α2-3α)*x - (1-α2)*log1p(-η) - 3α)
end

function a_disp(model::BACKSAFTModel, V, T, z,_data = @f(data))
    _d, n, η, m̄ = _data
    c = model.params.c.values
    u0 = model.params.epsilon.values
    m = model.params.segment.values
    τ = 0.7404804896930611 #sqrt(2)/6 * pi
    D1 = BACKSAFT_consts.D1
    D2 = BACKSAFT_consts.D2
    D3 = BACKSAFT_consts.D3
    D4 = BACKSAFT_consts.D4
    nc = length(model)
    ∑mxx = zero(Base.promote_eltype(model,V,T,z))
    ∑muxx = zero(∑mxx)
    for i in 1:nc
        zi,mi,u0i,ci = z[i],m[i],u0[i],c[i]
        mxx = mi*mi*zi*zi
        ∑mxx += mxx
        uii = u0i*(1+ci/T)
        ∑muxx += mxx*uii
        for j in 1:i-1
            zj,mj,u0j,cj = z[j],m[j],u0[j],c[j]
            mxy = mi*mj*zi*zj
            ∑mxx += 2*mxy
            ujj = u0j*(1+cj/T)
            uij = sqrt(uii*uij)
            ∑muxx += 2*mxy*uij
        end
    end
    u = ∑muxx/∑mxx
    u1 = u/T
    u2 = u1*u1
    u3 = u2*u1
    u4 = u2*u2
    η̄ = (η/τ)
    A1 = u1*evalpoly(η̄,D1)*η̄
    A2 = u2*evalpoly(η̄,D2)*η̄
    A3 = u3*evalpoly(η̄,D3)*η̄
    A4 = u4*evalpoly(η̄,D4)*η̄
    return m̄*(A1+A2+A3+A4)
end

d(model::BACKSAFTModel, V, T, z) = ck_diameter(model, T, z)

function a_chain(model::BACKSAFTModel, V, T, z,_data = @f(data))
    _d, n, η, m̄ = _data
    m = model.params.segment.values
    α = model.params.alpha.values
    res = zero(Base.promote_eltype(model,V,T,z))
    for i in 1:length(z)
        res += z[i]*(1 - m[i])*log(@f(g_hcb,α[i],_data))
    end
    return res/n
end
#1/(1-ζ3) + di*dj/(di+dj)*3ζ2/(1-ζ3)^2 + (di*dj/(di+dj))^2*2ζ2^2/(1-ζ3)^3
function g_hcb(model::BACKSAFTModel, V, T, z, α, _data = @f(data))
    _d, m̄, η, _ = _data
    x = 1/(1 - η)
    return x + x*x*3*(1+α)*α*η/(1+3α) + x^3*2*η^2*α^2/(1+3α)
end

const BACKSAFT_consts = (
    D1 = [-8.8043,4.164627,-48.203555,140.4362,-195.23339,113.515],
    D2 = [2.9396,-6.0865383,40.137956,-76.230797,-133.70055,860.25349,-1535.3224,1221.4261,-409.10539],
    D3 = [-2.8225,4.7600148,11.257177,-66.382743,69.248785],
    D4 = [0.34,-3.1875014,12.231796,-12.110681],
)