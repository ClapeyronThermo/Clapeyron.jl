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
default_references(::Type{BACKSAFT}) = ["10.1016/s0378-3812(02)00093-6"]
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
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `vol`: Single Parameter (`Float64`) - Segment Volume [`dm^3`]
- `epsilon`: Single Parameter (`Float64`) - Reduced dispersion energy  `[K/mol]`
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Paramater (no units)
- `c`: Single Parameter (`Float64`) - Adjustable parameter (no units)
- `alpha`: Single Parameter (`Float64`) - Non-spherical deviation (no units)

## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Pair Parameter (`Float64`) - Mixed segment Diameter `[m]`
- `epsilon`: Pair Parameter (`Float64`) - Mixed reduced dispersion energy`[K]`
- `c`: Single Parameter (`Float64`) - Adjustable parameter (no units)
- `alpha`: Single Parameter (`Float64`) - Non-spherical deviation (no units)

## Input models
- `idealmodel`: Ideal Model

## Description

BACKSAFT

## References
1. Mi, J.-G., Chen, J., Gao, G.-H., & Fei, W.-Y. (2002). Equation of state extended from SAFT with improved results for polar fluids across the critical point. Fluid Phase Equilibria, 201(2), 295–307. [doi:10.1016/s0378-3812(02)00093-6](https://doi.org/10.1016/s0378-3812(02)00093-6)
"""
BACKSAFT

recombine_impl!(model::BACKSAFTModel) = recombine_saft!(model)


function lb_volume(model::BACKSAFTModel,z)
    α = model.params.alpha.values[1]
    pol(x) = evalpoly(x,(1.0,3α-2,3α*α - 3α +1 , -α*α))
    k = Solvers.ad_newton(pol,1.81)
    seg = model.params.segment.values
    σᵢᵢ = model.params.sigma.values
    val = π/6*N_A*sum(z[i]*seg[i]*σᵢᵢ[i,i]^3 for i in 1:length(z)) #limit at η -> 0
    return k*val #only positive root of η 
end

function x0_volume_liquid(model::BACKSAFTModel,T,z)
    v_lb = lb_volume(model,z)
    return v_lb*1.01
end

function x0_crit_pure(model::BACKSAFTModel)
    lb_v = lb_volume(model)
    (2.0, log10(lb_v/0.4))
end

function data(model::BACKSAFTModel,V, T, z)
    _d = @f(d)
    m̄ = dot(z,model.params.segment.values)
    ζi = @f(ζ0123,_d)
    return _d, m̄, ζi
end

function a_res(model::BACKSAFTModel ,V, T, z,_data = @f(data))
    a_hcb_ = @f(a_hcb,_data)
    a_disp_ = @f(a_disp,_data)
    a_chain_ = @f(a_chain,_data)
    return  a_hcb_ + a_chain_ + (1.75*(a_chain_/a_hcb_)+1)*a_disp_
end

function ζ0123(model::BACKSAFTModel, V, T, z,_d=@f(d))
    #N_A*π/6/V * sum(z[i]*m[i]*@f(d,i)^n for i ∈ @comps)
    m = model.params.segment
    _0 = zero(V+T+first(z)+one(eltype(model)))
    ζ0,ζ1,ζ2,ζ3 = _0,_0,_0,_0
    for i ∈ 1:length(z)
        di =_d[i]
        xS = z[i]*m[i]
        ζ0 += xS
        ζ1 += xS*di
        ζ2 += xS*di*di
        ζ3 += xS*di*di*di
    end
    c = π/6*N_A/V
    ζ0,ζ1,ζ2,ζ3 = c*ζ0,c*ζ1,c*ζ2,c*ζ3
    return ζ0,ζ1,ζ2,ζ3
end

#=

z = 1/(1-y) + 3ay/(1-y)2 + (3a2y2-a2y3)/(1-y)^3
3a2y2(1-y)
=#
function a_hcb(model::BACKSAFTModel, V, T, z,_data = @f(data))
    _d, m̄, ζi = _data
    _,_,_,η = ζi
    α = model.params.alpha.values[1]
    m = model.params.segment.values[1]
    return m*(α^2/(1-η)^2-(α^2-3α)/(1-η)-(1-α^2)*log(1-η)-3α)
end

function a_disp(model::BACKSAFTModel, V, T, z,_data = @f(data))
    _d, m̄, ζi = _data
    _,_,_,η = ζi
    m = model.params.segment.values[1]
    c = model.params.c.values[1]
    u = model.params.epsilon.values[1]*(1+c/T)
    τ = .740480
    D1 = BACKSAFT_consts.D1
    D2 = BACKSAFT_consts.D2
    D3 = BACKSAFT_consts.D3
    D4 = BACKSAFT_consts.D4
    A1 = sum(D1[j]*(u/T)*(η/τ)^j for j ∈ 1:6)
    A2 = sum(D2[j]*(u/T)^2*(η/τ)^j for j ∈ 1:9)
    A3 = sum(D3[j]*(u/T)^3*(η/τ)^j for j ∈ 1:5)
    A4 = sum(D4[j]*(u/T)^4*(η/τ)^j for j ∈ 1:4)
    return m*(A1+A2+A3+A4)
end

d(model::BACKSAFTModel, V, T, z) = ck_diameter(model, T, z)

function ζ(model::BACKSAFTModel, V, T, z, n)
    m = model.params.segment.values
    _d = @f(d)
    return N_A*π/6/V * sum(z[i]*m[i]*_d[i]^n for i ∈ @comps)
end

function a_chain(model::BACKSAFTModel, V, T, z,_data = @f(data))
    _d, m̄, ζi = _data
    _,_,_,η = ζi
    m = model.params.segment.values[1]
    return (1-m)*log(@f(g_hcb,_data))
end
#1/(1-ζ3) + di*dj/(di+dj)*3ζ2/(1-ζ3)^2 + (di*dj/(di+dj))^2*2ζ2^2/(1-ζ3)^3
function g_hcb(model::BACKSAFTModel, V, T, z,_data = @f(data))
    _d, m̄, ζi = _data
    _,_,_,η = ζi
    α = model.params.alpha.values[1]
    return 1/(1-η)+3*(1+α)*α*η/((1-η)^2*(1+3α))+2*η^2*α^2/((1-η)^3*(1+3α))
end

const BACKSAFT_consts = (
    D1 = [-8.8043,4.164627,-48.203555,140.4362,-195.23339,113.515],
    D2 = [2.9396,-6.0865383,40.137956,-76.230797,-133.70055,860.25349,-1535.3224,1221.4261,-409.10539],
    D3 = [-2.8225,4.7600148,11.257177,-66.382743,69.248785],
    D4 = [0.34,-3.1875014,12.231796,-12.110681],
)
