struct PRParam <: EoSParam
    a::PairParam{Float64}
    b::PairParam{Float64}
    acentricfactor::SingleParam{Float64}
    Tc::SingleParam{Float64}
end

abstract type PRModel <: CubicModel end
@newmodel PR PRModel PRParam

export PR
function PR(components::Array{String,1}; userlocations::Array{String,1}=String[], verbose=false)
    params = getparams(components, ["properties/critical.csv", "SAFT/PCSAFT/PCSAFT_unlike.csv"]; userlocations=userlocations, verbose=verbose)
    k  = params["k"]
    pc = params["pc"].values
    Tc = params["Tc"]
    Tc_ = Tc.values
    acentricfactor = params["w"]
    a = epsilon_LorentzBerthelot(SingleParam(params["pc"], @. 0.457235*R̄^2*Tc_^2/pc/1e6), k)
    b = sigma_LorentzBerthelot(SingleParam(params["pc"], @. 0.077796*R̄*Tc_/pc/1e6))

    packagedparams = PRParam(a, b, acentricfactor, Tc)
    return PR(packagedparams)
end

function a_tot(model::PRModel, V, T, z)
    x = z/sum(z)
    n = sum(z)
    a = model.params.a.values
    b = model.params.b.values
    ω = model.params.acentricfactor.values
    Tc = model.params.Tc.values

    α = @. (1+(0.37464+1.54226*ω-0.26992*ω^2)*(1-√(T/Tc)))^2

    āᾱ = sum(a .* .√(α * α') .* (x * x'))
    b̄ = sum(b .* (x * x'))
    return -log(V-n*b̄) + āᾱ/(R̄*T*b̄*2^(3/2)) * log((2*V-2^(3/2)*b̄*n+2*b̄*n)/(2*V+2^(3/2)*b̄*n+2*b̄*n))
end


#=
function cubic_α(model::PRFamily,t,i,j)

    m_poly = (0.37464,1.54226,0.26992)
    _1 = one(t)
    aᵢ =model._a[i]
    tcᵢ = model.tc[i]
    sqrt_trᵢ = min(sqrt(t/tcᵢ),_1) #α(t) is one for supercritical values
    ωᵢ = model.ω[i]
    mᵢ = evalpoly(ωᵢ,m_poly)
    if i === j
        return  aᵢ* (_1+mᵢ * √(_1-sqrt_trᵢ))^2
    else
        aⱼ =model._a[j]
        tcⱼ = model.tc[j]
        sqrt_trⱼ = min(sqrt(t/tcⱼ),_1)
        ωⱼ = model.ω[j]
        mⱼ = evalpoly(ωⱼ,m_poly)
        sqrt_αᵢ= (_1+mᵢ * √(_1-sqrt_trᵢ))
        sqrt_αⱼ= (_1+mⱼ * √(_1-sqrt_trⱼ))
        return sqrt(aᵢ*aⱼ)*sqrt_αᵢ*sqrt_αⱼ
    end
end



function cubic_ab(model::PRFamily,p,t,x)
    #two options to introduce alpha:
    #here: it will allocate, but less ops
    bi = model._b
    b = dot(bi,x)
    #here: it will not allocate, but more ops
    sss= (i,j)->cubic_aα(model,t,i,j)
    a = cubic_mixing_rule(sss, x, model.aij)
    return a,b
end

function cubic_abp(mt::SingleVT,model::PengRobinson{SINGLE},v,t)
    a,b = cubic_ab(QuickStates.pt(),model,v,t) #v is ignored
    _1 = one(b)
    denom = evalpoly(v,(-b*b,2*b,_1))
    p = RGAS*t/(v-b) - a/denom
    return a,b,p
end

function cubic_abp(mt::MultiVT,model::PengRobinson{MULTI},v,t,x)
    a,b = cubic_ab(QuickStates.ptx(),model,v,t,x) #v is ignored
    _1 = one(b)
    denom = evalpoly(v,(-b*b,2*b,_1))
    p = RGAS*t/(v-b) - a/denom
    return a,b,p
end

function fugacity_coeff_impl(mt::SingleVT,model::PengRobinson{SINGLE},v,t)
    a,b,p =  cubic_abp(mt,model,v,t)
    RTinv = 1/(RGAS*t)
    A = a*p*RTinv*RTinv
    B = b*p*RTinv
    z = p*v*RTinv
     _1 = one(z)
     logϕ = z - _1 - log(z-B) - A/z
end

const PRΔ1 = 1+√2
const PRΔ2 = 1-√2
const ΔPRΔ = 2*√2

function  αR_impl(mt::MultiVT,model::PengRobinson{MULTI},rho,t,x)
    R = RGAS
    RTinv = 1/(RGAS*t)
    v = inv(rho)
    a,b,p =  cubic_abp(mt,model,v,t,x)
    -log(1-b*rho) - a*RTinv*log((PRΔ1*b*rho+1)/(PRΔ2*b*rho+1))/(ΔPRΔ*b)
end


function cubic_poly(mt::SinglePT,model::PengRobinson{SINGLE},p,t)
    a,b = cubic_ab(QuickStates.pt(),model,p,t)
    RTinv = 1/(RGAS*t)
    A = a*p*RTinv*RTinv
    B = b*p*RTinv
    k0 = B*(B*(B+1.0)-A)
    k1 = -B*(2*B+1.0) + A
    k2 = -1.0
    k3 = 1.0
    return (k0,k1,k2,k3)
end

function cubic_poly(mt::MultiPT,model::PengRobinson{MULTI},p,t,x)
    a,b = cubic_ab(QuickStates.ptx(),model,p,t,x)
    RTinv = 1/(RGAS*t)
    A = a*p*RTinv*RTinv
    B = b*p*RTinv
    k0 = B*(B*(B+1.0)-A)
    k1 = -B*(2*B+1.0) + A
    k2 = -1.0
    k3 = 1.0
    return (k0,k1,k2,k3)
end

=#
