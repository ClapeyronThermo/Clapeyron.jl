struct RKParam <: EoSParam
    a::PairParam{Float64}
    b::PairParam{Float64}
    Tc::SingleParam{Float64}
    Pc::SingleParam{Float64}
    Mw::SingleParam{Float64}
    Tbarc::Float64 # Not sure if we want to allow this
end

abstract type RKModel <: ABCubicModel end
@newmodel RK RKModel RKParam

export RK
function RK(components::Array{String,1}; userlocations::Array{String,1}=String[], verbose=false,idealmodel=BasicIdeal)
    params = getparams(components, ["properties/critical.csv", "properties/molarmass.csv","SAFT/PCSAFT/PCSAFT_unlike.csv"]; userlocations=userlocations, verbose=verbose)
    k  = params["k"]
    _pc = params["pc"]
    pc = _pc.values
    Mw = params["Mw"]
    _Tc = params["Tc"]
    Tc = _Tc.values
    T̄c = sum(sqrt(Tc*Tc'))
    a = epsilon_LorentzBerthelot(SingleParam(params["pc"], @. 1/(9*(2^(1/3)-1))*R̄^2*Tc^2.5/pc/√(T̄c)), k)
    b = sigma_LorentzBerthelot(SingleParam(params["pc"], @. (2^(1/3)-1)/3*R̄*Tc/pc))

    packagedparams = RKParam(a, b, params["Tc"],_pc,Mw,T̄c)
    model = RK(packagedparams,idealmodel)
    return model
end
function ab_consts(::Type{<:RKModel})
    Ωa =  1/(9*(2^(1/3)-1))
    Ωb = (2^(1/3)-1)/3
    return Ωa,Ωb
end
function cubic_ab(model::RKModel,T,z=SA[1.0],n=sum(z))
    invn2 = (one(n)/n)^2
    a = model.params.a.values
    b = model.params.b.values
    Tc = model.params.Tc.values
    T̄c = model.params.Tbarc
    āᾱ  = dot(z,Symmetric(a),z) * invn2*sqrt(T̄c/T)
    b̄ = dot(z,Symmetric(b),z) * invn2
    return āᾱ ,b̄
end

function cubic_abp(model::RKModel, V, T, z)
    n = sum(z)
    a,b = cubic_ab(model,T,z,n)
    v = V/n
    p =  R̄*T/(v-b) - a/((v+b)*v)
    return a,b,p
end

function cubic_poly(model::RKModel,p,T,z)
    n = sum(z)
    a,b = cubic_ab(model,T,z,n)
    RT⁻¹ = 1/(R̄*T)
    A = a*p* RT⁻¹* RT⁻¹
    B = b*p* RT⁻¹
    _1 = one(a)
    return [-A*B, -B*(B+_1) + A, -_1, _1]
end


function a_res(model::RKModel, V, T, z)
    n = sum(z)
    ā,b̄ = cubic_ab(model,T,z,n)
    ρ = n/V
    RT⁻¹ = 1/(R̄*T)
    return -log(1-b̄*ρ) - ā*RT⁻¹*log(b̄*ρ+1)/b̄
    #return -log(V-n*b̄) - ā/(R̄*T*b̄*√(T/T̄c))*log(1+n*b̄/V)
end


cubic_zc(::RKModel) = 1/3
