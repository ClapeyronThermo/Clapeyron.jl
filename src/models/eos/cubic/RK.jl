struct RKParam <: EoSParam
    a::PairParam{Float64}
    b::PairParam{Float64}
    Tc::SingleParam{Float64}
    Tbarc::Float64 # Not sure if we want to allow this
end

abstract type RKModel <: ABCubicModel end
@newmodel RK RKModel RKParam

export RK
function RK(components::Array{String,1}; userlocations::Array{String,1}=String[], verbose=false,idealmodel=BasicIdeal)
    params = getparams(components, ["properties/critical.csv", "SAFT/PCSAFT/PCSAFT_unlike.csv"]; userlocations=userlocations, verbose=verbose)
    k  = params["k"]
    pc = params["pc"].values
    Tc = params["Tc"].values
    T̄c = sum(sqrt(Tc*Tc'))
    idealmodel = idealmodelselector(idealmodel, components,verbose=verbose)
    a = epsilon_LorentzBerthelot(SingleParam(params["pc"], @. 1/(9*(2^(1/3)-1))*R̄^2*Tc^2.5/pc/1e6/√(T̄c)), k)
    b = sigma_LorentzBerthelot(SingleParam(params["pc"], @. (2^(1/3)-1)/3*R̄*Tc/pc/1e6))

    packagedparams = RKParam(a, b, params["Tc"],T̄c)
    return RK(packagedparams,idealmodel)
end

function cubic_ab(model::RKModel,T,x)
    a = model.params.a.values
    b = model.params.b.values
    Tc = model.params.Tc.values
    T̄c = model.params.Tbarc
    āᾱ  = sum(a .* (x * x'))*sqrt(T̄c/T)
    b̄ = sum(b .* (x * x'))
    return āᾱ ,b̄
end

function cubic_abp(model::RKModel, V, T, x)
    a,b = cubic_ab(model,T,x)
    p =  R̄*T/(V-b) - a/((V+b)*V)
    return a,b,p
end

function cubic_poly(model::RKModel,p,t,x)
    a,b = cubic_ab(model,T,x)
    RT⁻¹ = 1/(R̄*T)
    A = a*p* RT⁻¹* RT⁻¹
    B = b*p* RT⁻¹
    _1 = one(a)
    return (-A*B, -B*(B+_1) + A, -_1, _1)
end

function a_resx(model::RKModel, v, T, x)
    ā,b̄ = cubic_ab(model,T,x)
    ρ = 1/v
    RT⁻¹ = 1/(R̄*T)
    return -log(1-b̄*ρ) - ā*RT⁻¹*log(b̄*ρ+1)/b̄
    #return -log(V-n*b̄) - ā/(R̄*T*b̄*√(T/T̄c))*log(1+n*b̄/V)
end
