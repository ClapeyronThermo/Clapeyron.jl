struct RKParam <: EoSParam
    a::PairParam{Float64}
    b::PairParam{Float64}
    Tbarc::Float64 # Not sure if we want to allow this
end

abstract type RKModel <: CubicModel end
@newmodel RK RKModel RKParam

export RK
function RK(components::Array{String,1}; userlocations::Array{String,1}=String[], verbose=false)
    params = getparams(components, ["properties/critical.csv", "SAFT/PCSAFT/PCSAFT_unlike.csv"]; userlocations=userlocations, verbose=verbose)
    k  = params["k"]
    pc = params["pc"].values
    Tc = params["Tc"].values
    T̄c = sum(sqrt(Tc*Tc'))

    a = epsilon_LorentzBerthelot(SingleParam(params["pc"], @. 1/(9*(2^(1/3)-1))*R̄^2*Tc^2.5/pc/1e6/√(T̄c)), k)
    b = sigma_LorentzBerthelot(SingleParam(params["pc"], @. (2^(1/3)-1)/3*R̄*Tc/pc/1e6))

    packagedparams = RKParam(a, b, T̄c)
    return RK(packagedparams)
end

function a_tot(model::RKModel, V, T, z)
    x = z/sum(z)
    n = sum(z)
    a = model.params.a.values
    b = model.params.b.values
    T̄c = model.params.Tbarc
    ā = sum(a .* (x * x'))
    b̄ = sum(b .* (x * x'))
    return -log(V-n*b̄) - ā/(R̄*T*b̄*√(T/T̄c))*log(1+n*b̄/V)
end
