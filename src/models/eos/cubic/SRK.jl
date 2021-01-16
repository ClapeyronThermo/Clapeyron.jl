struct SRKParam <: EoSParam
    a::PairParam{Float64}
    b::PairParam{Float64}
    acentricfactor::SingleParam{Float64}
    Tc::SingleParam{Float64}
end

abstract type SRKModel <: CubicModel end
@newmodel SRK SRKModel SRKParam

export SRK
function SRK(components::Array{String,1}; userlocations::Array{String,1}=String[], verbose=false)
    params = getparams(components, ["properties/critical.csv", "SAFT/PCSAFT/PCSAFT_unlike.csv"]; userlocations=userlocations, verbose=verbose)
    k  = params["k"]
    pc = params["pc"].values
    Tc = params["Tc"]
    Tc_ = Tc.values
    acentricfactor = params["w"]
    a = epsilon_LorentzBerthelot(SingleParam(params["pc"], @. 1/(9*(2^(1/3)-1))*R̄^2*Tc_^2/pc/1e6), k)
    b = sigma_LorentzBerthelot(SingleParam(params["pc"], @. (2^(1/3)-1)/3*R̄*Tc_/pc/1e6))

    packagedparams = SRKParam(a, b, acentricfactor, Tc)
    return SRK(packagedparams)
end

function a_tot(model::SRKModel, V, T, z)
    x = z/sum(z)
    n = sum(z)
    a = model.params.a.values
    b = model.params.b.values
    ω = model.params.acentricfactor.values
    Tc = model.params.Tc.values

    α = @. (1+(0.480+1.547*ω-0.176*ω^2)*(1-√(T/Tc)))^2

    āᾱ = sum(a .* .√(α * α') .* (x * x'))
    b̄ = sum(b .* (x * x'))
    return -log(V-n*b̄) - āᾱ/(R̄*T*b̄)*log(1+n*b̄/V)
end
