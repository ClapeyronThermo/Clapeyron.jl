struct vdWParam <: EoSParam
    a::PairParam{Float64}
    b::PairParam{Float64}
end

abstract type vdWModel <: CubicModel end
@newmodel vdW vdWModel vdWParam

export vdW
function vdW(components::Array{String,1}; userlocations::Array{String,1}=String[], verbose=false)
    params = getparams(components, ["properties/critical.csv", "SAFT/PCSAFT/PCSAFT_unlike.csv"]; userlocations=userlocations, verbose=verbose)

    k = params["k"]
    pc = params["pc"].values
    Tc = params["Tc"].values

    a = epsilon_LorentzBerthelot(SingleParam(params["pc"], @. 27/64*R̄^2*Tc^2/pc/1e6), k)
    b = sigma_LorentzBerthelot(SingleParam(params["pc"], @. 1/8*R̄*Tc/pc/1e6))

    packagedparams = vdWParam(a, b)
    return vdW(packagedparams)
end

function a_tot(model::vdWModel, V, T, z)
    x = z/sum(z)
    n = sum(z)
    a = model.params.a.values
    b = model.params.b.values
    ā = sum(a .* (x * x'))
    b̄ = sum(b .* (x * x'))
    return -log(V-n*b̄) - ā*n/(R̄*T*V)
end

function a_res(model::vdWModel, V, T, z)
    return @f(a_tot) + log(V)  # + f(x)
end
