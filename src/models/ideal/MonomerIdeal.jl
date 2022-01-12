struct MonomerIdealParam <: EoSParam
    Mw::SingleParam{Float64}
end

abstract type MonomerIdealModel <: IdealModel end
@newmodelsimple MonomerIdeal MonomerIdealModel MonomerIdealParam

export MonomerIdeal
function MonomerIdeal(components::Array{String,1}; userlocations::Array{String,1}=String[], verbose=false)
    params = getparams(components, ["properties/molarmass"]; userlocations=userlocations, verbose=verbose)
    Mw = params["Mw"]
    packagedparams = MonomerIdealParam(Mw)
    return MonomerIdeal(packagedparams)
end

function a_ideal(model::MonomerIdealModel, v, T, z)
    x = z/sum(z)
    Mw = model.params.Mw.values
    Λ = @. h/√(k_B*T*Mw/N_A)
    return sum(@. x*log(z*N_A/v*Λ^3))-1.
end
