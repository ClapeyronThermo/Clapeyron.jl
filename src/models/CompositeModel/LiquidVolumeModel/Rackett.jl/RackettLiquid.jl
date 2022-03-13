abstract type RackettLiquidModel <: LiquidVolumeModel end


struct RackettLiquidParam <: LiquidVolumeModel
    Tc::SingleParam{Float64}
    Pc::SingleParam{Float64}
    Zc::SingleParam{Float64}
end

@newmodelsimple RackettLiquid RackettLiquidModel RackettLiquidParam


function RackettLiquid(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false)
    params = getparams(components, ["properties/critical.csv"]; userlocations=userlocations, verbose=verbose)
    acentricfactor = SingleParam(params["w"],"acentric factor")
    Tc = params["Tc"]
    Pc = params["pc"]
    if haskey(params,"Zc") || haskey(params,"zc")
        Zc = params["Zc"]
    else
        vc = params["vc"]
        _zc = Pc.values .* vc.values ./ R̄ .* Tc.values
        Zc = SingleParam("Critical Compressibility factor",components,_zc)
    end
    packagedparams = RackettLiquidParam(Tc,Pc,Zc)
    model = RackettLiquid(packagedparams, verbose=verbose)
    return model
end

export RackettLiquid 

function volume(model::RackettLiquidModel,p,T,z=SA[1.0];phase=:unknown,threaded=false)
    for i ∈ @comps
    
    end
end