abstract type DIPPR105LiquidModel <: LiquidVolumeModel end


struct DIPPR105LiquidParam <: EoSParam
    A::SingleParam{Float64}
    B::SingleParam{Float64}
    C::SingleParam{Float64}
    D::SingleParam{Float64}
    Tmin::SingleParam{Float64}
    Tmax::SingleParam{Float64}
end

@newmodelsimple DIPPR105Liquid DIPPR105LiquidModel DIPPR105LiquidParam

function DIPPR105Liquid(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false)
    params = getparams(components, ["Correlations/volume_correlations/dippr105_like.csv"]; userlocations=userlocations, verbose=verbose)
    A = params["A"]
    B = params["B"]
    C = params["C"]
    D = params["D"]
    Tmin = params["Tmin"]
    Tmax = params["Tmax"]
    packagedparams = DIPPR105LiquidParam(A,B,C,D,Tmin,Tmax)
    model = DIPPR105Liquid(packagedparams)
    return model
end

function volume_impl(model::DIPPR105LiquidModel,p,T,z=SA[1.0],phase=:unknown,threaded=false,vol0 = nothing)
    A,B,C,D = model.params.A.values,model.params.B.values,model.params.C.values,model.params.D.values
    #Tmin,Tmax = model.params.Tmin.values,model.params.Tmax.values
    res = zero(T + first(z))
    for i in 1:length(z)
        bexp = 1 + (1 - T/C[i])^D[i]
        vi = (B[i]^bexp)/A[i]
        res += z[i]*vi
    end
    return res
end

export DIPPR105Liquid