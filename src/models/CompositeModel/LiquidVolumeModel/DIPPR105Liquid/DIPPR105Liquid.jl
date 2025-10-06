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
default_locations(::Type{DIPPR105Liquid}) = ["Correlations/volume_correlations/dippr105_like.csv"]

function volume_impl(model::DIPPR105LiquidModel,p,T,z,phase,threaded,vol0)
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