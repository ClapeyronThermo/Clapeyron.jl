abstract type DIPPR101SatModel <: SaturationModel end

struct DIPPR101SatParam <: EoSParam 
    Tc::SingleParam{Float64}
    Pc::SingleParam{Float64}
    A::SingleParam{Float64}
    B::SingleParam{Float64}
    C::SingleParam{Float64}
    D::SingleParam{Float64}
    E::SingleParam{Float64}
    Tmin::SingleParam{Float64}
    Tmax::SingleParam{Float64}
end

@newmodelsimple DIPPR101Sat DIPPR101SatModel DIPPR101SatParam

function DIPPR101Sat(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false)
    params = getparams(components, ["properties/critical.csv","Correlations/saturation_correlations/dippr101_like.csv"]; userlocations=userlocations, verbose=verbose)
    Tc = params["Tc"]
    Pc = params["pc"]
    A = params["A"]
    B = params["B"]
    C = params["C"]
    D = params["D"]
    E = params["E"]
    Tmin = params["Tmin"]
    Tmax = params["Tmax"]
    packagedparams = DIPPR101SatParam(Tc,Pc,A,B,C,D,E,Tmin,Tmax)
    model = DIPPR101Sat(packagedparams, verbose=verbose)
    return model
end 

function crit_pure(model::DIPPR101SatModel)
    tc = only(model.params.Tc.values)
    pc = only(model.params.Pc.values)
    return (tc,pc,NaN)
end

function saturation_pressure_impl(model::DIPPR101SatModel,T,method::SaturationCorrelation)
    nan = zero(T)/zero(T)
    tc = only(model.params.Tc.values)
    A = only(model.params.A.values)
    B = only(model.params.B.values)
    C = only(model.params.C.values)
    D = only(model.params.D.values)
    E = only(model.params.E.values)
    Tmin = only(model.params.Tmin.values)
    Tmax = only(model.params.Tmax.values)

    T > tc && (return nan,nan,nan)
    Tmin <= T <= Tmax || (return nan,nan,nan)
    psat = exp(A + B/T + C*log(T) + D*T^E)
    return psat,nan,nan
end

export DIPPR101Sat
