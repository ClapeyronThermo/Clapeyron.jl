abstract type LeeKeslerSatModel <: SaturationModel end

struct LeeKeslerSatParam <: EoSParam 
    Tc::SingleParam{Float64}
    Pc::SingleParam{Float64}
    acentricfactor::SingleParam{Float64}
end

@newmodelsimple LeeKeslerSat LeeKeslerSatModel LeeKeslerSatParam

function LeeKeslerSat(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false)
    params = getparams(components, ["properties/critical.csv"]; userlocations=userlocations, verbose=verbose)
    acentricfactor = SingleParam(params["w"],"acentric factor")
    Tc = params["Tc"]
    Pc = params["pc"]
    packagedparams = LeeKeslerSatParam(Tc,Pc,acentricfactor)
    model = LeeKeslerSat(packagedparams, verbose=verbose)
    return model
end 

function crit_pure(model::LeeKeslerSatModel)
    tc = only(model.params.Tc.values)
    pc = only(model.params.Pc.values)
    return (tc,pc,NaN)
end

function saturation_pressure(model::LeeKeslerSatModel,T,v0=nothing)
    nan = zero(T)/zero(T)
    ω = only(model.params.acentricfactor.values)
    tc = only(model.params.Tc.values)
    pc = only(model.params.Pc.values)
    T > tc && (return nan,nan,nan)
    
    tr = T/tc
    trinv = inv(tr)
    lntr = log(tr)
    tr6 = tr^6
    
    f0 = 5.92714 - 6.09648*trinv - 1.28862*lntr + 0.169347*tr6
    f1 = 15.2518 - 15.6875*trinv - 13.4721*lntr + 0.43577*tr6
    lnpr = f0 + ω*f1
    psat = exp(lnpr)*pc
    return psat,nan,nan
end

export LeeKeslerSat
