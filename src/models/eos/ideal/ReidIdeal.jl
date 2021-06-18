struct ReidIdealParam <: EoSParam
    a::SingleParam{Float64}
    b::SingleParam{Float64}
    c::SingleParam{Float64}
    d::SingleParam{Float64}
end

abstract type ReidIdealModel <: IdealModel end
@newmodelsimple ReidIdeal ReidIdealModel ReidIdealParam

export ReidIdeal
function ReidIdeal(components::Array{String,1}; userlocations::Array{String,1}=String[], verbose=false)
    params = getparams(components, ["ideal/ReidIdeal.csv"]; userlocations=userlocations, verbose=verbose)
    a = params["a"]
    b = params["b"]
    c = params["c"]
    d = params["d"]
    packagedparams = ReidIdealParam(a, b, c, d)
    references = String[] #  Fill this up.
    return ReidIdeal(packagedparams; references=references)
end

function a_ideal(model::ReidIdealModel, V, T, z)
    x = z/sum(z)
    a = model.params.a.values
    b = model.params.b.values
    c = model.params.c.values
    d = model.params.d.values
    polycoeff = [a, b, c, d]
    return sum(x[i]*(log(z[i]/V) + 1/(R̄*T)*(sum(polycoeff[k][i]/k*(T^k-298^k) for k in 1:4)) -
        1/R̄*((a[i]-R̄)*log(T/298)+sum(polycoeff[k][i]/(k-1)*(T^(k-1)-298^(k-1)) for k in 2:4))) for i in @comps)
end
