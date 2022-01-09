abstract type BMAlphaModel <: AlphaModel end

struct BMAlphaParam <: EoSParam
    acentricfactor::SingleParam{Float64}
end

@newmodelsimple BMAlpha BMAlphaModel BMAlphaParam

export BMAlpha
function BMAlpha(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false)
    params = getparams(components, ["properties/critical.csv"]; userlocations=userlocations, verbose=verbose)
    acentricfactor = SingleParam(params["w"],"acentric factor")
    packagedparams = BMAlphaParam(acentricfactor)
    model = BMAlpha(packagedparams, verbose=verbose)
    return model
end

function α_function(model::RKModel,V,T,z,alpha_model::BMAlphaModel)
    Tc = model.params.Tc.values
    Tr = @. T/Tc
    ω  = alpha_model.params.acentricfactor.values
    m  = @. 0.480+1.547*ω-0.176*ω^2
    α  = @. (Tr>1)*(exp((1-2/(2+m))*(1-Tr^(1+m/2))))^2+(Tr<=1)*(1+m*(1-√(Tr)))^2
    return α
end

function α_function(model::PRModel,V,T,z,alpha_model::BMAlphaModel)
    Tc = model.params.Tc.values
    Tr = @. T/Tc
    ω  = alpha_model.params.acentricfactor.values
    #m  = @. 0.37464+1.54226*ω-0.26992*ω^2
    #α  = @. (Tr>1)*(exp((1-2/(2+m))*(1-Tr^(1+m/2))))^2+(Tr<=1)*(1+m*(1-√(Tr)))^2
    α = zeros(typeof(T),length(Tc))
    for i in @comps
        ωi = ω[i]
        m = evalpoly(ωi,(0.37464,1.54226,-0.26992))
        Tr = T/Tc[i]
        α[i] = ifelse(Tr>1,(exp((1-2/(2+m))*(1-Tr^(1+m/2))))^2,(1+m*(1-√(Tr)))^2)
    end
    return α
end