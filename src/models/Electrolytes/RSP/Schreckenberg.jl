abstract type SchreckenbergModel <: RSPModel end

struct SchreckenbergParam <: EoSParam
    d_T::SingleParam{Float64}
    d_V::SingleParam{Float64}
    charge::SingleParam{Float64}
end

struct Schreckenberg <: SchreckenbergModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    params::SchreckenbergParam
    references::Array{String,1}
end

@registermodel Schreckenberg
export Schreckenberg
function Schreckenberg(solvents,ions; userlocations::Vector{String}=String[], verbose::Bool=false)
    components = deepcopy(ions)
    prepend!(components,solvents)
    icomponents = 1:length(components)

    params = getparams(components, ["Electrolytes/RSP/Schreckenberg.csv","Electrolytes/properties/charges.csv"]; userlocations=userlocations, verbose=verbose, ignore_missing_singleparams=["d_T","d_V"])
    d_T = params["d_T"]
    d_V = params["d_V"]
    charge = params["charge"]
    packagedparams = SchreckenbergParam(d_T,d_V,charge)

    references = String[]
    
    model = Schreckenberg(components,icomponents,packagedparams,references)
    return model
end

function dielectric_constant(model::SchreckenbergModel,V,T,z,_data=nothing)
        d_T = model.params.d_T.values
        d_V = model.params.d_V.values
        Z = model.params.charge.values
        ineutral = model.icomponents[Z.==0]

        n_solv = sum(z[ineutral])
        ρ_solv = n_solv / V
        d̄ = zero(T+first(z))

        for i in ineutral
            di = d_V[i]*(d_T[i]/T-1)
            dij,zi = di,z[i]
            d̄ += dij*zi*zi
            for j in ineutral[ineutral.!=i]
                dj = d_V[j]*(d_T[j]/T-1)
                dij,zj = 0.5*(di+dj),z[j]
                d̄ += dij*zi*zj
            end
        end

    d̄ = d̄/(n_solv*n_solv)
    return 1+ρ_solv*d̄
end