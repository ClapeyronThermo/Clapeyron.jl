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

"""
    Schreckenberg(solvents::Array{String,1}, 
         ions::Array{String,1}; 
         userlocations::Vector{String}=[], 
         verbose::Bool=false)

## Input parameters
- `d_T::Float64`: Single Parameter - Temperature dependent dielectric constant `[-]`
- `d_V::Float64`: Single Parameter - Volume dependent dielectric constant `[-]`
- `charge::Float64`: Single Parameter - Charge `[-]`

## Description
This function is used to create a Schreckenberg model. The Schreckenberg term estimates the dielectric constant for a mixture of solvents.

## References
1. Schreckenberg, J., Dufal, S., Haslam, A.J., Adjiman, C.S., Jackson, G., Galindo, A. (2014). Modelling of the thermodynamic and solvation properties of electrolyte solutions with the statistical associating fluid theory for potentials of variable range. Molecular Physics, 112(17), 2339-2364.
"""
function Schreckenberg(solvents,ions; userlocations::Vector{String}=String[], verbose::Bool=false)
    components = deepcopy(ions)
    prepend!(components,solvents)
    components = format_components(components)
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

        if isempty(ineutral)
            return 1.0
        end
        
        n_solv = zero(first(z))
        for i in ineutral
            n_solv += z[i]
        end
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