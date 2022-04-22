abstract type SchreckenbergModel <: RSPModel end

struct SchreckenbergParam <: EoSParam
    d_T::SingleParam{Float64}
    d_V::SingleParam{Float64}
end

struct Schreckenberg <: SchreckenbergModel
    components::Array{String,1}
    solvents::Array{String,1}
    ions::Array{String,1}
    icomponents::UnitRange{Int}
    isolvents::UnitRange{Int}
    iions::UnitRange{Int}
    params::SchreckenbergParam
    absolutetolerance::Float64
    references::Array{String,1}
end

@registermodel Schreckenberg
export Schreckenberg
function Schreckenberg(solvents,salts; userlocations::Vector{String}=String[], verbose::Bool=false)
    ion_groups = GroupParam(salts, ["Electrolytes/properties/salts.csv"]; verbose=verbose)

    ions = ion_groups.flattenedgroups
    components = deepcopy(solvents)
    append!(components,ions)
    icomponents = 1:length(components)
    isolvents = 1:length(solvents)
    iions = (length(solvents)+1):length(components)

    params = getparams(solvents, ["Electrolytes/properties/Schreckenberg.csv"]; userlocations=userlocations, verbose=verbose)
    d_T = params["d_T"]
    d_V = params["d_V"]
    packagedparams = SchreckenbergParam(d_T,d_V)

    references = [""]
    
    model = Schreckenberg(components, solvents, ions, icomponents, isolvents, iions, packagedparams, 1e-12,references)
    return model
end

function RSP(electromodel::ElectrolyteModel,V,T,z,model::SchreckenbergModel)
        d_T = model.params.d_T.values
        d_V = model.params.d_V.values

        d = @. d_V*(d_T/T-1)

        d = (d .+d')/2

        n_solv = sum(z[i] for i ∈ model.isolvents)
        ρ_solv = n_solv / V

        x0 = z ./ n_solv 

        d̄ = sum(sum(x0[i]*x0[j]*d[i,j] for j ∈ model.isolvents) for i ∈ model.isolvents)
    return 1+ρ_solv*d̄
end

is_splittable(::Schreckenberg) = false