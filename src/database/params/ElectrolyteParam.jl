struct ElectrolyteParam <: ClapeyronParam
    components::Array{String,1}
    solvents::Union{Array{String,1},Array{Any,1}}
    salts::Array{String,1}
    ions::Array{String,1}
    icomponents::UnitRange{Int}
    isolvents::UnitRange{Int} #
    isalts::UnitRange{Int} #
    iions::UnitRange{Int} #
    stoic_coeff::Array{Float64,2}
end

function ElectrolyteParam(solventdata,salt_groups::GroupParam,iondata)
    
    if iondata isa GroupParam
        ions = iondata.components
    else
        ions = iondata
    end

    if solventdata isa GroupParam
        solvents = solventdata.components
    else
        solvents = solventdata
    end

    salts = salt_groups.components
    solvents = solvent_groups.components
    isolvents = 1:length(solvents)
    iions = (length(solvents)+1):(length(solvents)+length(ions))
    isalts = (length(solvents)+1):(length(solvents)+length(salts))
    
    components = vcat(solvents,salts)

    stoichiometric_coeff = zeros(length(salts),length(ions))
    for i in 1:length(salts)
        stoichiometric_coeff[i,:] = salt_groups.n_flattenedgroups[i]
    end

    return ElectrolyteParam(components,
    solvents,
    salts,
    ions,
    1:length(components),
    isolvents,
    isalts,
    iions,
    stoichiometric_coeff
    )
end