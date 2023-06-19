abstract type EmpiricIdealModel <: IdealModel end

struct EmpiricIdeal <: EmpiricIdealModel
    components::Vector{String}
    params::MultiFluidParam
    pures::Vector{SingleFluidIdeal}
    Rgas::Float64
    references::Vector{String}
end

struct EmpiricIdealfromMulti{𝔸} <: EmpiricIdealModel
    components::Vector{String}
    params::MultiFluidParam
    pures::Vector{SingleFluid{𝔸}}
    Rgas::Float64
    references::Vector{String}
end

Rgas(m::EmpiricIdealModel) = m.Rgas

function EmpiricIdeal(components;
    pure_userlocations = String[],
    estimate_pure = false,
    coolprop_userlocations = true,
    Rgas = R̄,
    verbose = false,
    )
    
    pures = [
        SingleFluidIdeal(comp;
        userlocations = pure_userlocations,
        verbose = verbose, 
        estimate_pure = estimate_pure, 
        coolprop_userlocations = coolprop_userlocations,
        ) 
        for comp in components]
    params = MultiFluidParam(components,pures)
    references = unique!(reduce(vcat,pure.references for pure in pures))
    model = EmpiricIdeal(components,params,pures,Rgas,references)
    return model
end

function idealmodel(m::MultiFluid)
    EmpiricIdealfromMulti(m.components,m.params,m.pures,m.Rgas,m.references)
end



