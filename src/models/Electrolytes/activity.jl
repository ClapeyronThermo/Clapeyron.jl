abstract type ElectrolyteActivityModel <: ElectrolyteModel end

struct ElectrolyteActivity{T<:IdealModel,c<:ActivityModel,i<:IonModel} <: ElectrolyteActivityModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    charge::Vector{Int64}
    idealmodel::T
    neutralmodel::c
    ionmodel::i
    references::Array{String,1}
end

function ElectrolyteActivity(solvents,ions; 
    idealmodel = BasicIdeal,
    neutralmodel = NRTL,
    ionmodel = PDH,
    RSPmodel = ConstRSP,
    userlocations=String[], 
    ideal_userlocations=String[],
     verbose=false)
    components = deepcopy(ions)
    prepend!(components,solvents)

    params = getparams(components, ["Electrolytes/properties/charges.csv"]; userlocations=userlocations, verbose=verbose)
    charge = params["charge"].values

    icomponents = 1:length(components)

    neutral_path = [DB_PATH*"/"*default_locations(neutralmodel)[1]]

    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    init_neutralmodel = neutralmodel(components;userlocations=userlocations,verbose=verbose)
    init_ionmodel = ionmodel(solvents,ions;RSPmodel=RSPmodel,userlocations=append!(userlocations,neutral_path),verbose=verbose)

    references = String[]
    model = ElectrolyteActivity(components,icomponents,charge,init_idealmodel,init_neutralmodel,init_ionmodel,references)
    return model
end

export ElectrolyteActivity

function excess_gibbs_free_energy(model::ElectrolyteActivityModel, V, T, z)
    return excess_gibbs_free_energy(model.neutralmodel,V,T,z)+excess_gibbs_free_energy(model.ionmodel,V,T,z)
end

function activity_coefficient(model::ElectrolyteActivityModel,p,T,z)
    ion = model.charge.!=0

    X = gradient_type(p,T,z)
    lnγres = (Solvers.gradient(x->excess_gibbs_free_energy(model.neutralmodel,p,T,x),z)/(R̄*T))::X
    lnγion = (Solvers.gradient(x->excess_gibbs_free_energy(model.ionmodel,p,T,x),z)/(R̄*T))::X

    z0 = ones(length(z))*1e-30
    @. z0[!ion] = z[!ion]
    z0 = z0./sum(z0)

    lnγres0 = (Solvers.gradient(x->excess_gibbs_free_energy(model.neutralmodel,p,T,x),z0)/(R̄*T))::X

    lnγ = @. ion*(lnγres-lnγres0)+!ion*lnγres+lnγion

    return exp.(lnγ)
end


function mean_ionic_activity_coefficient(model::ElectrolyteActivityModel,salts,p,T,m,zsolvent=[1.])
    isolvent = model.icomponents[model.charge.==0]
    iions = model.icomponents[model.charge.!=0]

    ν = salt_stoichiometry(model,salts)

    z = molality_to_composition(model,salts,m,zsolvent)

    γim = activity_coefficient(model,p,T,z)[iions].*sum(z[isolvent])/sum(z)
    γsm = (prod(γim'.^ν,dims=2)).^(1 ./sum(ν,dims=2))
    return γsm
end