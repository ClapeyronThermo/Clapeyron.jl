abstract type SAFTgammaEMieModel <: ElectrolyteSAFTModel end

struct SAFTgammaEMie{T<:IdealModel,c<:SAFTModel,i<:IonModel,b,r<:RSPModel} <: SAFTgammaEMieModel
    components::Array{String,1}
    solvents::Union{Array{String,1},Array{Any,1}}
    salts::Array{String,1}
    ions::Array{String,1}
    icomponents::UnitRange{Int}
    isolvents::UnitRange{Int}
    isalts::UnitRange{Int}
    iions::UnitRange{Int}
    stoic_coeff::Array{Float64}
    idealmodel::T
    saftmodel::c
    ionicmodel::i
    bornmodel::b
    rspmodel::r
    absolutetolerance::Float64
    references::Array{String,1}
end

@registermodel SAFTgammaEMie
export SAFTgammaEMie

function SAFTgammaEMie(solvents,salts,ions; saftmodel=SAFTgammaMie,
    idealmodel = BasicIdeal,
    ionicmodel = GCMSA,
    RSPmodel = Schreckenberg,
    bornmodel = GCBorn,
    userlocations=String[], 
    ideal_userlocations=String[],
     verbose=false)
    
    salt_groups = GroupParam(salts, ["Electrolytes/properties/salts.csv"]; verbose=verbose)
    solvent_groups = GroupParam(solvents, ["SAFT/SAFTgammaMie/SAFTgammaMie_groups.csv"]; verbose=verbose)

    nsolv = length(solvents)
    nsalts = length(salts)
    nions = length(ions)
    
    stoichiometric_coeff = zeros(length(salts),length(ions))
    for i in 1:length(salts)
        stoichiometric_coeff[i,:] = salt_groups.n_flattenedgroups[i]
    end

    components = deepcopy(salt_groups.components)
    prepend!(components,solvent_groups.components)

    isolvents = 1:nsolv
    iions = (nsolv+1):(nsolv+nions)
    isalts = (nsolv+1):(nsolv+nsalts)

    init_idealmodel = init_model(idealmodel,cat(solvents,ions,dims=1),ideal_userlocations,verbose)
    init_SAFTmodel = saftmodel(cat(solvents,ions,dims=1))
    init_Ionicmodel = ionicmodel(solvents,salts,ions;RSPmodel=nothing,SAFTlocations=["SAFT/"*string(saftmodel)])
    if bornmodel !== nothing
        init_bornmodel = bornmodel(solvents,salts,ions;RSPmodel=nothing,SAFTlocations=["SAFT/"*string(saftmodel)])
    else
        init_bornmodel = nothing
    end
    init_RSPmodel = RSPmodel(solvents,salts)

    solvents = solvent_groups.components
    salts = salt_groups.components
    ions = salt_groups.flattenedgroups
    components = deepcopy(salts)
    prepend!(components,solvents)
    icomponents = 1:length(components)

    references = String[]
    model = SAFTgammaEMie(components,solvents,salts,ions,icomponents,isolvents,isalts,iions,stoichiometric_coeff,init_idealmodel,init_SAFTmodel,init_Ionicmodel,init_bornmodel,init_RSPmodel,1e-12,references)
    return model
end

function a_res(model::SAFTgammaEMieModel, V, T, z,_data=@f(data))
    (data_saft,data_ion) = _data
    return a_res(model.saftmodel,V,T,z,data_saft)+a_res(model.ionicmodel,V,T,z,data_ion)+a_res(model.bornmodel,V,T,z,data_ion)
end

function data(model::SAFTgammaEMieModel, V, T, z)
    _data_saft = data(model.saftmodel,V,T,z)
    _data_msa = data_msa(model.ionicmodel,V,T,z)
    _data_rsp = dielectric_constant(model.rspmodel,V,T,z,_data_saft)
    return (_data_saft,(_data_msa,_data_rsp))
end

function x0_volume(model::SAFTgammaEMie,p,T,z; phase = :unknown)
    phase = Symbol(phase)
    if phase === :unknown || is_liquid(phase)
        return x0_volume_liquid(model.saftmodel,T,z)
    elseif is_vapour(phase)
        return x0_volume_gas(model.saftmodel,p,T,z)
    elseif is_supercritical(phase)
     else
        error("unreachable state on x0_volume")
    end
end

function x0_sat_pure(model::SAFTgammaEMieModel,T,z=[0.5,0.5])
    vl = x0_volume(model,1e5,T,z,phase=:l)
    return (log10(vl),10)
end
function scale_sat_pure(model::SAFTgammaEMieModel,z=[0.5,0.5])
    p    = 1/p_scale(model.saftmodel,z)
    μ    = 1/R̄/T_scale(model.saftmodel,z)
    return p,μ
end