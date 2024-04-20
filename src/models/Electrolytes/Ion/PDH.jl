struct PDHParam <: EoSParam
    charge::SingleParam{Float64}
    Mw::SingleParam{Float64}
end

abstract type PDHModel <: IonModel end

struct PDH{ϵ} <: PDHModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    params::PDHParam
    RSPmodel::ϵ
    solvent_density::Float64
    references::Array{String,1}
end

@registermodel PDH

"""
    PDH(solvents,ions; 
    RSPmodel=ConstW, 
    SAFTlocations=String[], 
    userlocations=String[], 
    ideal_userlocations=String[], 
    verbose=false)

Pitzer-Debye-Hückel (PDH) model for electrostatic interaction.
"""
PDH

export PDH
function PDH(solvents,ions; solvent_density = 1000., RSPmodel=ConstRSP, SAFTlocations=String[], userlocations=String[], verbose=false)
    components = deepcopy(ions)
    prepend!(components,solvents)
    icomponents = 1:length(components)

    params = getparams(components, ["Electrolytes/properties/charges.csv","properties/molarmass.csv"]; userlocations=userlocations, verbose=verbose)

    
    charge = params["charge"]
    Mw = params["Mw"]

    packagedparams = PDHParam(charge,Mw)

    references = String[]

    init_RSPmodel = init_electrolyte_model(RSPmodel,solvents,ions)

    model = PDH(components, icomponents, packagedparams, init_RSPmodel, solvent_density,references)
    return model
end

function data(model::PDHModel, V, T, z)
    return dielectric_constant(model.RSPmodel, V, T, z)
end

function excess_g_res(model::PDHModel, V, T, z, _data=@f(data))
    ϵ_r = _data[1]
    Z = model.params.charge.values
    Mw = model.params.Mw.values
    ρsol = model.solvent_density

    iions = model.icomponents[Z.!=0]
    isolv = model.icomponents[Z.==0]
    if length(iions) == 0
        return zero(V+T+first(z))
    end
    ρ = 14.9
    ∑z = sum(z)

    Ix = 0.5*sum(Z[i]^2*z[i] for i in iions)/∑z

    Msol = sum(Mw[i]*z[i] for i in isolv)/sum(z[i] for i in isolv)/1000
    Aϕ = 1/3*√(2π*N_A*ρsol)*(e_c^2/(4π*ϵ_0*ϵ_r*k_B*T))^1.5
    println("Aϕ = ", Aϕ)
    return -1/√(Msol)*4*Aϕ*Ix/ρ*log1p(ρ*√(Ix))*∑z*R̄*T
end

excess_gibbs_free_energy(model::PDHModel,p,T,z) = excess_g_res(model,p,T,z)