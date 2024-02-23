abstract type ePCSAFTModel <: ESElectrolyteModel end

struct ePCSAFT{T<:IdealModel,c<:EoSModel,i<:IonModel} <: ePCSAFTModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    charge::Vector{Int64}
    idealmodel::T
    neutralmodel::c
    ionmodel::i
    references::Array{String,1}
end

function ePCSAFT(solvents,ions; 
    idealmodel = BasicIdeal,
    neutralmodel = pharmaPCSAFT,
    ionmodel = DH,
    RSPmodel = ConstRSP,
    userlocations=String[], 
    ideal_userlocations=String[],
    assoc_options = nothing,
     verbose=false)
    components = deepcopy(ions)
    prepend!(components,solvents)

    params = getparams(format_components(components), ["Electrolytes/properties/charges.csv"]; userlocations=userlocations, verbose=verbose)
    _charge = params["charge"]
    charge = _charge.values

    icomponents = 1:length(components)

    neutral_path = DB_PATH.*["/SAFT/PCSAFT","/SAFT/PCSAFT/pharmaPCSAFT"]

    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    init_neutralmodel = neutralmodel(components;userlocations=userlocations,verbose=verbose,assoc_options=assoc_options)
    init_ionmodel = ionmodel(solvents,ions;RSPmodel=RSPmodel,userlocations=append!(userlocations,neutral_path),verbose=verbose)


    for i in ions
        init_neutralmodel.params.epsilon[i] = 0. #pure ion has ϵi 
        for j in ions
            if sign(_charge[i]) == sign(_charge[j]) #cation-cation and anion-anion interactions are neglected.
                init_neutralmodel.params.epsilon[i,j] = 0.
            end
        end
    end



    references = String[]
    components = format_components(components)
    model = ePCSAFT(components,icomponents,charge,init_idealmodel,init_neutralmodel,init_ionmodel,references)
    return model
end

function a_res(model::ePCSAFTModel, V, T, z)
    data_pcsaft = data(model.neutralmodel,V,T,z)
    data_ion = data(model.ionmodel,V,T,z)
    data_ion = (data_ion[1],data_pcsaft[1])

    return a_res(model.neutralmodel,V,T,z,data_pcsaft)+a_res(model.ionmodel,V,T,z,data_ion)
end

export ePCSAFT