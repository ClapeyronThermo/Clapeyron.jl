include("equations.jl")



struct ESElectrolyte{T<:IdealModel,c<:EoSModel,i<:IonModel} <: ESElectrolyteModel
    components::Array{String,1}
    charge::Vector{Int64}
    idealmodel::T
    neutralmodel::c
    ionmodel::i
    references::Array{String,1}
end

"""
    ESElectrolyte(solvents::Array{String,1},
        ions::Array{String,1};
        idealmodel::IdealModel = BasicIdeal,
        neutralmodel::EoSModel = pharmaPCSAFT,
        ionmodel::IonModel = DH,
        RSPmodel::RSPModel = ConstRSP,
        charge::Vector{Int} = Int[],
        ideal_userlocations = String[],
        neutralmodel_userlocations = String[],
        ionmodel_userlocations = String[],
        RSPmodel_userlocations = String[],
        verbose::Bool=false)

## Description
This function provides the necessary framework to create an electrolyte model by combining ideal, neutral and ion models:
```julia
model = ESElectrolyte(["water"],["sodium","chloride"];
            idealmodel = BasicIdeal,
            neutralmodel = pharmaPCSAFT,
            ionmodel = DH,
            RSPmodel = ConstRSP)
```
Any of the available models in Clapeyron can be combined in the above. Note that neutral (solvent) species and ions are defined separately. Within Clapeyron, we will only support ion-based electrolyte models; as such, any salt-based approach (i.e. where the salt is treated as a separate species) will not be supported.
"""
function ESElectrolyte(solvents,ions;
    idealmodel = BasicIdeal,
    neutralmodel = pharmaPCSAFT,
    ionmodel = DH,
    RSPmodel = ConstRSP,
    charge = nothing,
    ideal_userlocations = String[],
    neutralmodel_userlocations = String[],
    ionmodel_userlocations = String[],
    RSPmodel_userlocations = String[],
    assoc_options = AssocOptions(),
    verbose = false,
    reference_state = nothing)

    raw_components = vcat(solvents,ions)
    formatted_components = format_components(raw_components)

    if isnothing(charge)
        charge_params = getparams(formatted_components, ["Electrolytes/properties/charges.csv"]; verbose=verbose)
        init_charge = charge_params["charge"].values

    elseif charge isa Vector{String}
        charge_params = getparams(formatted_components, ["Electrolytes/properties/charges.csv"]; userlocations=charge, verbose=verbose)
        init_charge = charge_params["charge"].values
    else
        init_charge = charge
    end

    #path0 = default_locations(neutralmodel)
    #remove unused datapaths
    #neutral_path = joinpath.(DB_PATH,filter(∉(("properties/molarmass.csv","properties/molarmass_groups.csv,properties/critical_csv")),path0))
    init_idealmodel = init_model(idealmodel,raw_components,ideal_userlocations,verbose)
    init_RSP = @initmodel RSPmodel(solvents,ions,userlocations = RSPmodel_userlocations,verbose = verbose)
    if has_sites(neutralmodel)
        init_neutralmodel = neutralmodel(raw_components;userlocations=neutralmodel_userlocations,verbose,assoc_options)
    else
        init_neutralmodel = neutralmodel(raw_components;userlocations=neutralmodel_userlocations,verbose)
    end

    init_ionmodel = @initmodel ionmodel(solvents,ions;RSPmodel=init_RSP,userlocations=ionmodel_userlocations,verbose=verbose)

    #components = init_neutralmodel.components

    references = String[]
    model = ESElectrolyte(formatted_components,init_charge,init_idealmodel,init_neutralmodel,init_ionmodel,references)
    set_reference_state!(model,reference_state;verbose)
    return model
end


function init_preferred_method(method::typeof(bubble_pressure),model::ESElectrolyteModel,kwargs)
    Z = model.charge
    nonvolatiles = [model.components[i] for i in @iions]
    return FugBubblePressure(;nonvolatiles = nonvolatiles,kwargs...)
end

function init_preferred_method(method::typeof(bubble_temperature),model::ESElectrolyteModel,kwargs)
    Z = model.charge
    nonvolatiles = [model.components[i] for i in @iions]
    return FugBubbleTemperature(;nonvolatiles = nonvolatiles,kwargs...)
end

function tp_flash_K0!(K,model::ESElectrolyteModel,p,T,z)
    Z = model.charge
    neutral = iszero.(Z)
    pures = split_model(model,neutral)
    psat = first.(extended_saturation_pressure.(pures,T))
    K .= 0
    Kview = @view K[neutral]
    Kview .= psat ./ p
    return K
end

function recombine_impl!(model::ESElectrolyteModel)
    recombine!(model.neutralmodel)
    recombine_ion!(model,model.ionmodel)
end

recombine_ion!(model,ionmodel) = recombine!(ionmodel)

function show_comps_with_charge(io,components,Z)
    function f(x)
        comp,zi = x
        if iszero(zi)
            return '"' * comp* '"'
        elseif zi < 0
            return '"' * comp * '"' * low_color(" (-" * string(abs(zi)) * ")")
        else
            return '"' * comp * '"' * low_color(" (+" * string(abs(zi)) * ")")
        end
    end
    show_pairs(io,map(f,zip(components,Z)),quote_string = false)
end

function Base.show(io::IO,mime::MIME"text/plain",model::ESElectrolyteModel)
    neutralmodel = model.neutralmodel
    print(io,"Explicit Electrolyte Model")
    if hasfield(typeof(neutralmodel),:components)
        length(neutralmodel) == 1 && println(io, " with 1 component:")
        length(neutralmodel) > 1 && println(io, " with ", length(neutralmodel), " components:")
        if has_groups(neutralmodel)
            groups = neutralmodel.groups
            show_groups(io,groups)
            #println(io)
            #print(io,"Group Type: ",groups.grouptype)
        else
            show_comps_with_charge(io,neutralmodel.components,model.charge)
        end
    else
        print(io,"()")
    end
    print(io,'\n',"Neutral Model: ",typeof(model.neutralmodel))
    print(io,'\n',"Ion Model: ",typeof(model.ionmodel))

    if requires_rsp(model.ionmodel)
        print(io,'\n',"RSP Model: ",typeof(model.ionmodel.RSPmodel))
    end

    show_reference_state(io,model;space = true)
end


include("stability.jl")

export dielectric_constant, ESElectrolyte

