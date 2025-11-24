include("equations.jl")

abstract type ESElectrolyteModel <: ElectrolyteModel end

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
        charges = String[],
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
    charges = String[],
    ideal_userlocations = String[],
    neutralmodel_userlocations = String[],
    ionmodel_userlocations = String[],
    RSPmodel_userlocations = String[],
    assoc_options = AssocOptions(),
    verbose = false,
    reference_state = nothing)

    solvents = format_components(solvents)
    ions = format_components(ions)
    components = deepcopy(ions)
    prepend!(components,solvents)

    params = getparams(components, ["Electrolytes/properties/charges.csv"]; userlocations=charges, verbose=verbose)
    charge = params["charge"].values
    #path0 = default_locations(neutralmodel)
    #remove unused datapaths
    #neutral_path = joinpath.(DB_PATH,filter(∉(("properties/molarmass.csv","properties/molarmass_groups.csv,properties/critical_csv")),path0))
    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    init_RSP = @initmodel RSPmodel(solvents,ions,userlocations = RSPmodel_userlocations,verbose = verbose)
    if has_sites(neutralmodel)
        init_neutralmodel = neutralmodel(components;userlocations=neutralmodel_userlocations,verbose,assoc_options)
    else
        init_neutralmodel = neutralmodel(components;userlocations=neutralmodel_userlocations,verbose)
    end

    init_ionmodel = @initmodel ionmodel(solvents,ions;RSPmodel=init_RSP,userlocations=ionmodel_userlocations,verbose=verbose)

    components = init_neutralmodel.components

    references = String[]
    model = ESElectrolyte(components,charge,init_idealmodel,init_neutralmodel,init_ionmodel,references)
    set_reference_state!(model,reference_state;verbose)
    return model
end

#=
Taking an inspiration from the broadcast dispatch

Electrolyte models are deeply interwined.
We need to know the level at which they are interwined.

1. the charge parameter is shared between all models -> we need to pass charge to all inner models
2. some ion models use the EoS molecular size, some use their own -> we need to dispatch on that
3. some ion models require rsp, some not -> we need to dispatch on that too

I decided that there are three levels of interwining:

level 1: Electrolyte model provides charges to ion model, the ion model then can call itself
level 2: Electrolyte models provided

=#

abstract type IonDependency end
struct IndependentIonModel <: IonDependency end
struct DependentIonModel{T} <: IonDependency
    model::T
end

function IonDependency(model::ESElectrolyteModel)
    return IonDependency(model.ionmodel)
end

function IonDependency(model::IonModel)
    return IonDependency(model.RSPmodel)
end

IonDependency(model::RSPModel) = IndependentIonModel()

requires_rsp(::Type{T}) where T <: IonModel = _requires_rsp(T)
requires_rsp(model::IonModel) = _requires_rsp(typeof(model))
Base.@assume_effects :foldable function _requires_rsp(::Type{T}) where T
    return hasfield(T,:RSPmodel)
end

has_sigma(::Type{T}) where T <: IonModel = _has_sigma(T)
has_sigma(model::IonModel) = _has_sigma(typeof(model))
Base.@assume_effects :foldable function _has_sigma(::Type{T}) where T
    if hasfield(T,:params)
        P = fieldtype(T,:params)
        return hasfield(P,:sigma)
    else
        return false
    end
end

function get_sigma(ionmodel::IonModel, V, T, z, model, neutral_data = @f(data))
    if has_sigma(ionmodel)
        return ionmodel.params.sigma.values
    end

    if model isa CPAModel
        b = model.cubicmodel.params.b.values
        σ = similar(b,length(neutralmodel))
        for i in 1:length(neutralmodel)
            σ[i] = cbrt((3/2/N_A/π)*b[i,i])
        end
    else
        σ = diagvalues(model.params.sigma.values)
    end
    return σ
end

function a_res(model::ESElectrolyteModel, V, T, z)
    return a_res(model,V,T,z,IonDependency(model))
end

function iondata(model::ESElectrolyteModel,V,T,z)
    iondata(model::ESElectrolyteModel,V,T,z,IonDependency(model))
end

function iondata(model::ESElectrolyteModel, V, T, z, m::IndependentIonModel)
    neutralmodel = model.neutralmodel
    ionmodel = model.ionmodel
    neutral_data = data(neutralmodel,V,T,z)
    Z = model.charge
    σ = get_sigma(ionmodel, V, T, z, neutralmodel, neutral_data) #sigma is stored in the ionmodel
    if requires_rsp(ionmodel)
        ϵ_r = dielectric_constant(ionmodel, V, T, z, Z, m)
    else
        ϵ_r = one(Base.promote_eltype(ionmodel, V, T, z, Z))
    end
    return (Z, σ, ϵ_r)
end

function a_res(model::ESElectrolyteModel, V, T, z, m::IndependentIonModel)
    neutralmodel = model.neutralmodel
    ionmodel = model.ionmodel
    neutral_data = data(neutralmodel,V,T,z)
    Z = model.charge
    σ = get_sigma(ionmodel, V, T, z, neutralmodel, neutral_data) #sigma is stored in the ionmodel
    if requires_rsp(ionmodel)
        ϵ_r = dielectric_constant(ionmodel, V, T, z, Z, m)
    else
        ϵ_r = one(Base.promote_eltype(ionmodel, V, T, z, Z))
    end
    iondata = (Z, σ, ϵ_r)
    return a_res(neutralmodel, V, T, z, neutral_data) + a_res(ionmodel, V, T, z, iondata, neutralmodel, neutral_data)
end

function a_res(ionmodel::IonModel, V, T, z, iondata, neutralmodel, neutral_data)
    return a_res(ionmodel, V, T, z, iondata)
end

function lb_volume(model::ESElectrolyteModel,z)
    return lb_volume(model.neutralmodel,z)
end

function lb_volume(model::ESElectrolyteModel,T,z)
    return lb_volume(model.neutralmodel,T,z)
end

function x0_volume_liquid(model::ESElectrolyteModel,p,T,z)
    return x0_volume_liquid(model.neutralmodel,p,T,z)*1.15
end

function x0_volume_gas(model::ESElectrolyteModel,p,T,z)
    return x0_volume_gas(model.neutralmodel,p,T,z)
end

function mw(model::ElectrolyteModel)
    return mw(model.neutralmodel)
end

function p_scale(model::ElectrolyteModel,z)
    return p_scale(model.neutralmodel,z)
end

function T_scale(model::ElectrolyteModel,z)
    return T_scale(model.neutralmodel,z)
end

function a_born(model::ESElectrolyteModel,V,T,z)
    return a_born(model.ionmodel,V,T,z)
end

function debye_length(model::ESElectrolyteModel,V,T,z,ϵ_r = @f(dielectric_constant),∑z = sum(z))
    Z = model.charge
    return debye_length(V,T,z,ϵ_r,Z)
end

"""
    dielectric_constant(model::ElectrolyteModel, V, T, z)

Calculates the dielectric constant (also known as relative static permittivity) for a given electrolyte model.

## Examples
```julia
model = ConstRSP()
εr = dielectric_constant(model, 1.8e-5, 298.15, [1.0])
```
"""
function dielectric_constant end

function dielectric_constant(model::EoSModel,V, T, z)
    Z = model.charge
    return dielectric_constant(model, V, T, z, Z, IonDependency(model))
end

function dielectric_constant(model::ESElectrolyteModel,V, T, z, Z, ::IndependentIonModel)
    return dielectric_constant(model.ionmodel, V, T, z, Z, IonDependency(model.ionmodel))
end

function dielectric_constant(model::IonModel, V, T, z, Z, ::IndependentIonModel)
    return dielectric_constant(model.RSPmodel, V, T, z, Z, IonDependency(model.RSPmodel))
end

function dielectric_constant(model::RSPModel, V, T, z, Z, ::IndependentIonModel)
    return dielectric_constant(model, V, T, z, Z)
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

