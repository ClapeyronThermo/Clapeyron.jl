#wrapper used to cache results in case of activity models and CompositeModel
struct PTFlashWrapper{T,T2,R,S} <: EoSModel
    components::Vector{String}
    model::T
    pures::T2
    sat::Vector{R}
    fug::Vector{S}
    equilibrium::Symbol
end

Base.length(model::PTFlashWrapper) = length(model.model)

function tp_flash_K0(wrapper::PTFlashWrapper,p,T,z)
    first.(wrapper.sat) ./ p
end

function tp_flash_K0!(K,wrapper::PTFlashWrapper,p,T,z)
    K .=  first.(wrapper.sat) ./ p 
end

function PTFlashWrapper{TT}(model,equilibrium::Symbol,pures = split_pure_model(model)) where TT
    nc = length(model)
    sat = Vector{Tuple{TT,TT,TT}}(undef,nc)
    ϕpure = Vector{TT}(undef,nc)
    return PTFlashWrapper(component_list(model),model,pures,sat,ϕpure,equilibrium)
end

function PTFlashWrapper(model,equilibrium::Symbol,pures = split_pure_model(model))
    TT = Base.promote_eltype(model,Float64)
    return PTFlashWrapper{TT}(model,equilibrium,pures)
end

function PTFlashWrapper(model,p,T,equilibrium::Symbol)
    wrapper = PTFlashWrapper(model,equilibrium)
    update_temperature!(wrapper,T)
    return wrapper
end

function update_temperature!(model::PTFlashWrapper,T)
    pures = model.pures
    ϕ = model.fug
    sats = model.sat
    TT = eltype(ϕ)
    for i in 1:length(model)
        pure = pures[i]
        sat = saturation_pressure(pure,T)
        ps,vl,vv = sat
        sats[i] = sat
        if gas_model(pure) isa IdealModel
            ϕ[i] = 1.0
        else
            ϕ[i] = exp(VT_lnϕ_pure(gas_model(pure),vv,T,ps))
        end
    end
    return nothing
end

gas_model(model::PTFlashWrapper) = gas_model(model.model)

function volume_impl(model::PTFlashWrapper, p, T, z, phase, threaded, vol0)
    volume_impl(model.model, p, T, z, phase, threaded, vol0)
end

function Base.show(io::IO,mime::MIME"text/plain",wrapper::PTFlashWrapper)
    model = wrapper.model
    pure = wrapper.pures
    print(io,"PT-Flash Wrapper")
    length(model) == 1 && print(io, " with 1 component:")
    length(model) > 1 && print(io, " with ", length(model), " components:")
    println(io)
    show_pairs(io,wrapper.components)
    print(io,'\n',"Mixture model: ", typeof(model))
    print(io,'\n',"Pure model: ",eltype(pure))
    print(io,'\n',"Equilibrium type: :",wrapper.equilibrium)
    show_reference_state(io,model;space = true)
end

