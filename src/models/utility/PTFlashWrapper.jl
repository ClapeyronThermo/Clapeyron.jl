#wrapper used to cache results in case of activity models and CompositeModel
#=
struct PTFlashWrapper{T,T2,R,S} <: EoSModel
    components::Vector{String}
    model::T
    pures::T2
    sat::Vector{R}
    fug::Vector{S}
    equilibrium::Symbol
end =#

Base.length(Base.@specialize(model::PTFlashWrapper)) = length(model.model)
Base.eltype(Base.@specialize(model::PTFlashWrapper)) = Base.promote_eltype(model.model,model.fug)
__γ_unwrap(Base.@specialize(model::PTFlashWrapper)) = __γ_unwrap(model.model)
@inline gas_model(Base.@specialize(model::PTFlashWrapper)) = gas_model(model.model)

function tp_flash_K0(wrapper::PTFlashWrapper,p,T,z)
    first.(wrapper.sat) ./ p
end

function tp_flash_K0!(K,wrapper::PTFlashWrapper,p,T,z)
    K .= first.(wrapper.sat) ./ p
end

function K0_lle_init(wrapper::PTFlashWrapper,p,T,z,cache = tpd_cache(wrapper,p,T,z);reduced = true)
    return K0_lle_init(__γ_unwrap(wrapper),p,T,z,cache;reduced)
end

function PTFlashWrapper{TT}(model,equilibrium,pures = split_pure_model(model)) where TT
    nc = length(model)
    sat = Vector{Tuple{TT,TT,TT}}(undef,nc)
    nan = TT(NaN)
    fill!(sat,(nan,nan,nan))
    ϕpure = Vector{TT}(undef,nc)
    fill!(ϕpure,nan)
    return PTFlashWrapper(component_list(model),model,pures,sat,ϕpure,equilibrium)
end

function PTFlashWrapper(model,equilibrium,pures = split_pure_model(model))
    TT = Base.promote_eltype(model,Float64)
    return PTFlashWrapper{TT}(model,equilibrium,pures)
end

function PTFlashWrapper(model,p,T,z,equilibrium)
    TT = Base.promote_eltype(model,p,T,z)
    wrapper = PTFlashWrapper{TT}(model,equilibrium)
    update_temperature!(wrapper,T)
    return wrapper
end

split_pure_model(model::PTFlashWrapper,splitter) = model.pures[splitter]
idealmodel(model::PTFlashWrapper) = idealmodel(model.model)
fluid_model(model::PTFlashWrapper) = fluid_model(model.model)
molecular_weight(model::PTFlashWrapper,z) = molecular_weight(model.model,z)
reference_state(model::PTFlashWrapper) = reference_state(model.model)

function update_temperature!(model::PTFlashWrapper,T)
    isnan(T) && return nothing
    pures = model.pures
    lnϕ = model.fug
    sats = model.sat
    TT = eltype(lnϕ)
    for i in 1:length(model)
        pure = pures[i]
        sat = saturation_pressure(pure,T)
        ps,vl,vv = sat
        sats[i] = sat
        if gas_model(pure) isa IdealModel
            lnϕ[i] = 0.0
        else
            lnϕ[i] = VT_lnϕ_pure(gas_model(pure),vv,T,ps)
        end
    end
    return nothing
end

function _update_temperature_with_view!(model1::TT,model2::TT,T,_view) where TT <: PTFlashWrapper
    n1,n2 = length(model1),length(model2)
    if n1 > n2
        update_temperature!(model1,T)
        model2.sat .= @view model1.sat[_view]
        model2.fug .= @view model1.fug[_view]
    else
        update_temperature!(model2,T)
        model1.sat .= @view model2.sat[_view]
        model1.fug .= @view model2.fug[_view]
    end
    return nothing
end

_update_temperature_with_view!(model1,model2,T,_view) = nothing

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

function saturation_pressure_ad2(result,model,T)
    return saturation_pressure_ad(result,(model,T),(model,primalval(T)))
end

include("PTFlashWrapper/PT.jl")
include("PTFlashWrapper/fugacity.jl")
include("PTFlashWrapper/bubbledew.jl")