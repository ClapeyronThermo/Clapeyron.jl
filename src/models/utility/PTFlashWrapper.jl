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
Base.eltype(model::PTFlashWrapper) = Base.promote_eltype(model.model,fug)

function tp_flash_K0(wrapper::PTFlashWrapper,p,T,z)
    first.(wrapper.sat) ./ p
end

function tp_flash_K0!(K,wrapper::PTFlashWrapper,p,T,z)
    K .=  first.(wrapper.sat) ./ p 
end

function PTFlashWrapper{TT}(model,equilibrium,pures = split_pure_model(model)) where TT
    nc = length(model)
    sat = Vector{Tuple{TT,TT,TT}}(undef,nc)
    ϕpure = Vector{TT}(undef,nc)
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

function update_temperature!(model::PTFlashWrapper,T)
    isnan(T) && return nothing
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

function __x0_bubble_pressure(model::PTFlashWrapper,T,x,y0 = nothing,volatiles = FillArrays.Fill(true,length(model)),pure = nothing,crit = nothing)
    sat = model.sat #saturation, we do not approximate here.
    p0r = first.(sat)
    p0 = index_expansion(p0r,volatiles)
    xipi = p0 .* x ./ sum(x)
    p0 = sum(xipi)
    if isnothing(y0)
        yx = xipi
        yx ./= p0
    else
        yx = y0
    end
    p,_,_,y,vl0,vv0 = improve_bubbledew_suggestion(model,p0,T,x,yx,FugEnum.BUBBLE_PRESSURE,volatiles,high_conditions)
    return p,vl0,vv0,y
end

function __x0_dew_pressure(model::PTFlashWrapper,T,y,x0=nothing,condensables = FillArrays.Fill(true,length(model)),pure = nothing, crit = nothing)
    sat = model.sat #saturation, we do not approximate here.
    p0inv_r = 1. ./ first.(sat)
    p0inv = index_expansion(p0inv_r,condensables)
    yipi = y .* p0inv ./ sum(y)
    p0 = 1/sum(yipi)
    if isnothing(x0)
        xx = yipi
        xx .*= p0
    else
        xx = x0
    end
    p,_,x,_,vl0,vv0 = improve_bubbledew_suggestion(model,p0,T,xx,y,FugEnum.DEW_PRESSURE,condensables,high_conditions)
    return p,vl0,vv0,x
end

function improve_bubbledew_suggestion(model::PTFlashWrapper,p0,T0,x,y,method,in_media,high_conditions)
    vl = volume(model,p,T,x,phase = :l)/sum(x)
    vv = volume(model,p,T,y,phase = :v)/sum(y)
    return p,T,x,y,vl,vv
end