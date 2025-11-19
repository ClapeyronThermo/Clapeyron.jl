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

Base.length(model::PTFlashWrapper) = length(model.model)
Base.eltype(model::PTFlashWrapper) = Base.promote_eltype(model.model)
__γ_unwrap(model::PTFlashWrapper) = __γ_unwrap(model.model)
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
    p,_,_,y,vl0,vv0 = improve_bubbledew_suggestion(model,p0,T,x,yx,FugEnum.BUBBLE_PRESSURE,volatiles,false)
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
    p,_,x,_,vl0,vv0 = improve_bubbledew_suggestion(model,p0,T,xx,y,FugEnum.DEW_PRESSURE,condensables,false)
    return p,vl0,vv0,x
end

function __x0_bubble_temperature(model::PTFlashWrapper,p,x,Tx0 = nothing,volatiles = FillArrays.Fill(true,length(model)),pure = nothing,crit = nothing)
    x_r = @view x[volatiles]
    pure = @view model.pures[volatiles]
    sat = @view model.sat[volatiles]
    if Tx0 !== nothing
        T0 = Tx0
        for i in 1:length(pure)
            sat[i] = saturation_pressure(pure[i],T0)
        end
        p_i_r = first.(sat)
    else
        dPdTsat = extended_dpdT_temperature.(pure,p,crit)
        T0 = antoine_bubble_solve(dPdTsat,p,x_r)
        p_i_r = antoine_pressure.(dPdTsat,T0)
    end
    xipi_r = y_r = p_i_r .* x_r ./ sum(x_r)
    p = sum(xipi_r)
    y_r ./= p
    y0 = index_expansion(y_r,volatiles)
    _,T,_,y,vl0,vv0 = improve_bubbledew_suggestion(model,p,T0,x,y0,FugEnum.BUBBLE_TEMPERATURE,volatiles,false)
    update_temperature!(model,T)
    return T,vl0,vv0,y
end

function __x0_dew_temperature(model::EoSModel,p,y,Tx0 = nothing,condensables = FillArrays.Fill(true,length(model)),pure = split_pure_model(model,condensables),crit = nothing)
    y_r = @view y[condensables]
    pure = @view model.pures[condensables]
    sat = @view model.sat[condensables]
    if Tx0 !== nothing
        T0 = Tx0
        for i in 1:length(pure)
            sat[i] = saturation_pressure(pure[i],T0)
        end
        p0inv_r = 1 ./ first.(sat)
    else
        dPdTsat = extended_dpdT_temperature.(pure,p,crit)
        T0 = antoine_bubble_solve(dPdTsat,p,y_r)
        p0inv_r = 1 ./ antoine_pressure.(dPdTsat,T0)
    end
    yipi_r = x_r = y_r .* p0inv_r ./ sum(y_r)
    p_r = 1/sum(yipi_r)
    x_r .*= p_r
    x0 = index_expansion(x_r,condensables)
    update_temperature!(model,T0)
    _,T,x,_,vl0,vv0 = improve_bubbledew_suggestion(model,p,T0,x0,y,FugEnum.DEW_TEMPERATURE,condensables,false)
    return T,vl0,vv0,x
end

function improve_bubbledew_suggestion(model::PTFlashWrapper,p0,T0,x,y,method,in_media,high_conditions)
    TT = Base.promote_eltype(model,p0,T0,x,y)
    p,T = TT(p0),TT(T0)
    vl = volume(model,p,T,x,phase = :l)/sum(x)
    vv = volume(model,p,T,y,phase = :v)/sum(y)
    return p,T,x,y,vl,vv
end