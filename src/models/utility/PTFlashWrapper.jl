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

function modified_lnϕ_pure(wrapper::PTFlashWrapper,p,T,i;phase = :unknown)
    ps,vl,_ = wrapper.sat[i]
    lnϕsat = wrapper.fug[i]
    isnan(ps) && update_temperature!(wrapper,T)
    isnan(ps) && return ps
    new_phase = if is_unknown(phase)
        p > ps ? :l : :v
    else
        phase
    end
    if is_vapour(phase)
        RT = Rgas(wrapper)*T
        gasmodel = gas_model(wrapper.pures[i])
        vv = volume(gasmodel,p,T,phase = :v)
        lnϕv = VT_lnϕ_pure(gasmodel,vv,T,p)
        Δd = log(ps/p)
        (gas_model isa IdealModel) || (Δd += vl*(p - ps)/RT + lnϕsat)
        return lnϕv - Δd
    else
        return zero(Base.promote_eltype(wrapper,p,T))
    end

end

function tp_flash_fast_K0!(K,wrapper::PTFlashWrapper,p,T,z)
    K .= first.(wrapper.sat) ./ p
    return true
end

function suggest_K!(K,wrapper::PTFlashWrapper,p,T,z,cache = nothing,pure = wrapper.pures)
    phase = identify_phase(wrapper,p,T,z)
    sat = wrapper.sat
    lnϕsat = wrapper.fug
    lnϕz,v = modified_lnϕ(wrapper,p,T,z,cache,phase = phase)
    log∑z = log(sum(z))
    RT = Rgas(wrapper)*T
    for i in 1:length(z)
        ps,vl,_ = sat[i]
        di = lnϕz[i] + log(z[i]) -  log∑z
        lnϕv = modified_lnϕ_pure(wrapper,p,T,i,phase = :v)
        lnϕl = modified_lnϕ_pure(wrapper,p,T,i,phase = :l)
        tpd_v = lnϕv - di
        tpd_l = lnϕl - di
        if tpd_l < 0 && tpd_v < 0
            K[i] = exp(lnϕl)/exp(lnϕv)
        elseif tpd_l < 0 && tpd_v >= 0
            K[i] = exp(lnϕl)/exp(lnϕz[i])
        elseif tpd_l >= 0 && tpd_v < 0
            K[i] = exp(lnϕz[i])/exp(lnϕv)
        else #=tpd_l >= 0 && tpd_v >= 0=#
            K[i] = ps/p
        end
    end
    return K
end

function K0_lle_init(wrapper::PTFlashWrapper,p,T,z,cache = tpd_cache(wrapper,p,T,z);reduced = true)
    return K0_lle_init(__γ_unwrap(wrapper),p,T,z,cache;reduced)
end

function update_volume!(model::PTFlashWrapper,result,p = pressure(result),T = temperature(result))
    for i in 1:numphases(result)
        phase_i = identify_phase(result,i)
        if is_unknown(phase_i) || is_liquid(phase_i)
            result.volumes[i] = volume(model,p,T,result.compositions[i],phase = phase_i)
        end
    end
    return nothing
end

for xy in [:ph,:ps,:ts,:vt]
    xyz = Symbol(xy,:_flash)
    @eval begin 
        function init_preferred_method(method::typeof($xyz),model::PTFlashWrapper,kwargs)
            return RRXYFlash(;kwargs...)
        end
    end
end

for xy in [:qt, :qp]
    xyz = Symbol(xy,:_flash)
    @eval begin 
        function init_preferred_method(method::typeof($xyz),model::PTFlashWrapper,kwargs)
            return RRQXFlash(;kwargs...)
        end
    end
end
