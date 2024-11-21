"""
    Clapeyron.FlashResult

Structure used to contain the result of a flash.

"""
struct FlashResult{C,B,D}
    components::C
    fractions::B
    volumes::B
    data::D
end

struct FlashData{R}
    p::R
    T::R
    dG::R
end

function FlashData(p::R1,T::R2,dG::R2) where{R1,R2,R3}
    return FlashData(promote(p,T,dG)...)
end

FlashData(p,T) = FlashData(p,T,1.0*zero(p))

function Base.show(io::IO,mime::MIME"text/plain",obj::FlashResult)
    comps,β,volumes,data = obj
    np = length(comps)
    compact_io = IOContext(io, :compact => true)
    print(io,"Flash result at T = ")
    print(compact_io,data.T)
    print(io,", p = ")
    print(compact_io,data.p)
    print(io," with $np phase")
    if np > 1
        print(io,"s")
    end
    println(io,":")
    nt = map(comps,β,volumes) do xi,βi,vi
        (x = xi,β = βi,v = vi)
    end
    Base.print_matrix(IOContext(io, :compact => true),nt)
end

Base.getindex(x::Flash,i::Int) = getfield(x,i)
function Base.iterate(x::FlashResult)
    return (x[1],2)
end

function Base.iterate(x::FlashResult,state)
    if state > 4
        return nothing
    else
        return (x[state],state + 1)
    end
end

temperature(model::EoSModel,state::FlashResult) = state.T
pressure(model::EoSModel,state::FlashResult) = state.p
function volume(model::EoSModel,state::FlashResult)
    comps, β, volumes, data = state
    return dot(β,volumes)
end

function molar_density(model::EoSModel,state::FlashResult)
    comps, β, volumes, data = state
    v = volume(model,state)
    n = sum(β)
    return n/v
end

function molecular_weight(model::EoSModel,state::FlashResult)
    comps, β, volumes, data = state
    ∑mi = zero(eltype(comps))
    mw = mw(model)
    for i in 1:length(comps)
        mwi = comp_molecular_weight(mw,comps[i])
        ∑mi = β[i]*mwi
    end
    return ∑mi
end

function mass_density(model::EoSModel,state::FlashResult)
    V = volume(model,state)
    molar_weight = molecular_weight(model,state)
    return molar_weight/V
end

function gibbs_free_energy(model::EoSModel,state::FlashResult)
    comps, β, volumes, data = state
    T = data.T
    p = data.p
    res = zero(Base.promote_eltype(comps[1],volumes[1],T,model))
    for i in 1:length(comps)
        xi,βi,vi = comps[i],β[i],volumes[i]
        res += βi*VT_gibbs_energy(model,vi,T,xi,p)
    end
    return res
end


for prop in [:enthalpy,:entropy,:internal_energy,:helmholtz_free_energy]
    @eval begin
            function $prop(model::EoSModel,state::FlashResult)
                comps, β, volumes, data = state
                T = data.T
                res = zero(Base.promote_eltype(comps[1],volumes[1],T,model))
                for i in 1:length(comps)
                    xi,βi,vi = comps[i],β[i],volumes[i]
                    res += βi*VT.$prop(model,vi,T,xi)
                end
                return res
            end
        end
    end


function _multiphase_gibbs(model,p,T,result)
    if model isa PTFlashWrapper && length(β) == 2 #TODO: in the future, PTFlashWrappers could be multiphase
        if model.model isa RestrictedEquilibriaModel
            return zero(eltype(β))
        end
    end
    comps = result[1]
    β = result[2]
    volumes = result[3]
    data = FlashResult(model,p,T,comps,β,volumes)
    return Rgas(model)*T*data.data.dG
end

function FlashResult(model::EoSModel,p,T,comps,β,volumes)
    data = FlashData(p,T)
    flash = FlashResult(comps,β,volumes,data)
    g_mix = gibbs_free_energy(model,flash)
    dg = g_mix/(Rgas(model)*T)
    newdata = PTFlashData(p,T,dg)
    return FlashResult(comps,β,volumes,newdata)
end


"""
    FlashMethod <: ThermodynamicMethod

Abstract type for flash routines. 

to add a new method, it is necessary to define the following functions, depending on the type of supported flash:

- P-T Flash: `tp_flash_impl(model,p,T,n,method)`
- P-H Flash: `ph_flash_impl(model,p,h,n,method)`
- P-S Flash: `ps_flash_impl(model,p,s,n,method)`

If the flash method supports more than 2 phases, then it requires defining `numphases(method)` 
If the method accept component-dependent inputs, it should also define `index_reduction(method,nonzero_indices)`
"""
abstract type FlashMethod <: ThermodynamicMethod


"""
    numphases(method::FlashMethod)

Return the number of phases supported by a flash method. By default it is set to 2.
If the method allows it, you can set the number of phases by doing `method(;numphases = n)`.
"""
numphases(method::FlashMethod) = 2

"""
    supports_reduction(method::FlashMethod)::Bool

Checks if a Flash method supports index reduction (the ability to prune model components with compositions lower than a threshold). 
All current Clapeyron.jl methods support index reduction, but some methods that alllow passing a cache could have problems.
"""
supports_reduction(method::FlashMethod) = true