"""
    FlashResult(compositions,fractions,volumes,data::FlashData)
    FlashResult(model,p,T,z,compositions,fractions,volumes,g = nothing;sort = true)
    FlashResult(model,p,T,compositions,fractions,volumes,g = nothing;sort = true)
    FlashResult(model,p,T,z;phase = :unknown)
    FlashResult(p,T,z,compositions,fractions,volumes,g = nothing;sort = true)
    FlashResult(p,T,compositions,fractions,volumes,g = nothing;sort = true)
    FlashResult(flash::FlashResult,g = nothing;sort = true)

Structure used to contain the result of a flash.
Contains a list of molar compositions, a list of molar amounts per phase, a list of molar volumes and an auxiliary struct, `FlashData`, containing the pressure, temperature and reduced gibbs energy.
when an `EoSModel` is used as an input for a `FlashResult`, the reduced molar gibbs energy (g = g/NRT) is calculated, if not provided.
By default, the phases are sorted by volume, this can be changed by passing the keyword argument `sort = false`
`FlashResult(model,p,T,z;phase)` constructs a single phase `FlashResult`.
If the bulk composition `z` is provided, it will be used to scale the fractions, forcing `sum(fractions) == sum(z)`
"""
struct FlashResult{C,B,D}
    compositions::C
    fractions::B
    volumes::B
    data::D
end

"""
    FlashData

Auxiliary struct that contains information about the current `FlashResult` object. It stores the pressure, temperature and reduced gibbs energy (`g = G/nRT`)
"""
struct FlashData{R}
    p::R
    T::R
    g::R
end

Base.show(io::IO,::MIME"text/plain",options::FlashData) = show_as_namedtuple(io,options)
Base.show(io::IO,options::FlashData) = show_as_namedtuple(io,options)

function FlashData(p::R1,T::R2,g::R3) where{R1,R2,R3}
    if g === nothing
        FlashData(p,T)
    else
        return FlashData(promote(p,T,g)...)
    end
end

FlashData(p,T) = FlashData(p,T,1.0*zero(p))

#mol checker, with gibbs
function FlashResult(model::EoSModel,p,T,z::Union{Number,AbstractVector{<:Number}},comps,β,volumes,g = nothing;sort = true)
    ∑β = sum(β)
    ∑z = sum(z)
    if !isapprox(∑z,∑β)
        _β = β * (∑z/∑β)
    else
        _β = β * one(∑z/∑β)
    end
    return FlashResult(model,p,T,comps,_β,volumes,g;sort)
end

#constructor that fills the gibbs energy automatically
function FlashResult(model::EoSModel,p,T,comps,β,volumes,g = nothing;sort = true)
    if g == nothing
        flash = FlashResult(p,T,comps,β,volumes,sort = false)
        Gmix = gibbs_free_energy(model,flash)
        _g = Gmix/(sum(β)*Rgas(model)*T)
    else
        _g = g
    end
    return FlashResult(p,T,comps,β,volumes,_g;sort)
end

#mol checker, without gibbs
function FlashResult(p::Number,T::Number,z::Union{Number,AbstractVector{<:Number}},comps,β,volumes,g = nothing;sort = true)
    ∑β = sum(β)
    ∑z = sum(z)
    if !isapprox(∑z,∑β)
        _β = β * (∑z/∑β)
    else
        _β = β * one(∑z/∑β)
    end
    return FlashResult(p,T,comps,_β,volumes,g;sort)
end

#constructor that does not fill the gibbs field
function FlashResult(p::Number,T::Number,comps,β,volumes,g = nothing;sort = true)
    data = FlashData(p,T,g)
    if !sort || issorted(volumes)
        return FlashResult(comps,β,volumes,data)
    else
        idx = sortperm(volumes)
        return FlashResult(comps[idx],β[idx],volumes[idx],data)
    end
end

#flash remaker
function FlashResult(x::FlashResult,g = nothing;sort = true)
    comps,β,volumes,data = x.compositions,x.fractions,x.volumes,x.data
    if g !== nothing
        _g = g
    else
        _g = data.g
    end
    return FlashResult(data.p,data.T,comps,β,volumes,_g;sort)
end

#constructor for single phase
function FlashResult(model::EoSModel,p::Number,T::Number,z;phase = :unknown)
    ∑z = sum(z)
    β = [∑z]
    comps = [z ./ ∑z]
    volumes = [volume(model,p,T,z;phase = phase)]
    return FlashResult(model,p,T,comps,β,volumes;sort = false)
end

#nan constructor
function FlashResultInvalid(nc::Int,val::Number)
    nan = zero(val)/zero(val)
    β = [nan]
    comps = [fill(nan,nc)]
    volumes = [nan]
    data = FlashData(nan,nan,nan)
    return FlashResult(comps,β,volumes,data)
end

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

Base.getindex(x::FlashResult,i::Int) = getfield(x,i)
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

function index_expansion(x::FlashResult,idr::AbstractVector)
    if length(idr) == length(x.compositions[1])
        return x
    end
    newcomps = map(Base.Fix2(index_expansion,idr),x.comps)
    return FlashResult(newcomps,x.fractions,x.volumes,x.data)
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

function __molecular_weight(model,state::FlashResult)
    comps, β, volumes, data = state
    ∑mi = zero(eltype(comps[1]))
    for i in 1:length(comps)
        mwi = molecular_weight(model,comps[i])
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
abstract type FlashMethod <: ThermodynamicMethod end


"""
    numphases(method::FlashMethod)

Return the number of phases supported by a flash method. By default it is set to 2.
If the method allows it, you can set the number of phases by doing `method(;numphases = n)`.
"""
function numphases end

numphases(method::FlashMethod) = 2
numphases(result::FlashResult) = length(result.compositions)
"""
    supports_reduction(method::FlashMethod)::Bool

Checks if a Flash method supports index reduction (the ability to prune model components with compositions lower than a threshold). 
All current Clapeyron.jl methods support index reduction, but some methods that alllow passing a cache could have problems.
"""
supports_reduction(method::FlashMethod) = true

include("flash/general_flash.jl")
include("flash/PT.jl")
include("flash/PH.jl")
include("flash/PS.jl")
