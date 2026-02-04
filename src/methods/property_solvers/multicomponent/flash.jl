"""
    FlashResult(compositions,fractions,volumes,data::FlashData)
    FlashResult(model,p,T,z,compositions,fractions,volumes,g = nothing;sort = true)
    FlashResult(model,p,T,compositions,fractions,volumes,g = nothing;sort = true)
    FlashResult(model,p,T,z;phase = :unknown)
    FlashResult(p,T,z,compositions,fractions,volumes,g = nothing;sort = true)
    FlashResult(p,T,compositions,fractions,volumes,g = nothing;sort = true)
    FlashResult(flash::FlashResult,g = nothing;sort = true)

Structure used to contain the result of a flash.
Contains a list of molar compositions, a list of molar amounts per phase, a list of molar volumes and an auxiliary struct, `FlashData`, containing the pressure, temperature and reduced Gibbs energy.
when an `EoSModel` is used as an input for a `FlashResult`, the reduced molar Gibbs energy (g = g/NRT) is calculated, if not provided.
By default, the phases are sorted by volume, this can be changed by passing the keyword argument `sort = false`
`FlashResult(model,p,T,z;phase)` constructs a single phase `FlashResult`.
If the bulk composition `z` is provided, it will be used to scale the fractions, forcing `sum(fractions) == sum(z)`
"""
struct FlashResult{C,B,V,D}
    compositions::C
    fractions::B
    volumes::V
    data::D
end

"""
    FlashData

Auxiliary struct that contains information about the current `FlashResult` object. It stores the pressure, temperature and reduced Gibbs energy (`g = G/nRT`)
"""
struct FlashData{R}
    p::R
    T::R
    g::R
end

Base.show(io::IO,::MIME"text/plain",options::FlashData) = show_as_namedtuple(io,options)
Base.show(io::IO,options::FlashData) = show_as_namedtuple(io,options)

Solvers.primalval(data::FlashData) = FlashData(primalval(data.p),primalval(data.T),primalval(data.g))
Solvers.primalval(result::FlashResult) = FlashResult(primalval.(result.compositions),primalval(result.fractions),primalval(result.volumes),primalval(result.data))

function Solvers.recursive_fd_extract_derivative(X::XX,result::FlashResult) where XX
    comps = Solvers.recursive_fd_extract_derivative.(X,result.compositions)
    β = Solvers.recursive_fd_extract_derivative(X,result.fractions)
    vols = Solvers.recursive_fd_extract_derivative(X,result.volumes)
    T = Solvers.recursive_fd_extract_derivative(X,result.data.T)
    p = Solvers.recursive_fd_extract_derivative(X,result.data.p)
    g = Solvers.recursive_fd_extract_derivative(X,result.data.g)
    data = FlashData(p,T,g)
    return FlashResult(comps,β,vols,data)
end

function Solvers.recursive_fd_value(result::FlashResult)
    comps = Solvers.recursive_fd_value.(result.compositions)
    β = Solvers.recursive_fd_value(result.fractions)
    vols = Solvers.recursive_fd_value(result.volumes)
    T = Solvers.recursive_fd_value(result.data.T)
    p = Solvers.recursive_fd_value(result.data.p)
    g = Solvers.recursive_fd_value(result.data.g)
    data = FlashData(p,T,g)
    return FlashResult(comps,β,vols,data)
end

function Base.copyto!(dest::FlashResult,src::FlashResult)
    @assert numphases(dest) == numphases(src)
    @assert length(dest.compositions[1]) == length(src.compositions[1])
    copyto!(dest.volumes,src.volumes)
    copyto!(dest.fractions,src.fractions)
    dest_comps,src_comps = dest.compositions,src.compositions
    for i in 1:numphases(dest)
        copyto!(dest_comps[i],src_comps[i])
    end
    return dest
end

function FlashData(p::R1,T::R2,g::R3) where{R1,R2,R3}
    if g === nothing
        FlashData(promote(p,T)...)
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

#constructor that fills the Gibbs energy automatically
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
    comps = [z ./ ∑z]
    volumes = [volume(model,p,T,z;phase = phase)/∑z]
    β = [∑z*one(eltype(volumes))]
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

function FlashResultInvalid(nc::SVector{N,T},val::Number) where {N,T}
    nan = zero(T)/zero(T)
    xx = nc .* nan
    comps = [xx]
    volumes = [nan]
    β = [nan]
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


Base.iterate(result::FlashResult) = (result.compositions, Val(:β))
Base.iterate(result::FlashResult, ::Val{:β}) = (result.fractions, Val(:v))
Base.iterate(result::FlashResult, ::Val{:v}) = (result.volumes, Val(:data))
Base.iterate(result::FlashResult, ::Val{:data}) = (result.data, Val(:done))
Base.iterate(result::FlashResult, ::Val{:done}) = nothing

Base.eltype(result::FlashResult) = Base.promote_eltype(result.compositions[1],result.volumes,result.fractions,result.data.T)

function index_expansion(x::FlashResult,idr::AbstractVector)
    if length(idr) == length(x.compositions[1])
        return x
    end
    newcomps = map(Base.Fix2(index_expansion,idr),x.compositions)
    return FlashResult(newcomps,x.fractions,x.volumes,x.data)
end

"""
    eachphase(result::FlashResult)

Iterates over the values of (V,T,z,β) for each phase.
"""
function eachphase(x::FlashResult)
    return Iterators.zip(x.volumes,FillArrays.fill(x.data.T,numphases(x)),x.compositions,x.fractions)
end

temperature(state::FlashResult) = state.data.T
pressure(state::FlashResult) = state.data.p
volume(state::FlashResult) = dot(state.fractions,state.volumes)
molar_density(state::FlashResult) = sum(state.fractions)/volume(state)

pressure(model::EoSModel,state::FlashResult) = pressure(state)
temperature(model::EoSModel,state::FlashResult) = temperature(state)
volume(model::EoSModel,state::FlashResult) = volume(state)
molar_density(model::EoSModel,state::FlashResult) = molar_density(state)

function __molecular_weight(model,state::FlashResult)
    comps, β, volumes, data = state
    ∑mi = zero(eltype(comps[1]))
    for i in 1:length(comps)
        mwi = molecular_weight(model,comps[i])
        ∑mi += β[i]*mwi
    end
    return ∑mi
end

function mass_density(model::EoSModel,state::FlashResult)
    V = volume(model,state)
    molar_weight = molecular_weight(model,state)
    return molar_weight/V
end

function mass_density(model::EoSModel,state::FlashResult, i::Integer)
    vi,T,xi,βi = state.volumes[i],state.data.T,state.compositions[i],state.fractions[i]
    molar_weight = molecular_weight(model,xi)
    return molar_weight/vi
end

function volume(model::EoSModel,state::FlashResult, i::Integer)
    vi,T,xi,βi = state.volumes[i],state.data.T,state.compositions[i],state.fractions[i]
    return vi*βi
end

function gibbs_free_energy(model::EoSModel,state::FlashResult)
    p = pressure(state)
    res = zero(Base.promote_eltype(model,state))
    for (vi,T,xi,βi) in eachphase(state)
        res += βi*VT_gibbs_energy(model,vi,T,xi,p)
    end
    return res
end

function gibbs_free_energy(model::EoSModel,state::FlashResult, i)
    p = pressure(state)
    vi,T,xi,βi = state.volumes[i],state.data.T,state.compositions[i],state.fractions[i]
    return βi*VT_gibbs_energy(model,vi,T,xi,p)
end

for prop in [:enthalpy,:entropy,:internal_energy,:helmholtz_free_energy]
    @eval begin
            function $prop(model::EoSModel,state::FlashResult)
                res = zero(Base.promote_eltype(model,state))
                for (vi,T,xi,βi) in eachphase(state)
                    res += βi*VT0.$prop(model,vi,T,xi)
                end
                return res
            end

            function $prop(model::EoSModel,state::FlashResult, i::Integer)
                res = zero(Base.promote_eltype(model,state))
                vi,T,xi,βi = state.volumes[i],state.data.T,state.compositions[i],state.fractions[i]
                res += βi*VT0.$prop(model,vi,T,xi)
                return res
            end
        end
end

mass_entropy(model::EoSModel,state::FlashResult) = entropy(model,state)/molecular_weight(model,state)
mass_enthalpy(model::EoSModel,state::FlashResult) = mass_enthalpy(model,state)/molecular_weight(model,state)
mass_internal_energy(model::EoSModel,state::FlashResult) = mass_internal_energy(model,state)/molecular_weight(model,state)
mass_gibbs_free_energy(model::EoSModel,state::FlashResult) = mass_gibbs_free_energy(model,state)/molecular_weight(model,state)
mass_helmholtz_free_energy(model::EoSModel,state::FlashResult) = mass_helmholtz_free_energy(model,state)/molecular_weight(model,state)

for prop in [:mass_enthalpy,:mass_entropy,:mass_internal_energy,:mass_helmholtz_free_energy,:mass_gibbs_free_energy]
    @eval begin
            function $prop(model::EoSModel,state::FlashResult, i::Integer)
                res = zero(Base.promote_eltype(model,state))
                vi,T,xi = state.volumes[i],state.data.T,state.compositions[i]
                res += VT0.$prop(model,vi,T,xi)
                return res
            end
        end
end

function assert_only_phase_index(state::FlashResult)
    np = numphases(state)
    if isone(np)
        return 1
    elseif np > 1 #on some systems, there could be multiple phases, but only one fraction is nonzero
        βmax,imax = findmax(state.fractions)
        isfinite(βmax) || return imax #non finite value, return NaN, it will fail anyway.
        ∑β = sum(state.fractions)
        if βmax ≈ ∑β && all(>=(0),state.fractions)
            return imax
        else
            return 0
        end
    end
end

for prop in [:isochoric_heat_capacity, :isobaric_heat_capacity, :adiabatic_index,
    :mass_isochoric_heat_capacity, :mass_isobaric_heat_capacity,
    :isothermal_compressibility, :isentropic_compressibility, :speed_of_sound,
    :isobaric_expansivity, :joule_thomson_coefficient, :inversion_temperature,
    #higher :derivative :order :properties
    :fundamental_derivative_of_gas_dynamics,
    #volume :properties
    :compressibility_factor,:identify_phase]
    @eval begin
        function $prop(model::EoSModel,state::FlashResult)
            i = assert_only_phase_index(state::FlashResult)
            T = temperature(state)
            p = pressure(state)
            if iszero(i)
                invalid_property_multiphase_error($prop,numphases(state),p,T)
            end

            x,v = state.compositions[i],state.volumes[i]
            return VT0.$prop(model,v,T,x)
        end

        function $prop(model::EoSModel,state::FlashResult, i::Int)
            x,v = state.compositions[i],state.volumes[i]
            return VT0.$prop(model,v,T,x)
        end
    end
end

function _multiphase_gibbs(model,result,vapour_phase_index = 0)
    gibbs_energy(model,result)/Rgas(model)/result.data.T
end

function _multiphase_gibbs(model::PTFlashWrapper,result,vapour_phase_index = 0)
    modified_gibbs(model,result;vapour_phase_index)/Rgas(model)/result.data.T
end

function __mpflash_phase(vapour_phase_index,i) 
    if vapour_phase_index != 0
        phase = vapour_phase_index == i ? :vapour : :liquid
    else
        phase = :unknown
    end
    return phase
end

function modified_gibbs(model,result::FlashResult;vapour_phase_index = 0)
    np = numphases(result)
    g = zero(Base.promote_eltype(result.compositions[1],result.fractions,result.volumes,result.data.p,result.data.T,model))
    p,T = result.data.p,result.data.T
    v = result.volumes
    β = result.fractions
    x = result.compositions
    for i in 1:np
        phase = __mpflash_phase(vapour_phase_index,i)
        gi,_ = modified_gibbs(model,p,T,x[i],phase,v[i])
        g += β[i]*gi
    end
    return g
end

#utilities to add/remove phases from an existing FlashResult

function split_phase!(result::FlashResult,i::Integer,wj,βj,vj)
    β = result.fractions
    comps = result.compositions
    volumes = result.volumes
    data = result.data
    p,T = data.p,data.T
    n = sum(β)

    wi0,vi0,βi0 = comps[i],volumes[i],β[i]
    _wj = copy(wi0)
    nj = sum(wj)
    _wj .= wj ./ nj
    βj = Solvers.positive_linesearch(wi0,_wj,βj,s = -1,decay = 0.95)
    βi = (1 - βj/n) * βi0
    β[i] = βi0 - βi0*βj
    #add new phase

    push!(β,βi0*βj)
    push!(volumes,vj/nj)
    push!(comps,_wj)

    #remove moles from phase i
    wi0 .= wi0 .-  βj .* _wj
    wi0 ./= sum(wi0)

    return result
end

#merges phases i and j, leaving i and removing j
function merge_phase!(result::FlashResult,i,j)
    comps,β,volumes = result.compositions,result.fractions,result.volumes
    wi,vi,βi = comps[i],volumes[i],β[i]
    wj,vj,βj = comps[j],volumes[j],β[j]
    βk = βi + βj
    wk = wi
    wk .= (wi .* βi .+ wj .* βj) ./ βk
    vk = (vi * βi + vj * βj) / βk
    β[i] = βk
    volumes[i] = vk

    #delete phase j
    return delete_phase!(result,j)
end

function delete_phase!(result::FlashResult,i)
    comps,β,volumes = result.compositions,result.fractions,result.volumes
    Base.deleteat!(comps,i)
    Base.deleteat!(volumes,i)
    Base.deleteat!(β,i)
    return result
end

function findfirst_duplicate_phases(comps,β,volumes,ignore_zeros = true)
    equal_phases = (0,0)
    for i in 1:length(comps)
        xi,vi,βi = comps[i],volumes[i],β[i]
        iszero(βi) && ignore_zeros && continue
        for j in (i+1):length(comps)
            xj,vj,βj = comps[j],volumes[j],β[j]
            iszero(βj) && ignore_zeros && continue
            #equality criteria used in the HELD algorithm
            if isnan(vi) && isnan(vj)
                equal_v = true
            else
                equal_v = abs(1/vi - 1/vj) <= 1e-5
            end
            if dnorm(xi,xj,Inf) <= 1e-5 && equal_v
                return minmax(i,j)
            end
        end
    end
    return (0,0)
end

function findfirst_duplicate_phases(result::FlashResult,ignore_zeros = true)
    comps,β,volumes = result.compositions,result.fractions,result.volumes
    return findfirst_duplicate_phases(comps,β,volumes,ignore_zeros)
end

function merge_duplicate_phases!(result::FlashResult;ignore_zeros = true)
    nc = numphases(result)
    for i in 1:(nc*nc)
        i,j = findfirst_duplicate_phases(result,ignore_zeros)
        if i == j == 0
            break
        end
        merge_phase!(result,i,j)
    end
    return result
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

is_vle(method::FlashMethod) = is_vle(method.equilibrium)
is_lle(method::FlashMethod) = is_lle(method.equilibrium)
is_unknown(method::FlashMethod) = is_unknown(method.equilibrium)

@noinline function incorrect_np_flash_error(method,result)
    np = numphases(result)
    s = np == 1 ? "" : "s"
    throw(ArgumentError("$method does not support an input with $np phase$s as an initial point. Got the following input: \n\n $result"))
end

include("flash/general_flash.jl")

function xy_flash_ad(result,tup,tup_primal,spec1,spec2)
    if any(has_dual,tup)
        np = numphases(result)
        if isone(np) || isone(np - count(iszero,result.fractions))
            return xy_flash_ad1(result,tup,tup_primal,spec1,spec2)
        end

        function f(input,tups)
            model0,_val1,_val2,zbulk = tups
            TT = Base.promote_eltype(model0,_val1,_val2,zbulk,input)   
            output = similar(input,TT)
            spec = FlashSpecifications(spec1,_val1,spec2,_val2)
            xy_flash_neq(output,model0,zbulk,np,input,spec,nothing)
            return output
        end
        model,val1,val2,z = tup
        nc = length(model)
        λx = vcat(reduce(vcat,result.compositions),result.volumes,result.fractions,result.data.T)
        ∂spec = FlashSpecifications(spec1,val1,spec2,val2)
        ∂x = __gradients_for_root_finders(λx,tup,tup_primal,f)
        ∂comps,∂β,∂volumes,∂T = xy_input_to_result(∂spec,∂x,np,nc,z)

        if spec1 == pressure
            ∂p = oftype(∂T,val1)
        elseif spec2 == pressure
            ∂p = oftype(∂T,val2)
        else
            ∂p = pressure(model,∂volumes[end],∂T,∂comps[end])
        end

        if result.data.g isa Number && !isnan(result.data.g)
            return FlashResult(model,∂p,∂T,∂comps,∂β,∂volumes,sort = false)
        end
        return FlashResult(∂comps,∂β,∂volumes,FlashData(∂p,∂T))
    end
    return result
end

function __xy_flash_ad1_fill1(orig::SVector{N,T},val::V) where {N,T,V}
    v = ntuple(Returns(val),Val{N}())
    return SVector{N,V}(v)
end

function __xy_flash_ad1_fill1(orig::AbstractVector{T},val::V) where {T,V}
    dest = similar(orig,V)
    fill!(dest,val)
    return dest
end

function __xy_flash_ad1_fillβ(orig::SVector{N,T},β::B,ix) where {N,T,B}
    v = ntuple(i -> i == ix ? β : zero(β),Val{N}())
    return SVector{N,V}(v)
end

function __xy_flash_ad1_fillβ(orig::AbstractVector{T},β::B,ix) where {T,B}
    dest = similar(orig,B)
    fill!(dest,zero(β))
    dest[ix] = β
    return dest
end

function xy_flash_ad1(result,tup,tup_primal,spec1,spec2)
    
    function f(input,tups)
        model0,_val1,_val2,zbulk = tups
        v0,T0 = input
        f1 = spec_to_vt(model0,v0,T0,zbulk,spec1) - _val1
        f2 = spec_to_vt(model0,v0,T0,zbulk,spec2) - _val2
        return SVector(f1,f2)
    end
    i = findfirst(!iszero,result.fractions)
    λT = result.data.T
    λv = result.volumes[i]
    λx = SVector(λv,λT)
    ∂x = __gradients_for_root_finders(λx,tup,tup_primal,f)
    ∂v,∂T = ∂x[1],∂x[2]
    model,val1,val2,z = tup
    ∂β1 = sum(z)
    ∂comp1 = z ./ ∂β1
    ∂β = __xy_flash_ad1_fillβ(result.fractions,∂β1,i)
    ∂comps = __xy_flash_ad1_fill1(result.compositions,∂comp1)
    ∂volumes = __xy_flash_ad1_fill1(result.volumes,∂v)

    if spec1 == pressure
        ∂p = oftype(∂T,val1)
    elseif spec2 == pressure
        ∂p = oftype(∂T,val2)
    else
        ∂p = pressure(model,∂v,∂T,∂comp1)
    end

    if result.data.g isa Number && !isnan(result.data.g)
        return FlashResult(model,∂p,∂T,∂comps,∂β,∂volumes,sort = false)
    end
    
    return FlashResult(∂comps,∂β,∂volumes,FlashData(∂p,∂T))
end

include("flash/PT.jl")
include("flash/PH.jl")
include("flash/PS.jl")
include("flash/VT.jl")
include("flash/TS.jl")
include("flash/QT.jl")
include("flash/QP.jl")
include("flash/flash_HSU.jl")


export FlashResult, FlashData
