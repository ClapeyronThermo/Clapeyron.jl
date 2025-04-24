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

Solvers.primalval(data::FlashData) = FlashData(primalval(data.p),primalval(data.T),primalval(data.g))
Solvers.primalval(result::FlashResult) = FlashResult(primalval.(result.compositions),primalval(result.fractions),primalval(result.volumes),primalval(result.data))

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
    volumes = [volume(model,p,T,z;phase = phase)/∑z]
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
    p = pressure(state)
    res = zero(Base.promote_eltype(model,state))
    for (vi,T,xi,βi) in eachphase(state)
        res += βi*VT_gibbs_energy(model,vi,T,xi,p)
    end
    return res
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

function assert_only_phase_index(state::FlashResult)
    if isone(numphases(prop))
        return 1
    elseif numphases(prop) > 1 #on some systems, there could be multiple phases, but only one fraction is nonzero
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

@noinline function __multiphase_onephase_function_error(f,np,p,T)
    throw(ArgumentError("The state at p = $p, T = $T has $np phases, it cannot be used to evaluate $f"))
end



for prop in [:isochoric_heat_capacity, :isobaric_heat_capacity, :adiabatic_index,
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
                    __multiphase_onephase_function_error($prop,np,p,T)
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

function _multiphase_gibbs(model,p,T,result)
    if result isa FlashResult
        return _multiphase_gibbs(model,p,T,(result.compositions,result.fractions,result.volumes))
    end
    if model isa PTFlashWrapper && length(β) == 2 #TODO: in the future, PTFlashWrappers could be multiphase
        if model.model isa RestrictedEquilibriaModel
            return zero(eltype(β))
        end
    end
    comps = result[1]
    β = result[2]
    volumes = result[3]
    data = FlashResult(model,p,T,comps,β,volumes)
    return Rgas(model)*T*data.data.g
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
    @show β
    @show sum(β)
    #remove moles from phase i
    wi0 .= wi0 .-  βj .* _wj
    @show wi0
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

function findfirst_duplicate_phases(comps,β,volumes)
    equal_phases = (0,0)
    for i in 1:length(comps)
        xi,vi,βi = comps[i],volumes[i],β[i]
        iszero(βi) && continue
        for j in (i+1):length(comps)
            xj,vj,βj = comps[j],volumes[j],β[j]
            iszero(βj) && continue
            #equality criteria used in the HELD algorithm
            if dnorm(xi,xj,Inf) <= 1e-5 && abs(1/vi - 1/vj) <= 1e-5
                return minmax(i,j)
            end
        end
    end
    return (0,0)
end

function findfirst_duplicate_phases(result::FlashResult)
    comps,β,volumes = result.compositions,result.fractions,result.volumes
    return findfirst_duplicate_phases(comps,β,volumes)
end

function merge_duplicate_phases!(result::FlashResult)
    nc = numphases(result)
    for i in 1:(nc*nc)
        i,j = findfirst_duplicate_phases(result)
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

function tp_flash_1phase(model,p,T,z,result_primal)
    if has_dual(model) || has_dual(p) || has_dual(T) || has_dual(z)
        n1 = similar(z,Base.promote_eltype(model,p,T,z))
        β = similar(n1,1)
        v1 = similar(β)
        ∑z = sum(z)
        β[1] = ∑z
        n1 .= z
        n1 ./= ∑z
        v1[1] = volume_ad(model,result_primal.volumes[1],T,n1,p)
        nx = [n1]
        result = FlashResult(nx,β,v1,FlashData(p,T))
        if data_primal.g isa Number && !isnan(data_primal.g)
            return FlashResult(model,p,T,nx,β,vols,sort = false)
        end
        return FlashResult(nx,β,vols,FlashData(p,T))
    else
        return result_primal
    end
end

function tp_flash_ad(model,p,T,z,result_primal)
    if has_dual(model) || has_dual(p) || has_dual(T) || has_dual(z)
        data_primal = result_primal.data
        T_primal,p_primal = data_primal.T,data_primal.p
        v_primal = result_primal.volumes
        np,nc = numphases(result_primal),length(model)

        if np == 1
            return tp_flash_1phase(model,p,T,z,result_primal)
        end

        comps_primal = result_primal.compositions
        β_primal = result_primal.fractions
        βmax,jmax = findmax(β_primal) #find biggest fraction, use that as an anchor to update the rest of compositions
        v1_primal = v_primal[jmax]
        w1 = comps_primal[jmax]
        ∂ϕ1 = ∂lnϕ_cache(model, p_primal, T_primal, w1, Val{true}())
        ∂ϕj = ∂lnϕ_cache(model, p_primal, T_primal, w1, Val{true}())
        lnϕ1, ∂lnϕ∂n1, ∂lnϕ∂P1, ∂lnϕ∂T1, _ = ∂lnϕ∂n∂P∂T(model,p_primal,T_primal,w1,∂ϕ1,vol0 = v1_primal)
        ∂lnϕ∂n1 .-= 1/βmax
        for i in 1:nc
            ∂lnϕ∂n1[i,i] += 1/(w1[i]*βmax)
        end
        Hcache = similar(∂lnϕ∂n1)
        piv = zeros(Int,nc)
        s = similar(lnϕ1)
       
        n1 = similar(w1,Base.promote_eltype(model,p,T,z))
        vols = similar(n1)
        β = similar(n1)
        n1 .= primalval.(z)
        nx = fill(n1,0)
        for j in 1:np
            wj = comps_primal[j]
            vj_primal = v_primal[j]
            βj = β_primal[j]
            if j != jmax
                nj = similar(n1)
                nj .= wj .* βj
                lnϕj, ∂lnϕ∂nj, ∂lnϕ∂Pj, ∂lnϕ∂Tj, _ = ∂lnϕ∂n∂P∂T(model,p_primal,T_primal,wj,∂ϕj,vol0 = vj_primal)
                ∂lnϕ∂nj .-= 1/βj
                for i in 1:nc
                    ∂lnϕ∂nj[i,i] += 1/(wj[i]*βj)
                end

                H = ∂lnϕ∂nj
                H .+= ∂lnϕ∂n1
                Hcache .= H
                
                #∂ϕxj = eye./nj .- 1/βj.+ ∂lnϕ∂nxj/βj
                #∂ϕx1 = eye./n1 .- 1/β1 .+ ∂lnϕ∂ny/β1
                #H .= ∂ϕxj .+ ∂ϕx1

                if has_dual(p)
                    Hcache .= H
                    lu = Solvers.unsafe_LU!(Hcache,piv)
                    s .= ∂lnϕ∂Pj .- ∂lnϕ∂P1
                    ldiv!(lu,s)
                    nj .-= (p - p_primal) .* s
                end

                if has_dual(T)
                    Hcache .= H
                    lu = Solvers.unsafe_LU!(Hcache,piv)
                    s .= ∂lnϕ∂Tj .- ∂lnϕ∂T1
                    ldiv!(lu,s)
                    dwjdT = H\(∂lnϕ∂Tj .- ∂lnϕ∂T1)
                    nj .-= (T - T_primal) .* s
                end
                n1 .-= nj
                nj ./= sum(nj)
                vols[j] = volume_ad(model,vj_primal,T,nj,p)
                push!(nx,nj)
            else
                push!(nx,n1)
            end
        end
        n1./= sum(n1)

        vols[jmax] = volume_ad(model,v_primal[jmax],T,n1,p)
        for j in 1:np
            wj = comps_primal[j]
            vj_primal = v_primal[j]
        end

        result = FlashResult(nx,β,vols,FlashData(p,T))
        if data_primal.g isa Number && !isnan(data_primal.g)
            return FlashResult(model,p,T,nx,β,vols,sort = false)
        end
        return FlashResult(nx,β,vols,FlashData(p,T))
        
    else
        return result_primal
    end
end

include("flash/general_flash.jl")
include("flash/PT.jl")
include("flash/PH.jl")
include("flash/PS.jl")
include("flash/VT.jl")
include("flash/TS.jl")
include("flash/QT.jl")
include("flash/QP.jl")

export FlashResult, FlashData
