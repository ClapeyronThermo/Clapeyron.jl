struct ∂lnϕTag end

struct ∂lnϕCache{TAG,T,N}
    result::DiffResults.MutableDiffResult{2, T, Tuple{Vector{T}, Matrix{T}}}
    aux::Vector{T}
    lnϕ::Vector{T}
    ∂lnϕ∂n::Matrix{T}
    ∂lnϕ∂P::Vector{T}
    ∂P∂n::Vector{T}
    ∂lnϕ∂T::Vector{T}
    hconfig::ForwardDiff.HessianConfig{TAG,T,N,Vector{ForwardDiff.Dual{TAG, ForwardDiff.Dual{TAG, T, N}, N}},NTuple{2,Vector{ForwardDiff.Dual{TAG, T, N}}}}
    jcache::Vector{T}
    dlnγdT_cache::Vector{ForwardDiff.Dual{TAG,T,1}}
end

function Base.show(io::IO,cache::∂lnϕCache)
    nc = length(cache.lnϕ)
    na = length(cache.aux)
    ∂T = na - nc == 1
    ∂γ∂T = !iszero(length(cache.dlnγdT_cache))
    print(io,typeof(cache))
    print(io,"(nc = ",nc)
    ∂T && print(io,", ∂T = ",∂T)
    ∂γ∂T && print(io,", ∂γ∂T = ",∂γ∂T)
    print(io,")")
end

@inline _cache_tuple(c::∂lnϕCache) = (c.result, c.aux, c.lnϕ, c.∂lnϕ∂n, c.∂lnϕ∂P, c.∂P∂n,
                                      c.∂lnϕ∂T, c.hconfig, c.jcache, c.dlnγdT_cache)
Base.length(::∂lnϕCache) = 10
Base.iterate(c::∂lnϕCache) = iterate(_cache_tuple(c))
Base.iterate(c::∂lnϕCache, state) = iterate(_cache_tuple(c), state)
Base.indexed_iterate(c::∂lnϕCache, i::Int, state = 1) = Base.indexed_iterate(_cache_tuple(c), i, state)
Base.getindex(c::∂lnϕCache,i::Int) = getfield(c,i)

function ∂lnϕ_cache(model::EoSModel, p, T, z, ::Val{B}) where B
    TT = Base.promote_eltype(model,p,T,z)
    nc = length(model)
    lnϕ = zeros(TT,nc)
    aux = zeros(TT,nc + B)
    ∂lnϕ∂n = lnϕ * transpose(lnϕ)
    result = DiffResults.HessianResult(aux)
    ∂lnϕ∂n = lnϕ * transpose(lnϕ)
    ∂lnϕ∂P = similar(lnϕ)
    ∂P∂n = similar(lnϕ)
    hconfig = ForwardDiff.HessianConfig(∂lnϕTag(),result,aux)
    ∂tag = ForwardDiff.tagtype(hconfig.jacobian_config.duals[1][1])
    if has_lnγ_impl(__γ_unwrap(model))
        jcache = similar(aux)
        dlnγdT_cache = Vector{ForwardDiff.Dual{∂tag,TT,1}}(undef,nc)
    else
        jcache = aux
        dlnγdT_cache = Vector{ForwardDiff.Dual{∂tag,TT,1}}(undef,0)
    end
    if B
        ∂lnϕ∂T = similar(lnϕ)
    else
        ∂lnϕ∂T = lnϕ
    end
    return ∂lnϕCache(result,aux,lnϕ,∂lnϕ∂n,∂lnϕ∂P,∂P∂n,∂lnϕ∂T,hconfig,jcache,dlnγdT_cache)
end

function __psi_res(model,T,rhovec)
    nc = length(model)
    _rhovec = @view rhovec[1:nc]
    return psi_res(model,T,_rhovec)
end

function __psi_res(model,data)
    nc = length(model)
    rhovec = @view data[1:nc]
    T = data[nc + 1]
    return psi_res(model,T,rhovec)
end

StaticForwardDiffTags.inner_function(::typeof(__psi_res)) = psi_res

#returns μ_res/RT and p_res/RT from isochoric derivatives of Ψ = A_res/(RT·V), V = 1.
function __μpr(model, T, ρ, cache)
    nc = length(model)

    #deferred forwarddiff tag for faster AD in some cases
    Ψ = @deferred_isochoric_rhovec(__psi_res,psi_res)
    if cache isa ∂lnϕCache
        aux = cache.aux
        ρᵢ = @view aux[1:nc]
        ρᵢ .= ρ
        gconfig = Solvers._GradientConfig(cache.hconfig)
        _result = ForwardDiff.gradient!(cache.result, Ψ, aux, gconfig, Val{false}())
    else
        if cache isa AbstractVector
            result = DiffResults.MutableDiffResult(zero(eltype(cache)), (cache,))
        else
            result = DiffResults.GradientResult(ρ)
        end
        
        _result = ForwardDiff.gradient!(result, Ψ, ρ)
        ρᵢ = ρ
    end

    Ψr = DiffResults.value(_result)
    μ_res = @view DiffResults.gradient(_result)[1:nc]
    p_res = dot(μ_res, ρᵢ) - Ψr
    return μ_res, p_res
end

# Function to compute fugacity coefficient
function lnϕ(model::EoSModel, p, T, z=SA[1.],cache = nothing;
            phase=:unknown,
            vol0=nothing,
            threaded = true,
            vol = volume(model,p,T,z;phase,vol0,threaded))

    ∑z = sum(z)
    RT = Rgas(model)*T
    logZ = log(p*vol/RT/∑z)
    nc = length(z)

    if cache isa Vector
        return lnϕ!(cache, model, p, T, z; vol)
    elseif cache isa ∂lnϕCache
        lnϕ_out = cache.lnϕ
        if nc == 1
            lnϕ_out[1] = VT_∑zlogϕ(model,vol/∑z,T,SA[1.0],p)
        else
            rhovec = @view cache.aux[1:nc]
            rhovec .= z ./ vol
            μ_res_RT,_ = __μpr(model,T,rhovec,cache)
            lnϕ_out .= μ_res_RT .- logZ
        end
    else
        μ_res = VT_chemical_potential_res(model, vol, T, z)
        if ismutable(μ_res)
            lnϕ_out = μ_res
            lnϕ_out .= μ_res ./ RT .- logZ
        else
            lnϕ_out = μ_res/RT .- logZ
        end
    end
    return lnϕ_out, vol
end

function lnϕ(model::IdealModel, p, T, z=SA[1.],cache = nothing;
            phase=:unknown,
            vol0=nothing,
            threaded = true,
            vol = volume(model,p,T,z;phase,vol0,threaded))

    lnϕ = FillArrays.Zeros(length(z))
    return lnϕ, vol
end

function lnϕ!(cache::∂lnϕCache, model::EoSModel, p, T, z=SA[1.];
            phase=:unknown,
            vol0=nothing,
            threaded = true,
            vol = volume(model,p,T,z;phase,vol0,threaded))

    is_idealmodel(model) && return (fill!(cache[3],0.0),vol)

    return lnϕ(model,p,T,z,cache;vol)
end

function lnϕ!(cache::Nothing, model::EoSModel, p, T, z=SA[1.];
            phase=:unknown,
            vol0=nothing,
            threaded = true,
            vol = volume(model,p,T,z;phase,vol0,threaded))
    return lnϕ(model,p,T,z;vol)
end

function lnϕ!(lnϕ_vec::AbstractVector, model::EoSModel, p, T, z=SA[1.],cache = nothing;
            phase=:unknown,
            vol0=nothing,
            threaded = true,
            vol = volume(model,p,T,z;phase,vol0,threaded))

    if is_idealmodel(model)
        lnϕ_vec .= 0
        return lnϕ_vec,vol
    end

    if cache !== nothing
        lnϕx, _, = lnϕ(model,p,T,z,cache;vol)
        lnϕ_vec .= lnϕx
    else
        RT = Rgas(model)*T
        μ_res = VT_chemical_potential_res!(lnϕ_vec, model, vol, T, z)
        Z = p*vol/RT/sum(z)
        lnϕ_vec .= μ_res ./ RT .- log(Z)
    end
    return lnϕ_vec, vol
end

function ∑zlogϕ(model::EoSModel, p, T, z=SA[1.],cache = nothing;
            phase=:unknown,
            vol0=nothing,
            threaded = true,
            vol = volume(model,p,T,z;phase,vol0,threaded))

    return VT_∑zlogϕ(model,vol,T,z,p),vol
end

#μ_res/RT, p_res/RT and the Hessian of Ψ over the nc densities.
function __μpr_hess(model, T, rhovec, cache)
    nc = length(model)
    
    Ψ = @deferred_isochoric_rhovec(__psi_res,psi_res)

    if cache === nothing
        result = DiffResults.HessianResult(rhovec)
        _result = ForwardDiff.hessian!(result, Ψ, rhovec)
    else
        aux = cache.aux
        rhovec_aux = @view aux[1:nc]
        rhovec_aux .= rhovec
        _result = ForwardDiff.hessian!(cache.result, Ψ, aux, cache.hconfig, Val{false}())
    end

    μ_res = @view DiffResults.gradient(_result)[1:nc]
    H = @view DiffResults.hessian(_result)[1:nc, 1:nc]
    p_res = dot(μ_res, rhovec) - DiffResults.value(_result)
    return μ_res, p_res, H, zero(p_res)
end

#μ_res/RT, p_res/RT, the Hessian of Ψ, and ∂(p_res/RT)/∂T; over the nc densities.
function __μprT_hess(model, T, rhovec, cache)
    nc = length(model)

    Ψ = WithContext(Base.Fix1(__psi_res,model),∂Tag{psi_res})

    if cache === nothing
        aux = vcat(rhovec, T)
        result = DiffResults.HessianResult(aux)
        _result = ForwardDiff.hessian!(result, Ψ, aux)
    else
        aux = cache.aux
        rhovec_aux = @view aux[1:nc]
        rhovec_aux .= rhovec
        aux[nc+1] = T
        _result = ForwardDiff.hessian!(cache.result, Ψ, aux, cache.hconfig, Val{false}())
    end

    ∇Ψ = DiffResults.gradient(_result)
    μ_res = @view ∇Ψ[1:nc]
    H = @view DiffResults.hessian(_result)[1:nc+1, 1:nc+1]
    p_res = dot(μ_res, rhovec) - DiffResults.value(_result)
    dpdT = -∇Ψ[nc+1]
    for i in 1:nc
        dpdT += rhovec[i]*H[i, nc+1]
    end
    return μ_res, p_res, H, dpdT
end

#fills lnϕ and a (in the cache, if available). returns p/RT, a·ρ, the Hessian view, ∂(p/RT)/∂T|ρ.
function __lnϕ_core!(model, T, rhovec, cache, mup_hess_fn::F) where F
    μ, p_res, H, dptotdT = mup_hess_fn(model, T, rhovec, cache)
    ρ = sum(rhovec)
    lnZ = log1p(p_res/ρ)
    if cache === nothing
        a = similar(μ)
        if ismutable(μ)
            μ .= μ .- lnZ
            lnϕ = μ
        else
            lnϕ = μ .- lnZ
        end
    else
        a = cache.∂P∂n
        lnϕ = cache.lnϕ
        lnϕ .= μ .- lnZ
    end

    #a = ∂(p/RT)/∂ρ = 1 + Hρ. returns a·ρ.
    aρ = zero(eltype(a))
    for i in eachindex(ρ)
        ai = one(eltype(a))
        for k in eachindex(rhovec)
            ai += H[i, k]*rhovec[k]
        end
        a[i] = ai
        aρ += ai*ρ[i]
    end

    return p_res, aρ, H, lnϕ, a, dptotdT
end

_dlnϕdT!(cache::∂lnϕCache, a, H, T, p_over_RT, aρ, dptotdT) = _dlnϕdT!(cache.∂lnϕ∂T, a, H, T, p_over_RT, aρ, dptotdT)
_dlnϕdT!(::Nothing, a, H, T, p_over_RT, aρ, dptotdT) = _dlnϕdT!(similar(a), a, H, T, p_over_RT, aρ, dptotdT)

function _dlnϕdT!(∂lnϕ∂T, a, H, T, p_over_RT, aρ, dptotdT)
    nc = length(a)
    s = -(p_over_RT/T + dptotdT)/aρ
    for i in 1:nc
        ∂lnϕ∂T[i] = H[i, nc+1] + 1/T + s*a[i]
    end
    return ∂lnϕ∂T
end

_dlnϕdnP!(cache::∂lnϕCache,a,H,RT,p_over_RT,ap,vol,n) = _dlnϕdnP!(cache.∂lnϕ∂n,cache.∂lnϕ∂P,a,H,RT,p_over_RT,ap,vol,n)

function _dlnϕdnP!(::Nothing,a,H,RT,p_over_RT,ap,vol,n)
    nc = length(a)
    ∂lnϕ∂n = similar(a,(nc,nc))
    ∂lnϕ∂P = similar(a)
    return _dlnϕdnP!(∂lnϕ∂n, ∂lnϕ∂P, a, H, RT, p_over_RT, aρ, vol, n)
end

function _dlnϕdnP!(∂lnϕ∂n, ∂lnϕ∂P, a, H, RT, p_over_RT, aρ, vol, n)
    ninv = 1/n
    nc = length(a)
    for i in 1:nc
        ai_aρ = a[i]/aρ
        ∂lnϕ∂P[i] = (ai_aρ - 1/p_over_RT)/RT
        for j in 1:nc
            ∂lnϕ∂n[i, j] = (H[i, j] - a[j]*ai_aρ)/vol + ninv
        end
    end
    return ∂lnϕ∂n, ∂lnϕ∂P
end

# Function to compute fugacity coefficient and its pressure and composition derivatives
function ∂lnϕ∂n∂P(model::EoSModel, p, T, z=SA[1.], cache = ∂lnϕ_cache(model,p,T,z,Val{false}());
            phase=:unknown,
            vol0=nothing,
            threaded = true,
            vol = volume(model,p,T,z;phase,vol0,threaded))

    nc = length(model)
    if cache === nothing
        rhovec = z ./ vol
    else
        aux = cache.aux
        rhovec = @view aux[1:nc]
        rhovec .= z ./ vol
    end

    n = sum(z)
    RT = Rgas(model)*T
    p_over_RT = p/RT
    _, aρ, H, lnϕ, a, _ = __lnϕ_core!(model, T, rhovec, cache, __μpr_hess)
    ∂lnϕ∂n, ∂lnϕ∂P = _dlnϕdnP!(cache, a, H, RT, p_over_RT, aρ, vol, n)
    return lnϕ, ∂lnϕ∂n, ∂lnϕ∂P, vol
end

# Function to compute fugacity coefficient pressure derivatives
function ∂lnϕ∂P(model::EoSModel, p, T, z=SA[1.], cache = ∂lnϕ_cache(model,p,T,z,Val{false}());
            phase=:unknown,
            vol0=nothing,
            threaded = true,
            vol = volume(model,p,T,z;phase,vol0,threaded))

    nc = length(model)
    RT = Rgas(model)*T
    pᵣ(x) = pressure_res(model, 1.0, T, @view(x[1:nc]))/RT
    ρt = sum(z)/vol

    if cache === nothing
        rhovec = zeros(Base.promote_eltype(model,p,T,z),nc)
        rhovec .= z ./ vol
        _result = ForwardDiff.gradient!(DiffResults.GradientResult(rhovec), pᵣ, rhovec)
        ∇pᵣ = @view DiffResults.gradient(_result)[1:nc]
        aρ = ρt + dot(rhovec, ∇pᵣ)
        ∂lnϕ∂P = rhovec
    else
        aux = cache.aux
        rhovec = @view aux[1:nc]
        rhovec .= z ./ vol
        gconfig = Solvers._GradientConfig(cache.hconfig)
        _result = ForwardDiff.gradient!(cache.result, pᵣ, aux, gconfig, Val{false}())
        ∇pᵣ = @view DiffResults.gradient(_result)[1:nc]
        aρ = ρt + dot(rhovec, ∇pᵣ)
        ∂lnϕ∂P = cache.∂lnϕ∂P
    end

    p_over_RT = p/RT
    for i in 1:nc
        ∂lnϕ∂P[i] = ((1 + ∇pᵣ[i])/aρ - 1/p_over_RT)/RT
    end

    return ∂lnϕ∂P, vol
end

function __dlnρdT(model, T, rhovec, p_over_RT = nothing)
    ρt = sum(rhovec)
    function ψ(_ρ, _T)
        _v = 1/_ρ
        return _ρ*ρt*a_res(model, _v, _T, rhovec)
    end
    ψx, g, Hψ = Solvers.∂2(ψ, one(T), T)
    aρ = ρt + Hψ[1, 1]
    dptotdT = Hψ[1, 2] - g[2]
    if p_over_RT === nothing
        _p_over_RT = ρt + g[1] - ψx 
    else
        _p_over_RT = p_over_RT
    end
    return -(_p_over_RT/T + dptotdT)/aρ
end

# Function to compute fugacity coefficient temperature derivatives
function ∂lnϕ∂T(model::EoSModel, p, T, z=SA[1.], cache = ∂lnϕ_cache(model,p,T,z,Val{false}());
            phase=:unknown,
            vol0=nothing,
            threaded = true,
            vol = volume(model,p,T,z;phase,vol0,threaded))

    nc = length(model)
    if cache === nothing
        rhovec = z ./ vol
    else
        rhovec = @view cache.aux[1:nc]
        rhovec .= z ./ vol
    end

    s = __dlnρdT(model, T, rhovec,p/(Rgas(model)*T))

    #T-derivative of ρₜ·a_res along the isobar, differentiated w.r.t. ρ
    function Φ(x)
        ρ = @view x[1:nc]
        ρt = sum(ρ)
        f(ε) = ρt*a_res(model, 1/(1 + s*ε), T + ε, ρ)
        return Solvers.derivative(f, zero(T))
    end

    if cache === nothing
        ∂lnϕ∂T = ForwardDiff.gradient(Φ, rhovec)
        ∂lnϕ∂T .+= 1/T + s
    else
        gconfig = Solvers._GradientConfig(cache.hconfig)
        _result = ForwardDiff.gradient!(cache.result, Φ, cache.aux, gconfig, Val{false}())
        ∂lnϕ∂T = cache.∂lnϕ∂T
        ∂lnϕ∂T .= @view(DiffResults.gradient(_result)[1:nc])
        ∂lnϕ∂T .+= 1/T + s
    end
    return ∂lnϕ∂T, vol
end

# Function to compute fugacity coefficient and its temperature, pressure and composition derivatives
function ∂lnϕ∂n∂P∂T(model::EoSModel, p, T, z=SA[1.],cache = ∂lnϕ_cache(model,p,T,z,Val{true}());
            phase=:unknown,
            vol0=nothing,
            threaded = true,
            vol = volume(model,p,T,z;phase,vol0,threaded))

    nc = length(model)
    if cache === nothing
        rhovec = z ./ vol
    else
        aux = cache.aux
        rhovec = @view aux[1:nc]
        rhovec .= z ./ vol
    end

    if cache isa ∂lnϕCache && nc == length(cache.aux)
        #not enough cache for one pass, so use two
        lnϕ, ∂lnϕ∂n, ∂lnϕ∂P, vol = ∂lnϕ∂n∂P(model,p,T,z,cache;vol)
        temp = cache.∂P∂n
        temp .= lnϕ
        _∂lnϕ∂T = ∂lnϕ∂T(model,p,T,z,cache;vol)
        for i in 1:nc #switch indice, use dpdn as cache for ∂lnϕ∂T
            ∂lnϕi∂T = _∂lnϕ∂T[i]
            lnϕi = temp[i]
            temp[i] = ∂lnϕi∂T
            lnϕ[i] = lnϕi
        end
        lnϕ, ∂lnϕ∂n, ∂lnϕ∂P, temp, vol
    end

    n = sum(z)
    RT = Rgas(model)*T
    p_over_RT = p/RT
    _, aρ, H, lnϕ, a, dptotdT = __lnϕ_core!(model, T, rhovec, cache, __μprT_hess)
    ∂lnϕ∂n, ∂lnϕ∂P = _dlnϕdnP!(cache, a, H, RT, p_over_RT, aρ, vol, n)
    ∂lnϕ∂T = _dlnϕdT!(cache, a, H, T, p_over_RT, aρ, dptotdT)
    return lnϕ, ∂lnϕ∂n, ∂lnϕ∂P, ∂lnϕ∂T, vol
end

#functions used to overload gamma-phi models
function modified_lnϕ(model, p, T, z, cache; phase = :unknown, vol0 = nothing)
    lnϕz,vz = lnϕ(model, p, T, z, cache; phase, vol0)
    if isnan(vz)
        lnϕz,vz = lnϕ(model, p, T, z, cache; phase)
    end
    return lnϕz,vz
end

modified_gibbs(model,p,T,w) = modified_gibbs(model,p,T,w,:unknown,oftype(zero(Base.promote_eltype(model,p,T,w)),NaN))
modified_gibbs(model,p,T,w,phase) = modified_gibbs(model,p,T,w,phase,oftype(zero(Base.promote_eltype(model,p,T,w)),NaN))

function modified_gibbs(model,p,T,w,phase,vol)
    if isnan(vol) || isnothing(vol)
        volw = volume(model,p,T,w,phase = phase)
    else
        volw = vol
    end
    RT = Rgas(model)*T
    g = PT_property(model,p,T,w,phase,volw,VT_gibbs_energy)
    return g/RT,volw
end

_dlnϕdn!(cache::∂lnϕCache,a,H,RT,p_over_RT,ap,vol,n) = _dlnϕdn!(cache.∂lnϕ∂n,a,H,RT,p_over_RT,ap,vol,n)

function _dlnϕdn!(::Nothing,a,H,RT,p_over_RT,ap,vol,n)
    nc = length(a)
    ∂lnϕ∂n = similar(a,(nc,nc))
    return _dlnϕdn!(∂lnϕ∂n, a, H, RT, p_over_RT, aρ, vol, n)
end

function _dlnϕdn!(∂lnϕ∂n, a, H, RT, p_over_RT, aρ, vol, n)
    ninv = 1/n
    nc = length(a)
    for i in 1:nc
        ai_aρ = a[i]/aρ
        for j in 1:nc
            ∂lnϕ∂n[i, j] = (H[i, j] - a[j]*ai_aρ)/vol + ninv
        end
    end
    return ∂lnϕ∂n
end

function modified_∂lnϕ∂n(model, p, T, z, cache; phase = :unknown, vol0 = nothing)
    vol = volume(model,p,T,z;phase,vol0,threaded = false)
    nc = length(model)
    if cache === nothing
        rhovec = z ./ vol
    else
        aux = cache.aux
        rhovec = @view aux[1:nc]
        rhovec .= z ./ vol
    end

    n = sum(z)
    RT = Rgas(model)*T
    p_over_RT = p/RT
    _, aρ, H, lnϕ, a, _ = __lnϕ_core!(model, T, rhovec, cache, __μpr_hess)
    ∂lnϕ∂n = _dlnϕdn!(cache, a, H, RT, p_over_RT, aρ, vol, n)
    return lnϕ, ∂lnϕ∂n, vol
end

function lnf(model, V, T, z,cache = nothing)
    RT = Rgas(model)*T
    n = sum(z)
    logZp = log(V/RT/n)
    nc = length(z)
    if cache isa ∂lnϕCache
        lnf = cache.lnϕ
        if nc == 1
            lnf1,p = VT_lnf_pure(model,V,T)
            lnf[1] = lnf1
        else
            rhovec = @view cache.aux[1:nc]
            rhovec .= z ./ V
            μ_res_RT,p_res_RT = __μpr(model, T, rhovec, cache)
            p = RT*(p_res_RT + n/V)
            lnf .= μ_res_RT .- logZp
        end
    else
        rhovec = z ./ V
        μ_res_RT,p_res_RT = __μpr(model, T, rhovec, cache)
        p = RT*(p_res_RT + n/V)
        if ismutable(μ_res_RT)
            lnf = μ_res_RT
            lnf .= μ_res_RT .- logZp
        else
            lnf = μ_res_RT .- logZp
        end
    end
    return lnf, p
end

function VT_lnf_pure(model,V,T)
    RT = Rgas(model)*T
    F,dFdV = f∂fdV_res(model,V,T,SA[1.0])
    p_res = -dFdV
    μ_res = F + p_res*V
    Zp = V/RT
    lnf = μ_res/RT - log(Zp)
    p = p_res + RT/V
    return lnf,p
end
