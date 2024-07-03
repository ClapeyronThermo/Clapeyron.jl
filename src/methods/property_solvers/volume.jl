#aproximates liquid volume at a known pressure and t,
#by using isothermal compressibility
#dP = (α/β)dT - (1/βV)dV, dT = 0
#dP = -(1/βV)dV

"""
    volume_compress(model,p,T,z=SA[1.0];V0=x0_volume(model,p,T,z,phase=:liquid),max_iters=100)
Main routine to calculate a volume, given a pressure, temperature, composition and intitial volume guess. each step is taken by locally aproximating the EoS as an isothermal compressibility process.
The new volume is calculated by the following recurrence formula:
```julia
v[i+1] = v[i]*exp(β[i]*(p-p(v[i])))
```
In the liquid root region, the iterations follow `v0 < v[i] < v[i+1] < v(p)`, allowing the calculation of the liquid root without entering the metastable region.
"""
function volume_compress(model,p,T,z=SA[1.0];V0=x0_volume(model,p,T,z,phase=:liquid),max_iters=100)
    p,T,V0 = promote(p,T,V0)
    return _volume_compress(model,p,T,z,V0,max_iters)
end

function _volume_compress(model,_p,_T,_z=SA[1.0],V0=x0_volume(model,p,T,z,phase=:liquid),max_iters=100)
    _0 = zero(Base.promote_eltype(model,_p,_T,_z))
    _1 = one(_0)
    isnan(V0) && return _0/_0
    p₀ = primalval(_1*_p)
    XX = typeof(p₀)
    T = primalval(_T)
    _nan = primalval(_0/_0)
    nan = primalval(_nan)
    logV0 = primalval(log(V0)*_1)
    z = primalval(_z)
    log_lb_v = log(primalval(lb_volume(model,T,z)))
    function logstep(logVᵢ::TT) where TT
        logVᵢ < log_lb_v && return TT(zero(logVᵢ)/zero(logVᵢ))
        Vᵢ = exp(logVᵢ)
        _pᵢ,_dpdVᵢ = p∂p∂V(model,Vᵢ,T,z)
        pᵢ,dpdVᵢ = primalval(_pᵢ),primalval(_dpdVᵢ) #ther could be rare cases where the model itself has derivative information.
        dpdVᵢ > 0 && return TT(zero(logVᵢ)/zero(logVᵢ)) #inline mechanical stability.
        abs(pᵢ-p₀) < 3eps(p₀) && return TT(zero(Vᵢ)) #this helps convergence near critical points.
        Δᵢ = (p₀-pᵢ)/(Vᵢ*dpdVᵢ) #(_p - pset)*κ
        return TT(Δᵢ)
    end
    function f_fixpoint(logVᵢ::TT) where TT
        return TT(logVᵢ + logstep(logVᵢ))
    end

    logV = @nan(Solvers.fixpoint(f_fixpoint,logV0,Solvers.SSFixPoint(),rtol = 1e-12,max_iters=max_iters)::XX,nan)
    #netwon step to recover derivative information:
    #V = V - (p(V) - p)/(dpdV(V))
    #dVdP = -1/dpdV
    #dVdT = dpdT/dpdV
    #dVdn = dpdn/dpdV
    Vsol = exp(logV)
    psol,dpdVsol = p∂p∂V(model,Vsol,_T,_z)
    return Vsol - (psol - _p)/dpdVsol
end

"""
    volume_virial(model::EoSModel,p,T,z=SA[1.0])
    volume_virial(B::Real,p,T,z=SA[1.0])
Calculates an aproximation to the gas volume at specified pressure, volume and composition, by aproximating:
```julia
Z(v) ≈ 1 + B(T)/v
```
where `Z` is the compressibility factor and `B` is the second virial coefficient.
If `B>0`, (over the inversion temperature) returns `NaN`. If the solution to the problem is complex (`Z = 1 + B/v` implies solving a quadratic polynomial), returns `-2*B`.
If you pass an `EoSModel` as the first argument, `B` will be calculated from the EoS at the input `T`. You can provide your own second virial coefficient instead of a model.
"""
function volume_virial end

function volume_virial(model::EoSModel,p,T,z=SA[1.0])
    B = second_virial_coefficient(model,T,z)
    return volume_virial(B,p,T,z,R = Rgas(model))
end

function volume_virial(B::Real,p,T,z=SA[1.0];R = R̄)
    _0 = zero(B)

    #=
    PV/RT∑z = 1 + B/V // a = P/RT∑z, .*= V
    aV2 = V + B 
    aV2 - V - B = 0 
    =#
    B > _0 && return _0/_0
    a = p/(R *T*sum(z))
    b = -1
    c = -B
    Δ = b*b-4*a*c
    if Δ <= 0
        #virial approximation could not be calculated
        #return value at spinodal
        return -2*B
    end
    #only the left root has physical meaning

    #stable way of calculating quadratics, seems to matter here
    if b >= 0
        return 2*c/(- b - sqrt(Δ))
    else
        return (-b + sqrt(Δ))/(2*a)
    end
end

#(z = pV/RT)
#(RT/p = V/z)
"""
    volume(model::EoSModel, p, T, z=SA[1.0]; phase=:unknown, threaded=true, vol0=nothing)

Calculates the volume (m³) of the compound modelled by `model` at a certain pressure, temperature and moles.
`phase` is a Symbol that determines the initial volume root to look for:
- If `phase =:unknown` (Default), it will return the physically correct volume root with the least gibbs energy.
- If `phase =:liquid`, it will return the volume of the phase using a liquid initial point.
- If `phase =:vapor`, it will return the volume of the phase using a gas initial point.
- If `phase =:solid`, it will return the volume of the phase using a solid initial point (only supported for EoS that support a solid phase)
- If `phase =:stable`, it will return the physically correct volume root with the least gibbs energy, and perform a stability test on the result.

All volume calculations are checked for mechanical stability, that is: `dP/dV <= 0`.

The calculation of both volume roots can be calculated in serial (`threaded=false`) or in parallel (`threaded=true`).

An initial estimate of the volume `vol0` can be optionally be provided.

!!! tip
    The volume computation may fail and return `NaN` because the default initial point is too far from the actual volume.
    Providing a value for `vol0` may help in these situations.
    Such a starting point can be found from physical knowledge, or by computing the volume using a different model for example.

!!! warning "Stability checks"
    The stability check is disabled by default. that means that the volume obtained just follows the the relation `P = pressure(model,V,T,z)`.
    For single component models, this is alright, but phase splits (with different compositions that the input) can and will occur, meaning that
    the volume solution does not correspond to an existing phase.
    For unknown multicomponent mixtures, it is recommended to use a phase equilibrium procedure (like `tp_flash`) to obtain a list of valid compositions, and then perform a volume calculation over those compositions.
    You can also pass `phase=:stable` to perform the stability test inside the volume solver. Finally, you can perform the stability test after the volume solver:
    ```julia
    v = volume(model,p,T,z)
    isstable(model,v,T,z)
    ```
"""
function volume(model::EoSModel,p,T,z=SA[1.0];phase=:unknown, threaded=true,vol0=nothing)
    #this is used for dispatch on symbolic variables
    return _volume(model,p,T,z,phase,threaded,vol0)
end

function _volume(model::EoSModel,p,T,z=SA[1.0],phase=:unknown, threaded=true,vol0=nothing)
    return volume_impl(model,p,T,z,phase,threaded,vol0)
end

fluid_model(model) = model
solid_model(model) = model

function volume_impl(model::EoSModel,p,T,z=SA[1.0],phase=:unknown, threaded=true,vol0=nothing)
    return default_volume_impl(model,p,T,z,phase,threaded,vol0)
end

function default_volume_impl(model::EoSModel,p,T,z=SA[1.0],phase=:unknown, threaded=true,vol0=nothing)
#Threaded version
    check_arraysize(model,z)
    TYPE = typeof(p+T+first(z)+oneunit(eltype(model)))
    nan = zero(TYPE)/zero(TYPE)
    #err() = @error("model $model Failed to converge to a volume root at pressure p = $p [Pa], T = $T [K] and compositions = $z")
    fluid = fluid_model(model)
    solid = solid_model(model)

    if !isnothing(vol0)
        if !isnan(vol0)
            V0 = vol0
            if is_solid(phase) #to allow specification of the model.
                return _volume_compress(solid,p,T,z,V0)
            end
            V = _volume_compress(fluid,p,T,z,V0)
            if solid !== fluid && isnan(V)
                return _volume_compress(solid,p,T,z,V0)
            end
            return V
        end
    end

    if phase != :unknown && phase != :stable
        V0 = x0_volume(model,p,T,z,phase=phase)
        if is_solid(phase)
            V = _volume_compress(solid,p,T,z,V0)
        else
            V = _volume_compress(fluid,p,T,z,V0)
        end
        return V
    end
    Vg0 = x0_volume(fluid,p,T,z,phase=:v)
    Vl0 = x0_volume(fluid,p,T,z,phase=:l)
    Vs0 = x0_volume_solid(solid,T,z) #Needs to be const-propagated.
    volumes0 = (Vg0,Vl0,Vs0)
    if threaded
        #=
        ch = Channel{TYPE}(3) do ys
            Threads.@sync for v0 in volumes0
                Threads.@spawn put!(ys, _volume_compress($model,$p,$T,$z,v0)) 
            end
        end
        v1::TYPE = take!(ch)
        v2::TYPE = take!(ch)
        v3::TYPE = take!(ch) 
        volumes = (v1,v2,v3)
        =#
        _Vg = StableTasks.@spawn _volume_compress($fluid,$p,$T,$z,$Vg0)
        _Vl = StableTasks.@spawn _volume_compress($fluid,$p,$T,$z,$Vl0)
        if !isnan(Vs0)
            _Vs = StableTasks.@spawn _volume_compress($solid,$p,$T,$z,$Vs0)
        else
            _Vs = nan
        end
        Vg = fetch(_Vg)::TYPE
        Vl = fetch(_Vl)::TYPE
        Vs = fetch(_Vs)::TYPE
        volumes = (Vg,Vl,Vs)
    else
        Vg = _volume_compress(fluid,p,T,z,Vg0)
        Vl = _volume_compress(fluid,p,T,z,Vl0)
        Vs = _volume_compress(solid,p,T,z,Vs0)
        volumes = (Vg,Vl,Vs)
    end
    
    idx,v,g = volume_label((fluid,fluid,solid),p,T,z,volumes)
    if phase == :stable
        !VT_isstable(model,v,T,z,false) && return nan
    end
    return v
end

function volume_label(models::F,p,T,z,vols) where F
    function gibbs(model,fV)
        isnan(fV) && return one(fV)/zero(fV)
        f(V) = eos(model,V,T,z)
        _f,_dV = Solvers.f∂f(f,fV)
        return ifelse(abs((p+_dV)/p) > 0.03,one(fV)/zero(fV),_f + p*fV)
    end

    idx = 0
    V = p
    model = models[1]
    _0 = zero(@f(Base.promote_eltype))
    g = one(_0)/_0
    v = _0/_0
    for (i,vi) in pairs(vols)
        gi = gibbs(models[i],vi)
        if gi < g
            g = gi
            idx = i
            v = vi
        end
    end
    return idx,v,g
end

#=
used by MultiComponentFlash.jl extension
=#
function _label_and_volumes(model::EoSModel,cond)
    #gibbs comparison, the phase with the least amount of gibbs energy is the most stable.
    p,T,z = cond.p,cond.T,cond.z
    Vl = volume(model,p,T,z,phase =:l)
    Vv = volume(model,p,T,z,phase =:v)
    function gibbs(fV)
        isnan(fV) && return one(fV)/zero(fV)
        _df,_f = ∂f(model,fV,T,z)
        dV,_ = _df
        return ifelse(abs((p+dV)/p) > 0.03,zero(dV)/one(dV),_f + p*fV)
    end
    isnan(Vl) && return 1,Vv,Vv #could not converge on gas volume, assuming stable liquid phase
    isnan(Vv) && return 0,Vl,Vl #could not converge on liquid volume, assuming stable gas phase
    gl,gv = gibbs(Vl),gibbs(Vv)
    V = gv < gl ? 1 : 0
    return V,Vl,Vv
end

export volume
