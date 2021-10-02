#aproximates liquid volume at a known pressure and t,
#by using isothermal compressibility


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
    _0 = zero(p+T+first(z))
    _nan = _0/_0
    logV0 = log(V0)
    lb_v = lb_volume(model,z)
    function logstep(_V)
        _V < log(lb_v) && return zero(_V)/zero(_V)
        _V = exp(_V)
        _p,dpdV = p∂p∂V(model,_V,T,z)
        _Δ = (p-_p)/(_V*dpdV)
        return _Δ
    end

    function f_fixpoint(_V)
        Δ = logstep(_V)
        _V + sign(Δ)*abs(Δ)^(one(Δ)+Δ)
    end

    Δ0 = logstep(logV0)
    Δ00 = Δ0
    Δ1 = logstep(logV0+Δ0)
    logV1 = logV0
    iter = 0

    #convergence hill. we do Successive substitution until we land on a
    #good initial point for aitken acceleration
    for i in 1:20
        if abs(Δ00) > abs(Δ1)
            break
        end
        if !isfinite(Δ1)
            return _nan
        end
        Δ0 = Δ1
        Δ1 = logstep(logV1+Δ0)
        logV1 = logV1 + Δ1
        iter +=1
    end
    if iter == 20
        return _nan
    end
    res = @nan(Solvers.fixpoint(f_fixpoint,logV1,Solvers.SSFixPoint(),rtol = 1e-12,max_iters=max_iters),_nan)

    return exp(res)
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
    return volume_virial(B,p,T,z)
end

function volume_virial(B::Real,p,T,z=SA[1.0])
    _0 = zero(B)
    B > _0 && return _0/_0
    a = p/(R̄*T)
    b = -1
    c = -B
    Δ = b*b-4*a*c
    n = sum(z)
    if Δ <= 0
        #virial approximation could not be calculated
        #return value at spinodal
        return -2*B
    end
    #only the left root has physical meaning
    return (-b + sqrt(b*b-4*a*c))/(2*a)
end

#(z = pV/RT)
#(RT/p = V/z)
"""
    volume(model::EoSModel,p,T,z=SA[1.0];phase=:unknown,threaded=true)

calculates the volume (m³) of the compound modelled by `model` at a certain pressure,temperature and moles.

`phase` is a Symbol that determines the initial volume root to look for:

- If `phase =:liquid` it will return the volume of the phase using a liquid initial point.

- If `phase =:vapor` it will return the volume of the phase using a gas initial point.

The default is `phase =:unknown`. with this, both liquid and volume roots will be calculated, and 
the phase with the least amount of energy is returned.

The calculation of both volume roots can be calculated in serial (`threaded=false`) or in parallel (`threaded=true`)

"""
function volume(model::EoSModel,p,T,z=SA[1.0];phase=:unknown,threaded=true)

    fp(_V) = log(pressure(model,_V,T,z)/p)

#Threaded version
    phase = Symbol(phase)
    if phase != :unknown
        V0 = x0_volume(model,p,T,z,phase=phase)
        #return Solvers.ad_newton(fp,Vg0)
        return volume_compress(model,p,T,z,V0=V0)
    end

    if threaded     
            Vg0 = x0_volume(model,p,T,z,phase=:v)
            Vl0 = x0_volume(model,p,T,z,phase=:l)
            _Vg = Threads.@spawn volume_compress(model,$p,$T,$z;V0=Vg0)
            _Vl = Threads.@spawn volume_compress(model,$p,$T,$z;V0=$Vl0)
            Vg = fetch(_Vg)     
            Vl = fetch(_Vl)
    else
            Vg0 = x0_volume(model,p,T,z,phase=:v)
            Vl0 = x0_volume(model,p,T,z,phase=:l)
            Vg =  volume_compress(model,p,T,z,V0=Vg0)
            Vl =  volume_compress(model,p,T,z,V0=Vl0)        
    end

# Serial version
#=
    Vg0 = volume_virial(model,p,T,z)
    Vg =Solvers.ad_newton(fp,Vg0,rtol=1e-08)

    Vl =  volume_compress(model,p,T,z)
=#
    function gibbs(V)
        _df,f =  ∂f(model,V,T,z)
        dV,dt = _df
        if abs((p+dV)/p) > 0.03
            return Inf
        else
            return f  +p*V
        end
    end
    #this catches the supercritical phase as well

    isnan(Vl) && return Vg
    isnan(Vg) && return Vl
    (isnan(Vl) & isnan(Vg)) && error("Failed to converge to a root")

    gg = gibbs(Vg)
    gl = gibbs(Vl)
    #@show Vg,Vl
    #@show gg,gl
    return ifelse(gg<gl,Vg,Vl)

end

function volume(model::IdealModel,p,T,z=SA[1.0];phase=:unknown,threaded=false)
    return sum(z)*R̄*T/p
end

export volume
