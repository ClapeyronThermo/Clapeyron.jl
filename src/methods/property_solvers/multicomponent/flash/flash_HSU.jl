function vx_flash_pure(model,x,v,z,spec::F,T0) where F
    crit = crit_pure(model)
    Tc,Pc,Vc = crit
    ∑z = sum(z)
    z1 = SA[oneunit(∑z)] 
    x̄,v̄ = x/∑z,v/∑z
    xc = spec_to_vt(model,v̄,Tc,z1,spec)
    _1 = oneunit(xc)
    if xc < x #we can get to the spec raising the temperature
        T = v_HSU_solve1(model,x,v̄,spec,Tc)
        p = pressure(model,v̄,T,z1)
        FlashResult([z1],[∑z],[_1*v̄],FlashData(p,T))
    end

    if T0 !== nothing
        res = v_HSU_solve2(model,x̄,v̄,spec,T0)
    else
        _T0,inside = v_HSU_T0(model,x̄,v̄,spec,crit)
        if !inside
            T = v_HSU_solve1(model,x̄,v̄,spec,_T0)
            p = pressure(model,v̄,T,z1)
            data = FlashData(p,T)
            return FlashResult([z1],[∑z],[_1*v̄],FlashData(p,T))
        else
            res = v_HSU_solve2(model,x̄,v̄,spec,_T0)
            T,vl,vv,p = res
            if vv == vl
                return FlashResult([z1],[∑z],[_1*v̄],FlashData(p,T))
            else
                βv = (v̄ - vl)/(vv - vl)
                volumes = [vl,vv]
                return FlashResult(model,p,T,z,[z1,z1],[∑z-∑z*βv,∑z*βv],volumes;sort = false)
            end
        end
    end
end

function v_HSU_T0(model,x,v,spec::F,crit) where F
    Tc,Pc,Vc = crit
    z1 = SA[1.0]
    tentative_phase = v <= Vc ? :liquid : :vapour
    fsat = Base.Fix1(saturation_pressure,model)
    #phase 1
    Tmax = Tc*one(Base.promote_eltype(model,v,x))
    Tmin = zero(Tmax)
    
    Ti = 0.9*Tmax
    inside = false
    for i in 1:20
        sat,dsat = Solvers.f∂f(fsat,Ti)
        ps,vl,vv = sat
        dpsdT,dvldT,dvvdT = dsat
        if is_liquid(tentative_phase)
            vx,dvxdT = vl,dvldT
        else
            vx,dvxdT = vv,dvvdT
        end
        xl = spec_to_vt(model,vl,Ti,z1,spec)
        xv = spec_to_vt(model,vv,Ti,z1,spec)
        xx = spec_to_vt(model, v,Ti,z1,spec)
        βxv = (x - xl)/(xv - xl)
        if vl <= v <= vv
            #note, that if the energy value is outside, we are sure that it is outside
            #but if the energy value is inside, maybe it is outside at a higher temp
            inside = (0 <= βxv <= 1)
            Tmin = Ti
            break
        end
        #TODO: find better updating scheme
        Tmax = Ti
        dT = (v - vx)/dvxdT
        Ti_newton = Ti + (v - vx)/dvxdT
        if Tmin <= Ti_newton <= Tmax && abs(dT) <= 0.1Tmax
            Ti = Ti_newton
        else
            Ti = 0.95Tmax
        end
    end
    return Ti,inside
end


function v_HSU_solve1(model,x,v,spec::F,T0) where F
    f(T) = spec_to_vt(model,v,T,SA[1.0],spec) - x
    prob = Roots.ZeroProblem(Solvers.to_newton(f),T0)
    return Roots.solve(prob,Roots.Newton())
end

function v_HSU_solve2(model,x,v,spec::S,T0) where S
    _,vl0,vv0 = saturation_pressure(model,T0)
    V0 = svec3(log(T0),log(vl0),log(vv0),log(T0))
    ps,mus = equilibria_scale(model)
    f(_x) = μpx_equality1(model,exp(_x[2]),exp(_x[3]),v,exp(_x[1]),ps,mus,x,spec)
    sol = Solvers.nlsolve2(f,V0,Solvers.Newton2Var(),NEqOptions())
    T,vl,vv = exp(sol[1]),exp(sol[2]),exp(sol[3])
    β = (v - vl)/(vv - vl)
    if 0 <= β <= 1
        return T,vl,vv,pressure(model,vv,T)
    else
        Tx = v_HSU_solve1(model,x,v,spec,T)
        return T,v,v,pressure(model,v,Tx)
    end
end


function μpx_equality1(model,v1,v2,v,T,ps,μs,s,spec::typeof(entropy))
    z = SA[1.0]
    a1,av1,at1 = ∂f_vec(model,v1,T,z)
    a2,av2,at2 = ∂f_vec(model,v2,T,z)
    g1,g2 = a1 - v1*av1,a2 - v2*av2
    Fμ = (g1 - g2)/μs
    Fp = (av2 - av1)/ps
    β1 = (v - v2)/(v1 - v2) #one on v = v1,zero on v = v2
    s1,s2 = -at1,-at2
    Fs = (β1*s1 + (1-β1)*s2 - s)/T
    F = SVector(Fμ,Fp,Fs)
    return F
end

function μpx_equality1(model,v1,v2,v,T,ps,μs,u,spec::typeof(internal_energy))
    z = SA[1.0]
    a1,av1,at1 = ∂f_vec(model,v1,T,z)
    a2,av2,at2 = ∂f_vec(model,v2,T,z)
    g1,g2 = a1 - v1*av1,a2 - v2*av2
    Fμ = (g1 - g2)/μs
    Fp = (av2 - av1)/ps
    β1 = (v - v2)/(v1 - v2) #one on v = v1,zero on v = v2
    u1,u2 = a1 - T*at1,a2 - T*at2
    Fu = (β1*u1 + (1-β1)*u2 - u)/μs
    return SVector(Fμ,Fp,Fu)
end

function μpx_equality1(model,v1,v2,v,T,ps,μs,h,spec::typeof(enthalpy))
    z = SA[1.0]
    a1,av1,at1 = ∂f_vec(model,v1,T,z)
    a2,av2,at2 = ∂f_vec(model,v2,T,z)
    g1,g2 = a1 - v1*av1,a2 - v2*av2
    Fμ = (g1 - g2)/μs
    Fp = (av2 - av1)/ps
    β1 = (v - v2)/(v1 - v2) #one on v = v1,zero on v = v2
    h1,h2 = a1 - T*at1 - v1*av1,a2 - T*at2 - v2*av2
    Fh = (β1*h1 + (1-β1)*h2 - h)/μs
    return SVector(Fμ,Fp,Fh)
end


"""
    result = uv_flash(model, u, v, n, method::FlashMethod = GeneralizedXYFlash())
    result = uv_flash(model, u, v, n; kwargs...)

Routine to solve non-reactive two-phase multicomponent flash problem. with U-V specifications.
Wrapper around [Clapeyron.xy_flash](@ref), with automatic initial point calculations. 
Inputs:
 - `u`, internal energy
 - `v`, total volume
 - `z`, vector of number of moles of each species

All keyword arguments are forwarded to [`GeneralizedXYFlash`](@ref).

 Outputs:
 - `result`, a [`FlashResult`](@ref) struct containing molar fractions, vapour fractions, molar volumes and the equilibrium temperature and pressure.
"""
function uv_flash end

function uv_flash(model::EoSModel,u,v,z =SA[1.0];kwargs...)
    method = init_preferred_method(uv_flash,model,kwargs)
    return uv_flash(model,u,v,z,method)
end

function init_preferred_method(method::typeof(uv_flash),model::EoSModel,kwargs) 
    GeneralizedXYFlash(;kwargs...)
end

function uv_flash(model,u,v,z,method::FlashMethod)
    check_arraysize(model,z)
    if supports_reduction(method)
        model_r,idx_r = index_reduction(model,z)
        z_r = z[idx_r]
        method_r = index_reduction(method,idx_r)
    else
        model_r,idx_r = model,trues(length(model))
        method_r,z_r = method,z
    end
    if length(model_r) == 1
        T0 = hasfield(typeof(method),:T0) ? method.T0 : nothing
        result1 = vx_flash_pure(model,u,v,z,internal_energy,T0)
        return index_expansion(result1,idx_r)
    end

    result = uv_flash_impl(model_r,u,v,z_r,method_r)
    if !issorted(result.volumes)
        #this is in case we catch a bad result.
        result = FlashResult(result)
    end
    ∑β = sum(result.fractions)
    result.fractions ./= ∑β
    result.fractions .*= sum(z)
    return index_expansion(result,idx_r)
end

function uv_flash_impl(model,u,v,z,method::GeneralizedXYFlash)
    flash0 = px_flash_x0(model,u,v,z,internal_energy,method)
    isone(numphases(flash0)) && return flash0
    spec = FlashSpecifications(internal_energy,u,volume,v)
    return xy_flash(model,spec,z,flash0,method)
end

export uv_flash
