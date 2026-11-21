module CriticalLineContinuation
using Clapeyron, LinearAlgebra
const FD = Clapeyron.ForwardDiff
include("NormalizedCriticality.jl")
using .NormalizedCriticality: binary_constraints
export ContinuationOptions, continue_curve, trace_binary_critical_line

"""Numerical controls for a local, Float64 pseudo-arclength continuation prototype."""
Base.@kwdef struct ContinuationOptions
    ds::Float64 = 0.025
    dsmin::Float64 = 1e-7
    dsmax::Float64 = 0.05
    atol::Float64 = 1e-10
    rank_rtol::Float64 = 1e-10
    max_steps::Int = 200
    max_newton::Int = 12
    max_retries::Int = 12
    min_alpha::Float64 = 1/128
end
function validate(o)
    all(isfinite, (o.ds,o.dsmin,o.dsmax,o.atol,o.rank_rtol,o.min_alpha)) || throw(ArgumentError("Nonfinite options"))
    0 < o.dsmin <= o.ds <= o.dsmax || throw(ArgumentError("Require 0 < dsmin <= ds <= dsmax"))
    o.atol > 0 && 0 < o.rank_rtol < 1 && 0 < o.min_alpha <= 1 || throw(ArgumentError("Invalid tolerances"))
    o.max_steps >= 0 && o.max_newton > 0 && o.max_retries >= 0 || throw(ArgumentError("Invalid iteration limits"))
end
function newton(g, initial, domain, o)
    u = copy(initial)
    domain(u) || return (;ok=false,u,iterations=0,reason=:domain,residual=Inf)
    for iteration in 0:o.max_newton
        r = try g(u) catch; return (;ok=false,u,iterations=iteration,reason=:residual_error,residual=Inf) end
        nr = norm(r, Inf)
        isfinite(nr) || return (;ok=false,u,iterations=iteration,reason=:nonfinite_residual,residual=nr)
        nr <= o.atol && return (;ok=true,u,iterations=iteration,reason=:converged,residual=nr)
        iteration == o.max_newton && break
        J = try FD.jacobian(g,u) catch; return (;ok=false,u,iterations=iteration,reason=:jacobian_error,residual=nr) end
        d = try -(J\r) catch; return (;ok=false,u,iterations=iteration,reason=:singular_corrector,residual=nr) end
        all(isfinite,d) || return (;ok=false,u,iterations=iteration,reason=:nonfinite_step,residual=nr)
        alpha=1.; accepted=false
        while alpha >= o.min_alpha
            v = u+alpha*d
            if domain(v)
                nv = try norm(g(v),Inf) catch; Inf end
                if isfinite(nv) && nv < nr
                    u=v; accepted=true; break
                end
            end
            alpha /= 2
        end
        accepted || return (;ok=false,u,iterations=iteration,reason=:line_search_failed,residual=nr)
    end
    return (;ok=false,u,iterations=o.max_newton,reason=:newton_limit,residual=norm(g(u),Inf))
end
function tangent(F,u,previous,direction,coordinate,o)
    J=FD.jacobian(F,u)
    all(isfinite,J) || return nothing
    size(J) == (length(u)-1,length(u)) || throw(ArgumentError("F must have one fewer equation than state variables"))
    s=svd(J;full=true)
    s.S[end] > o.rank_rtol*s.S[1] && s.S[1]>0 || return nothing
    t=s.V[:,end]
    if previous === nothing
        # A seed at a composition fold has no preferred composition orientation.
        # Resolve deterministically on the largest tangent coordinate instead.
        k=abs(t[coordinate])>sqrt(eps(Float64)) ? coordinate : argmax(abs.(t))
        t *= direction*sign(t[k])
    elseif dot(t,previous)<0
        t *= -1
    end
    return (;vector=t,smin=s.S[end],smax=s.S[1])
end

# Bracket a simple coordinate turning point on the corrected local arc. A change
# of endpoint tangent signs can hide an exit/re-entry even when both states lie
# inside the requested bounds. No global event guarantee for multiple unresolved
# oscillations within one step is asserted.
function fold_event(F,u,v,t,tend,coordinate,bounds,domain,o)
    t[coordinate]*tend[coordinate] < 0 || return (;ok=true,extremum=nothing,event=nothing)
    lo=0.;hi=dot(v-u,t)
    hi>0 || return (;ok=false,extremum=nothing,event=nothing)
    slope=t[coordinate];w=copy(v)
    for iteration in 1:40
        mid=(lo+hi)/2;pred=u+mid*t
        g(y)=vcat(F(y),dot(y-pred,t))
        corrected=newton(g,pred,domain,o)
        corrected.ok || return (;ok=false,extremum=nothing,event=nothing)
        w=corrected.u
        tg=try tangent(F,w,t,1,coordinate,o) catch; nothing end
        tg===nothing && return (;ok=false,extremum=nothing,event=nothing)
        sw=tg.vector[coordinate]
        abs(sw)<1e-10 && break
        if slope*sw>0;lo=mid;else;hi=mid;end
        hi-lo<1e-12*max(1.,norm(u)) && break
    end
    event=w[coordinate]<=bounds[1] ? bounds[1] : w[coordinate]>=bounds[2] ? bounds[2] : nothing
    return (;ok=true,extremum=w,event)
end

"""
Continue a one-dimensional regular root set F(u)=0, preserving trajectory order.
Coordinate bounds are user-requested stopping surfaces, not physical endpoints.
The seed must already satisfy F. Status and every rejected attempt are retained.
A rank-deficient tangent or exhausted corrector returns a partial result.
"""
function continue_curve(F,u0; direction=1, coordinate=length(u0), bounds=(-Inf,Inf),
                        domain=u->all(isfinite,u), options=ContinuationOptions(), initial_tangent=nothing)
    o=options; validate(o)
    direction in (-1,1) || throw(ArgumentError("direction must be -1 or +1"))
    u=Float64.(u0)
    1<=coordinate<=length(u) || throw(ArgumentError("Invalid coordinate"))
    bounds[1]<bounds[2] || throw(ArgumentError("Invalid bounds"))
    domain(u) && bounds[1]<=u[coordinate]<=bounds[2] || throw(ArgumentError("Seed outside domain/bounds"))
    r=F(u)
    length(r)==length(u)-1 || throw(ArgumentError("F must define a curve"))
    all(isfinite,r) && norm(r,Inf)<=o.atol || throw(ArgumentError("Seed does not satisfy criticality tolerance"))
    previous=initial_tangent===nothing ? nothing : Float64.(initial_tangent)
    if previous!==nothing
        length(previous)==length(u) && all(isfinite,previous) && norm(previous)>0 || throw(ArgumentError("Invalid initial tangent"))
        previous /= norm(previous)
    end
    states=[copy(u)]; diagnostics=NamedTuple[]; attempts=NamedTuple[]; ds=o.ds
    result(status)=(;states,diagnostics,attempts,status)
    for step in 1:o.max_steps
        tg=try tangent(F,u,previous,direction,coordinate,o) catch e
            e isa ArgumentError && rethrow()
            return result(:tangent_error)
        end
        tg===nothing && return result(:singular_tangent)
        t=tg.vector; accepted=nothing
        for retry in 0:o.max_retries
            ds<o.dsmin && break
            pred=u+ds*t
            g(v)=vcat(F(v),dot(v-pred,t))
            candidate=newton(g,pred,domain,o)
            push!(attempts,(;step,retry,ds,ok=candidate.ok,reason=candidate.reason,residual=candidate.residual))
            if candidate.ok
                accepted=candidate; break
            end
            ds/=2
        end
        accepted===nothing && return result(ds<o.dsmin ? :minimum_step : :corrector_failed)
        v=accepted.u; event=nothing; event_bracket_end=v
        if v[coordinate]<=bounds[1]; event=bounds[1]
        elseif v[coordinate]>=bounds[2]; event=bounds[2]
        end
        if any(isfinite,bounds)
            tg_end=try tangent(F,v,t,direction,coordinate,o) catch; nothing end
            tg_end===nothing && return result(:event_tangent_failed)
            fold=fold_event(F,u,v,t,tg_end.vector,coordinate,bounds,domain,o)
            fold.ok || return result(:event_bracket_failed)
            if fold.event!==nothing
                event=fold.event;event_bracket_end=fold.extremum
            end
        end
        if event!==nothing
            # Correct a bracketed coordinate crossing. Bounds are requested
            # stopping surfaces, not physical/pure endpoints. A simple fold is
            # bracketed above; multiple hidden oscillations are not certified.
            fraction=(event-u[coordinate])/(event_bracket_end[coordinate]-u[coordinate])
            guess=u+clamp(fraction,0,1)*(event_bracket_end-u)
            ge(w)=vcat(F(w),w[coordinate]-event)
            endpoint=newton(ge,guess,domain,o)
            push!(attempts,(;step,retry=-1,ds,ok=endpoint.ok,reason=endpoint.ok ? :endpoint_converged : endpoint.reason,residual=endpoint.residual))
            endpoint.ok || return result(:endpoint_correction_failed)
            v=endpoint.u
        end
        norm(v-u)>eps(Float64)*max(1,norm(u)) || return result(:no_progress)
        push!(states,copy(v))
        push!(diagnostics,(;step,ds=norm(v-u),newton_iterations=accepted.iterations,residual=norm(F(v),Inf),tangent=copy(t),smin=tg.smin,smax=tg.smax,event))
        u=v; previous=t
        event!==nothing && return result(:coordinate_boundary)
        ds=min(o.dsmax,ds*(accepted.iterations<=4 ? 1.2 : 1.0))
    end
    return result(:max_steps)
end

"""
Trace a binary local critical line from (T, p, V) at mole fraction x0.
Returns both directions in arclength order. Uses eigenvalue-normalized binary residuals by default; residual_mode=:raw
retains unscaled determinant equations. Normalized and raw tolerances differ.
stop_coordinate selects composition, temperature or density boundary correction.
Branches are oriented by that coordinate, not always composition. Other bounds
remain hard guards; this is not simultaneous earliest-boundary detection. This prototype does not certify global phase
stability, discover all branches, cross multiple critical modes, or include x=0/1.
"""
function trace_binary_critical_line(model; x0=0.5, seed=nothing, xbounds=(0.005,0.995),
                                    Tbounds=(80.,1500.), rhobounds=(1.,60000.), options=ContinuationOptions(), residual_mode=:normalized, stop_coordinate=:composition)
    stop_coordinate in (:composition,:temperature,:density) || throw(ArgumentError("stop_coordinate must be :composition, :temperature or :density"))
    residual_mode in (:normalized,:raw) || throw(ArgumentError("residual_mode must be :normalized or :raw"))
    length(model)==2 || throw(ArgumentError("Only binary models are supported"))
    0 < xbounds[1] < xbounds[2] < 1 || throw(ArgumentError("Composition bounds must be strictly interior"))
    xbounds[1] < x0 < xbounds[2] || throw(ArgumentError("Seed composition must be interior to bounds"))
    Tbounds[1]>0 && Tbounds[1]<Tbounds[2] && rhobounds[1]>0 && rhobounds[1]<rhobounds[2] || throw(ArgumentError("Invalid positive T/density bounds"))
    s=seed===nothing ? crit_mix(model,[x0,1-x0]) : seed
    length(s)==3 && all(isfinite,s) && s[1]>0 && s[3]>0 || throw(ArgumentError("Invalid critical seed"))
    Tscale=s[1]; rhoscale=1/s[3]
    constraints(u)=binary_constraints(model,1/(rhoscale*exp(u[2])),Tscale*exp(u[1]),[u[3],1-u[3]])
    F(u)=getproperty(constraints(u),residual_mode)
    # Selected bounds are stopping surfaces, not a hard domain guard. This lets
    # the local predictor bracket a crossing before correcting onto it exactly.
    # Unselected bounds remain numerical guards; no first-of-all-surfaces claim.
    coord=stop_coordinate==:temperature ? 1 : stop_coordinate==:density ? 2 : 3
    selected_bounds=coord==1 ? log.(Tbounds./Tscale) : coord==2 ? log.(rhobounds./rhoscale) : xbounds
    function domain(u)
        all(isfinite,u) && 0<u[3]<1 || return false
        coord==3 || xbounds[1]<u[3]<xbounds[2] || return false
        coord==1 || log(Tbounds[1]/Tscale)<u[1]<log(Tbounds[2]/Tscale) || return false
        coord==2 || log(rhobounds[1]/rhoscale)<u[2]<log(rhobounds[2]/rhoscale) || return false
        return true
    end
    u0=[0.,0.,Float64(x0)]
    branches=[continue_curve(F,u0;direction=d,coordinate=coord,bounds=selected_bounds,domain,options) for d in (-1,1)]
    function physical(u)
        c=constraints(u); T=Tscale*exp(u[1]); rho=rhoscale*exp(u[2])
        return (;T,rho,x=u[3],p=pressure(model,1/rho,T,[u[3],1-u[3]]),
                residual=getproperty(c,residual_mode),raw_residual=c.raw,
                normalized_residual=c.normalized,lambda_transverse=c.lambda_transverse)
    end
    return (;branches,points=[[physical(u) for u in b.states] for b in branches],Tscale,rhoscale,residual_mode,stop_coordinate)
end
end
