"""
    SSXYFlash{T}(;kwargs...)

Method to solve non-reactive multicomponent, two-phase flash problem, using a generalized succesive substitution formulation.

Only two phases are supported. If `K0` is `nothing`, it will be calculated via fugacity coefficients at p,T conditions.

### Keyword Arguments:
- `equilibrium` (optional) = equilibrium type ":vle" for liquid vapor equilibria, ":lle" for liquid liquid equilibria, `:unknown` if not specified.
- `p0` (optional), initial guess pressure, ignored if pressure is one of the flash specifications.
- `T0` (optional), initial guess temperature, ignored if temperature is one of the flash specifications.
- `K0` (optional), initial guess for the K-values.
- `x0` (optional), initial guess for the composition of phase x.
- `y0` = optional, initial guess for the composition of phase y.
- `vol0` = optional, initial guesses for phase x and phase y volumes.
- `tol_xy` = absolute composition tolerance
- `tol_pT` = relative temperature or pressure tolerance.
- `tol_of` = absolute objective property tolerance.
- `max_iters` = maximum number of iterations.
- `ss_iters` = maximum number of inner sucessive substitution iterations.
- `flash_result::FlashResult`: can be provided instead of `x0`,`y0` and `vol0` for initial guesses.
"""
struct SSXYFlash{P,T} <: FlashMethod
    equilibrium::Symbol
    T0::Union{P,Nothing}
    p0::Union{P,Nothing}
    K0::Union{Vector{T},Nothing}
    x0::Union{Vector{T},Nothing}
    y0::Union{Vector{T},Nothing}
    v0::Union{Tuple{T,T},Nothing}
    tol_xy::Float64
    tol_pT::Float64
    tol_of::Float64
    max_iters::Int
    ss_iters::Int
    verbose::Bool
end

function Solvers.primalval(method::SSXYFlash{P,T}) where {P,T}
    if P == Nothing
        λP = Nothing
    else
        λP = Solvers.primal_eltype(P)
    end

    if T == Nothing
        λT = Nothing
    else
        λT = Solvers.primal_eltype(P)
    end
    return SSXYFlash{λP,λT}(method.equilibrium,primalval(method.T0),primalval(method.p0),
                    primalval(method.K0),primalval(method.x0),primalval(method.y0),primalval(method.v0),
                    method.tol_xy,method.tol_pT,method.tol_of,
                    method.max_iters,method.ss_iters,method.verbose)
end

Base.eltype(method::SSXYFlash{T}) where T = T

function index_reduction(m::SSXYFlash,idx::AbstractVector)
    K0,x0,y0 = method.K0,method.x0,method.y0
    K0 !== nothing && (K0 = K0[idx])
    x0 !== nothing && (x0 = x0[idx])
    y0 !== nothing && (y0 = y0[idx])
    return SSXYFlash(;equilibrium,m.T0,m.p0,K0,x0,y0,m.v0,m.tol_xy,m.tol_pT,m.tol_of,m.max_iters,m.ss_iters,m.verbose)
end

index_reduction(m::SSXYFlash{Nothing,Nothing},idx::AbstractVector) = m

numphases(::SSXYFlash) = 2

function SSXYFlash(;equilibrium = :unknown,
                        T0 = nothing,
                        p0 = nothing,
                        K0 = nothing,
                        x0 = nothing,
                        y0 = nothing,
                        v0 = nothing,
                        tol_xy = 1e-10,
                        tol_pT = 1e-12,
                        tol_of = 1e-10,
                        max_iters = 50,
                        ss_iters = 10,
                        flash_result = nothing,
                        verbose = false)

    !(is_vle(equilibrium) | is_lle(equilibrium) | is_unknown(equilibrium))  && throw(error("invalid equilibrium specification for SSXYFlash"))
    if flash_result isa FlashResult
        comps,β,volumes = flash_result.compositions,flash_result.fractions,flash_result.volumes
        np = numphases(flash_result)
        np != 2 && incorrect_np_flash_error(SSXYFlash,flash_result)
        w1,w2 = comps[1],comps[2]
        v = (volumes[1],volumes[2])
        P00 = flash_result.data.p
        T00 = flash_result.data.T
        return SSXYFlash(;equilibrium = equilibrium,T0 = T00,p0 = P00,x0 = w1,y0 = w2,v0 = v,
                        tol_xy = tol_xy,tol_pT = tol_pT,tol_of = tol_of,
                        max_iters = max_iters,ss_iters = ss_iters,verbose = verbose)
    end

    if K0 == x0 == y0 === nothing #nothing specified
        #is_lle(equilibrium)
        T = Nothing
    else
        if !isnothing(K0) & isnothing(x0) & isnothing(y0) #K0 specified
            T = eltype(K0)
        elseif isnothing(K0) & !isnothing(x0) & !isnothing(y0)  #x0, y0 specified
            T = eltype(x0)
        else
            throw(error("invalid specification of initial points"))
        end
    end

    if T == Nothing && v0 !== nothing
        TT = Base.promote_eltype(v0[1],v0[2])
        _v0 = (v0[1],v0[2])
    elseif T != nothing && v0 !== nothing
        TT = Base.promote_eltype(one(T),v0[1],v0[2])
        _v0 = (v0[1],v0[2])
    else
        TT = T
        _v0 = v0
    end

    if T0 === nothing && p0 === nothing
        S = Nothing
    elseif T0 !== nothing && p0 !== nothing
        S = typeof(T0*p0)
    else
        S = typeof(something(T0,p0))
    end
    return SSXYFlash{S,TT}(equilibrium,T0,p0,K0,x0,y0,_v0,
                            tol_xy,tol_pT,tol_of,
                            max_iters,ss_iters,verbose)
end


function update_volume!(model,result,p = pressure(result),T = temperature(result))
    return result
end

function xy_flash(model::EoSModel,spec::FlashSpecifications,z,flash0::FlashResult,method::SSXYFlash)
    #we suppose model is already cached

    verbose = method.verbose
    verbose && @info "start of SS-XY flash with $spec"
    ∑z = sum(z)
    val1,val2 = spec.val1,spec.val2
    spec1,spec2 = spec.spec1,spec.spec2

    iter_type = if ((spec1 == temperature) & (spec2 == pressure)) || ((spec2 == temperature) & (spec1 == pressure))
        :pt #PT-flash
    elseif ((spec1 == temperature) & (spec2 != temperature)) || ((spec2 == temperature) & (spec1 != temperature))
        :tx #TX-flash: p is updated
    elseif ((spec1 == pressure) & (spec2 != pressure)) || ((spec2 == pressure) & (spec1 != pressure))
        :px #PX-flash: T is updated
    else
        throw(ArgumentError("ss_xy_flash does not support $spec"))
    end

    TT = Base.promote_eltype(model,val1,val2,z,flash0)
    spec_index = 2

    if iter_type == :pt
        if spec1 == temperature
            p,T = convert(TT,val1),convert(TT,val2)
        else
            p,T = convert(TT,val2),convert(TT,val1)
        end
        spec_obj = zero(TT)
    elseif iter_type == :tx
        p = convert(TT,pressure(flash0))
        if spec1 == temperature
            T,spec_obj = convert(TT,val1),convert(TT,val2)
        else
            T,spec_obj = convert(TT,val2),convert(TT,val1)
            spec_index = 1
        end
    else
        T = convert(TT,temperature(flash0))
        if spec1 == pressure
            p,spec_obj = convert(TT,val1),convert(TT,val2)
        else
            p,spec_obj = convert(TT,val2),convert(TT,val1)
            spec_index = 1
        end
    end

    #PTFlashWrapper needs to know this to update the liquid volume
    is_volume = spec1 == volume || spec2 == volume

    T_old = T
    p_old = p

    x = similar(z,TT)
    y = similar(z,TT)
    nc = length(model)
    x .= flash0.compositions[1]
    y .= flash0.compositions[2]
    volx = convert(TT,flash0.volumes[1])
    voly = convert(TT,flash0.volumes[2])
    phasex = identify_phase(model,flash0,1)
    phasey = identify_phase(model,flash0,2)
    non_inx = FillArrays.Fill(false,nc)
    non_inw = (non_inx,non_inx)
    phases = (phasex,phasey)
    lnK = similar(z,TT)
    lnK_old = similar(z,TT)

    K = similar(z,TT)
    K .= y ./ x
    K_old = similar(K)

    β = convert(TT,flash0.fractions[2])/sum(flash0.fractions)
    OF_old = convert(TT,NaN)
    OF = convert(TT,NaN)
    phasez = :unknown
    iz = 0
    volz = convert(TT,NaN)

    ss_status = RREq #we suppose equilibria here
    outer_status = :working
    #=
    :working
    :failure
    :trivial
    :maxiter
    =#
    dlnϕ_cache = ∂lnϕ_cache(model, val1, val2, x, Val{false}())
    max_iters = method.max_iters
    ss_iters = method.ss_iters
    tol_of = convert(TT,method.tol_of)
    tol_pT = convert(TT,method.tol_pT)
    tol_xy = convert(TT,method.tol_xy)
    lle = is_liquid(phasex) && is_liquid(phasey)
    vapour_idx = lle ? -1 : 2
    flash_result = FlashResult([x,y],[β,β],[volx,voly],FlashData(p,T,zero(TT),vapour_idx))

    verbose && @info "iter  ss_iter  status    p                T                OF_spec          pT_OF"
    for i in 1:max_iters
        error_lnK = convert(TT,Inf)
        outer_status == :failure && break
        ss_converged = false
        ss_count = 0
        for j in 1:ss_iters
            ss_status == RRTrivial && break
            ss_status == RRFailure && break
            ss_count += 1
            if error_lnK < tol_xy
                ss_converged = true
                break
            end

            lnK_old .= lnK
            lnK,volx,voly,_ = update_K!(lnK,model,p,T,x,y,z,β,(volx,voly),phases,non_inw,dlnϕ_cache)
            K_old .= K
            K .= exp.(lnK)

            ss_status = rachfordrice_status(K,z; K_tol = tol_xy)

            #verbose && @info "$it    $status   $β  $(round(error_lnK,sigdigits=4)) $K"

            if isnan(β) && ss_status != RRTrivial
                #try to save K? basically damping
                lnK .= 0.5 * lnK .+ 0.5 * lnK_old
                K .= exp.(lnK)
            end

            β_old = β
            β = rachfordrice(K, z; β0 = β, K_tol = tol_xy, verbose)
            x,y = update_rr!(K,β,z,x,y)

            error_lnK = dnorm(lnK,lnK_old,1)
        end

        flash_result.volumes[1] = volx
        flash_result.volumes[2] = voly
        flash_result.fractions[1] = ∑z*(1 - β)
        flash_result.fractions[2] = ∑z*β

        ss_status == RRTrivial && break
        ss_status == RRFailure && break
        OF_old = OF
        
        #calculate liquid volumes with PTFlashWrapper
        is_volume && update_volume!(model,result,p,T)

        if spec_index == 1
            spec_i = spec1(model,flash_result,iz)
        elseif spec_index == 2
            spec_i = spec2(model,flash_result,iz)
        else
            spec_i = zero(OF)
        end
        OF = spec_i - spec_obj
        OF_pT = zero(OF)

        if iter_type == :pt
            #do nothing
        elseif iter_type == :tx
            if i == 1
                p_old = p
                #starting point, we just move inside the flash
                if is_vapour(phasey) && β > 0.5
                    p += sqrt(tol_pT)*p
                else
                    p -= sqrt(tol_pT)*p
                end
            else
                dOF = (OF - OF_old)/(p - p_old) #just secant, maybe another method could be better?
                p = p - OF/dOF
            end
            OF_pT = p - p_old
        else #px    
            if i == 1
                T_old = T
                #starting point, we just move inside the flash
                if is_vapour(phasey) && β > 0.5
                    T -= sqrt(tol_pT)*T
                else
                    T += sqrt(tol_pT)*T
                end
            else
                dOF = (OF - OF_old)/(T - T_old) #just secant, maybe another method could be better?
                T_old = T
                T = T - OF/dOF
            end
            OF_pT = T - T_old
            update_temperature!(model,T)
        end
        !isfinite(OF) && (outer_status = :failure)
        !isfinite(T) && (outer_status = :failure)
        !isfinite(p) && (outer_status = :failure)
        ss_status == RRFailure && (outer_status = :failure)

        abs(OF) < tol_of && ss_converged && break
        max(abs(T - T_old),abs(p - p_old)) < tol_pT && iter_type != :pt && ss_converged && break
        i == max_iters && (outer_status = :maxiter)

        if ss_status == RRTrivial && outer_status != :failure
            if is_unknown(phasez) && outer_status != :trivial
                phasez = identify_phase(model,p,T,z)
                x .= z
                y .= z
                x ./= ∑z
                y ./= ∑z
                vapour_idx = 2
                if is_liquid(phasez)
                    flash_result.fractions[1] = 0
                    flash_result.fractions[2] = ∑z
                else
                    flash_result.fractions[1] = ∑z
                    flash_result.fractions[2] = 0
                end
                outer_status = :trivial
            end
            volz = volume(model,p,T,z,phase = phasez,vol0 = volz)/∑z
            flash_result.volumes .= volz
        end

        if verbose
            bb = flash_result.fractions[2]/sum(flash_result.fractions)
            @info "$(__pad_val(i,4))  $(__pad_val(ss_count,4))     $outer_status   $(__pad_val(p,16)) $(__pad_val(T,16)) $(__pad_val(OF,16)) $(__pad_val(OF_pT,16))"
        end

        flash_result = FlashResult(flash_result.compositions,flash_result.fractions,flash_result.volumes,FlashData(p,T,zero(TT),vapour_idx))
    end
    if outer_status == :failure
        flash_result.compositions[1] .= NaN
        flash_result.compositions[2] .= NaN
        flash_result.volumes .= NaN
        flash_result.fractions .= NaN
        g = convert(TT,NaN)
    else
        g,_ = modified_gibbs(model,flash_result)
    end

    flash_result = FlashResult(flash_result.compositions,flash_result.fractions,flash_result.volumes,FlashData(p,T,g,vapour_idx))
    update_volume!(model,flash_result)
    return flash_result
end

export SSXYFlash
