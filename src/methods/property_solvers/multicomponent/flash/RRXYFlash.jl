"""
    RRXYFlash{T}(;kwargs...)

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
struct RRXYFlash{P,T} <: FlashMethod
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

function Solvers.primalval(method::RRXYFlash{P,T}) where {P,T}
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
    return RRXYFlash{λP,λT}(method.equilibrium,primalval(method.T0),primalval(method.p0),
                    primalval(method.K0),primalval(method.x0),primalval(method.y0),primalval(method.v0),
                    method.tol_xy,method.tol_pT,method.tol_of,
                    method.max_iters,method.ss_iters,method.verbose)
end

Base.eltype(method::RRXYFlash{T}) where T = T

function index_reduction(m::RRXYFlash,idx::AbstractVector)
    K0,x0,y0 = m.K0,m.x0,m.y0
    K0 !== nothing && (K0 = K0[idx])
    x0 !== nothing && (x0 = x0[idx])
    y0 !== nothing && (y0 = y0[idx])
    return RRXYFlash(;m.equilibrium,m.T0,m.p0,K0,x0,y0,m.v0,m.tol_xy,m.tol_pT,m.tol_of,m.max_iters,m.ss_iters,m.verbose)
end

index_reduction(m::RRXYFlash{Nothing,Nothing},idx::AbstractVector) = m

numphases(::RRXYFlash) = 2

function RRXYFlash(;equilibrium = :unknown,
                        T0 = nothing,
                        p0 = nothing,
                        K0 = nothing,
                        x0 = nothing,
                        y0 = nothing,
                        v0 = nothing,
                        tol_xy = 1e-10,
                        tol_pT = 1e-12,
                        tol_of = 1e-10,
                        max_iters = 100,
                        ss_iters = 5,
                        flash_result = nothing,
                        verbose = false)

    !(is_vle(equilibrium) | is_lle(equilibrium) | is_unknown(equilibrium))  && throw(error("invalid equilibrium specification for RRXYFlash"))
    if flash_result isa FlashResult
        comps,β,volumes = flash_result.compositions,flash_result.fractions,flash_result.volumes
        np = numphases(flash_result)
        np != 2 && incorrect_np_flash_error(RRXYFlash,flash_result)
        w1,w2 = comps[1],comps[2]
        v = (volumes[1],volumes[2])
        P00 = flash_result.data.p
        T00 = flash_result.data.T
        return RRXYFlash(;equilibrium = equilibrium,T0 = T00,p0 = P00,x0 = w1,y0 = w2,v0 = v,
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
    return RRXYFlash{S,TT}(equilibrium,T0,p0,K0,x0,y0,_v0,
                            tol_xy,tol_pT,tol_of,
                            max_iters,ss_iters,verbose)
end

function update_volume!(model,result,p = pressure(result),T = temperature(result))
    return result
end

function ss_xy_flash_standardize_specs(specs::FlashSpecifications{typeof(pressure),<:Any,T2,<:Any}) where T2
    if specs.spec2 == temperature
        return :pt,specs
    else
        return :px,specs
    end
end

function ss_xy_flash_standardize_specs(specs::FlashSpecifications{typeof(temperature),<:Any,T2,<:Any}) where T2
    if specs.spec2 == pressure
        return :pt,FlashSpecifications(specs.spec2,specs.val2,specs.spec1,specs.val1)
    else
        return :tx,specs
    end
end

function ss_xy_flash_standardize_specs(specs::FlashSpecifications)
    return ss_xy_flash_standardize_specs(FlashSpecifications(specs.spec2,specs.val2,specs.spec1,specs.val1))
end

function xy_flash(model::EoSModel,spec::FlashSpecifications,z,flash0::FlashResult,method::RRXYFlash)
    #we suppose model is already cached

    verbose = method.verbose
    verbose && @info "start of SS-XY flash with $spec"
    ∑z = sum(z)

    TT = Base.promote_eltype(model,spec.val1,spec.val2,z,flash0)
    nan = convert(TT,NaN)
    iter_type,norm_spec = ss_xy_flash_standardize_specs(spec)
    if iter_type == :pt
        p,T = convert(TT,norm_spec.val1),convert(TT,norm_spec.val2)
    elseif iter_type == :tx
        p,T = convert(TT,pressure(flash0)),convert(TT,norm_spec.val1)
    else #px
        p,T = convert(TT,norm_spec.val1),convert(TT,temperature(flash0))
    end

    spec_obj = norm_spec.val2
    spec_function = norm_spec.spec2

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
    outer_lnK_old = similar(z,TT)

    β = convert(TT,flash0.fractions[2])/sum(flash0.fractions)
    OF_old = nan
    OF = nan

    #used if the phase is trivial
    phasez = :unknown
    iz = 0
    volz = nan

    #we suppose equilibria, the solver then will see if we are in eq or not.
    ss_status = RREq

    outer_status = :initial
    #=
    iteration states:

    :initial : initial iterations: dampened secant updates to p/T
    :failure : NaN appeared
    :trivial : one phase
    :maxiter : last iteration
    :bounded : solver has found bounds, so solution is guaranteed
    =#

    #initial root-finding state
    OF_state = Solvers.solve1_initial_state(TT)

    dlnϕ_cache = ∂lnϕ_cache(model, p, T, x, Val{false}())
    max_iters = method.max_iters
    ss_iters = method.ss_iters
    total_ss_iters = 0
    total_outer_iters = 0

    tol_of = convert(TT,method.tol_of)
    tol_pT = convert(TT,method.tol_pT)
    tol_xy = convert(TT,method.tol_xy)
    lle = is_liquid(phasex) && is_liquid(phasey)
    vapour_idx = lle ? -1 : 2
    flash_result = FlashResult([x,y],[β,β],[volx,voly],FlashData(p,T,zero(TT),vapour_idx))

    verbose && @info "________________________________________________________________________________________
      iter  ss_iter  status    p                T                OF_spec          pT_OF"
    for i in 1:max_iters
        error_lnK = convert(TT,Inf)
        outer_status == :failure && break
        ss_converged = false
        ss_count = 0
        total_outer_iters += 1
        outer_lnK_old .= lnK

        #inner Rachford-Rice iterations. adapted from MichelsenTPFlash
        for j in 1:ss_iters
            ss_status == RRTrivial && break
            ss_status == RRFailure && break

            if error_lnK < tol_xy
                ss_converged = true
                break
            end
            ss_count += 1
            total_ss_iters += 1

            lnK_old .= lnK
            lnK,volx,voly,_ = update_K!(lnK,model,p,T,x,y,z,β,(volx,voly),phases,non_inw,dlnϕ_cache)
            K_old .= K
            K .= exp.(lnK)

            ss_status = rachfordrice_status(K,z; K_tol = tol_xy)

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

        #calculate liquid volumes with PTFlashWrapper, if necessary
        spec_function === volume && update_volume!(model,result,p,T)

        #Objective function
        OF = spec_function(model,flash_result,iz) - spec_obj
        α_lnK = zero(OF)
        w̄,w̄_old,w,w_old,w̄_new = nan,nan,nan,nan,nan

        if iter_type == :tx
            lnp,lnp_old = log(p),log(p_old)
            w̄,w̄_old,w,w_old = lnp,lnp_old,p,p_old
        elseif iter_type == :px
            τ,τ_old = 1/T,1/T_old
            w̄,w̄_old,w,w_old = τ,τ_old,T,T_old
        end
        p_old,T_old = p,T

        if i == 1
            #direction to move in first iter. move towards middle of vapour fraction, if possible.
            #by default, it is defined in pressure terms (high pressure -> less vapour)
            #we invert it if the temperature is the variable being iterated.
            _sign = (is_vapour(phasey) && β > 0.5) ? 1 : -1
            iter_type == :tx && (_sign *= -1)
            OF_state = Solvers.solve1_update_state(OF_state,w̄,OF)
            w_new = w + sqrt(tol_pT)*_sign*w
            w̄_new = w_new #not used here
        else
            w̄_new,OF_state = Solvers.solve1_new_iter(OF_state,w̄,OF,full_iter = ss_converged)
            α_lnK = (w̄_new - w̄_old)/(w̄ - w̄_old)
        end

        if iter_type == :px
            w_old = T
            T = w = i == 1 ? w̄_new : 1/w̄_new
            update_temperature!(model,T)
        elseif iter_type == :tx
            w_old = p
            p = w = i == 1 ? w̄_new : exp(w̄_new)
        end

        OF_pT = (w - w_old)/w

        #linear interpolation to suggest new K.
        if outer_status == :bounded
            lnK .= α_lnK .* lnK .+ (1 - α_lnK) .* outer_lnK_old
        end

        #outer iteration status determination
        !isfinite(OF) && (outer_status = :failure) #failure in objective
        !isfinite(w) && (outer_status = :failure) #failure in p/T
        ss_status == RRFailure && (outer_status = :failure) #failure in SS iteration
        ss_converged && i > 3  && abs(OF_pT) < sqrt(tol_pT) && (outer_status = :bounded)
        abs(OF) < tol_of && ss_converged && break
        abs(OF_pT) < tol_pT && iter_type != :pt && ss_converged && break
        i == max_iters && (outer_status = :maxiter)

        #if the iteration results in a trivial phase, keep going until the specifications are satisfied.
        if ss_status == RRTrivial && outer_status != :failure
            if is_unknown(phasez) && outer_status != :trivial
                phasez = identify_phase(model,p,T,z)
                x .= z ./ ∑z
                y .= z ./ ∑z
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

verbose &&
@info "________________________________________________________________________________________
      Final K values:        $K
      Final vapour fraction: $β
      Outer iterations:      $total_outer_iters
      SS iterations:         $total_ss_iters
      Final temperature:     $T
      Final pressure:        $p

"

    ss_status = rachfordrice_status(K,z,K_tol = tol_xy)
    verbose && ss_status != RREq && @info "result is single-phase (does not satisfy Rachford-Rice constraints)."
    volx,voly = flash_result.volumes[1],flash_result.volumes[2]
    #maybe azeotrope, do nothing in this case
    if abs(volx - voly) > sqrt(max(abs(volx),abs(voly))) && ss_status != RREq
        verbose && @info "trivial result but different volumes (maybe azeotrope?)"
        ss_status = RREq
    elseif ss_status == RREq && β <= eps(eltype(β))
        ss_status = RRLiquid
    elseif ss_status == RREq && β >= one(β) - eps(eltype(β))
        ss_status = RRVapour
    elseif !material_balance_rr_converged((x,y),z,β) #material balance failed
        verbose && @info "material balance failed."
        ss_status = RRFailure
        outer_status = :failure
    end

    verbose && ss_status == RRLiquid && @info "procedure converged to a single liquid phase."
    verbose && ss_status == RRVapour && @info "procedure converged to a single vapour phase."

    if outer_status == :failure
        flash_result.compositions[1] .= nan
        flash_result.compositions[2] .= nan
        flash_result.volumes .= nan
        flash_result.fractions .= nan
    end

    if ss_status != RREq && outer_status != :trivial #trivial system detected after the iterations
        _0 = zero(eltype(x))
        _1 = one(eltype(x))
        x .= z ./ ∑z
        y .= z ./ ∑z
        if ss_status == RRLiquid
            β = zero(TT)
            vz = volume(model,p,T,z,phase = :l)/∑z
        elseif ss_status == RRVapour
            β = one(TT)
            vz = volume(model,p,T,z,phase = :v)/∑z
        else
            β = nan
            vz = nan
        end
        flash_result.volumes[1] = vz
        flash_result.volumes[2] = vz
        flash_result.fractions[1] = ∑z*(1-β)
        flash_result.fractions[2] = ∑z*β
    end

    if outer_status != :failure
        g,_ = modified_gibbs(model,flash_result)
    else
        g = nan
    end

    flash_result = FlashResult(flash_result.compositions,flash_result.fractions,flash_result.volumes,FlashData(p,T,g,vapour_idx))
    update_volume!(model,flash_result)
    return flash_result
end

export RRXYFlash
