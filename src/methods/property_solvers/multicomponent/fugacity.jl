##general successive substitution method to solve bubble/dew problems via fugacity coefficients
module FugEnum
    @enum BubbleDew begin
        BUBBLE_PRESSURE
        BUBBLE_TEMPERATURE
        DEW_PRESSURE
        DEW_TEMPERATURE
        LLE_PRESSURE
        LLE_TEMPERATURE
    end
    is_bubble(x::BubbleDew) = (x == BUBBLE_PRESSURE || x == BUBBLE_TEMPERATURE)
    is_dew(x::BubbleDew) = (x == DEW_PRESSURE || x == DEW_TEMPERATURE)
    is_temperature(x::BubbleDew) = (x == BUBBLE_TEMPERATURE  || x == DEW_TEMPERATURE || x == LLE_TEMPERATURE)
    is_pressure(x::BubbleDew) = !is_temperature(x)
    function phases(x::BubbleDew)
        if x == LLE_PRESSURE || x == LLE_TEMPERATURE
            return :liquid,:liquid
        else
            return :liquid,:vapour
        end
    end
end

is_lle(x::FugEnum.BubbleDew) = x == FugEnum.LLE_PRESSURE || x == FugEnum.LLE_TEMPERATURE

Base.@kwdef struct FugData{T}
    method::FugEnum.BubbleDew
    itmax_ss::Int = 5
    itmax_newton::Int = 10
    tol_pT::T = 1e-8
    tol_xy::T = 1e-8
    tol_of::T = 1e-8
    second_order::Bool = true
    verbose::Bool = false
end

FugEnum.is_pressure(data::FugData) = FugEnum.is_pressure(data.method)
FugEnum.is_bubble(data::FugData) = FugEnum.is_bubble(data.method)

function invalid_bd_error(data)
    pot = FugEnum.is_pressure(data) ? " pressure" : " temperature"

    spec = if is_lle(data)
        "LLE" * pot
    elseif FugEnum.is_bubble(data)
        "bubble" * pot
    else
        "dew" * pot
    end
    throw(error("invalid specification for $spec"))
end

function fug_bubbledew_cache(modelx,modely,p,T,x,y,val::Val{B}) where B
    TT = Base.promote_eltype(modelx,p,T,x,y)
    n1 = length(modelx)
    n2 = length(modely)
    nmin,nmax = minmax(n1,n2)
    w1 = similar(x,TT,nmin)
    w2 = similar(x,TT,nmin)
    w3 = similar(x,TT,nmin)
    w4 = similar(x,TT,nmin)
    w5 = similar(x,TT,nmin)
    w6 = similar(x,TT,nmin)
    if nmin == nmax
        w7 = w6
    else
        w7 = similar(x,TT,nmax)
    end
    volcache = Base.RefValue{Tuple{TT,TT}}()
    Hϕx = ∂lnϕ_cache(modelx, p, T, x, val)
    if n1 != n2
        Hϕy = ∂lnϕ_cache(modely, p, T, y, val)
    else
        Hϕy = Hϕx
    end
    return w1,w2,w3,w4,w5,w6,volcache,Hϕx,Hϕy,w7
end

function _fug_OF_ss(model::EoSModel,p,T,x,y,vol0,data::FugData,cache)
    volx,voly = vol0
    method = data.method
    itmax_ss = data.itmax_ss
    itmax_newton = data.itmax_newton
    tol_xy = data.tol_xy
    tol_of = data.tol_of
    tol_pT = data.tol_pT
    second_order = data.second_order
    phasex,phasey = FugEnum.phases(method)
    verbose = data.verbose
    _bubble_or_lle = FugEnum.is_bubble(method) || is_lle(method)
    _pressure = FugEnum.is_pressure(method)
    converged = false
    tol_stability = abs2(cbrt(tol_xy))
    #caches for ∂lnϕ∂n∂P∂T/∂lnϕ∂n∂P
    lnK,K,w,w_old,w_calc,w_restart,_,Hϕx = cache
    OF = NaN*zero(eltype(lnK))
    
    if _bubble_or_lle
        w .= y
        _x,_y = x,w
    else
        w .= x
        _x,_y = w,y
    end

    valid_iter = true
    T_old,p_old = T,p

    if _pressure
    verbose && @info "___________________________________________________________________________________
      newton_iters  ss_iters     p                OF               K"    
    else
    verbose && @info "___________________________________________________________________________________
      newton_iters  ss_iters     T                OF               K"   
    end

    for j in 1:itmax_newton
        w_restart .= w_calc
        volx_restart,voly_restart = volx,voly
        valid_iter = true
        if isnan(volx) || isnan(voly)
            break
        end
        error = Inf*one(eltype(lnK))
        ss_count = 0
        for i in 1:itmax_ss
            error < tol_xy && break
            ss_count += 1
            lnϕx, volx = modified_lnϕ(model, p, T, _x, Hϕx, vol0=volx, phase = phasex)
            if isnan(volx)
                lnϕx, volx = modified_lnϕ(model, 1.1p, T, _x, Hϕx, phase = phasex)
            end
            lnK .= lnϕx

            lnϕy, voly = modified_lnϕ(model, p, T, _y, Hϕx, vol0=voly, phase = phasey)
            if isnan(voly)
                lnϕy, voly = modified_lnϕ(model, p, T, _y, Hϕx, phase = phasey)
            end
            lnK .-= lnϕy

            if isnan(volx) || isnan(voly)
                break
            end

            K .= exp.(lnK)
            w_old .=  w

            if _bubble_or_lle
                w .= _x .* K
                w_calc .= w
            else
                w .= _y ./ K
                w_calc .= w
            end
            w ./= sum(w)
            error = dnorm(w,w_old,Inf) #||x-x_old||∞

            #if _bubble
            #    tpd = 1 + @sum(_y[i]*(lnϕy[i] + log(_y[i]) - log(_x[i]) - lnϕx[i] - 1))
            #else
            #    tpd = 1 + @sum(_x[i]*(lnϕx[i] + log(_x[i]) - log(_y[i]) - lnϕy[i] - 1))
            #end
            if dnorm(_x,_y,Inf) < tol_stability #the interation procedure went wrong. perform a T/P movement first
                w_calc .= w_restart
                w .= w_restart
                w ./= sum(w_restart)
                valid_iter = false
                volx,voly = volx_restart,voly_restart
                K .= y ./ x
                lnK .= log.(K)
                break
            end
        end

        OF_old = OF
        OF = sum(w_calc) - 1.0
    
        if _pressure && second_order
            ∂lnϕ∂Px, volx = ∂lnϕ∂P(model, p, T, _x, Hϕx, phase=phasex, vol0=volx)
            ∂OF = dot(∂lnϕ∂Px,w_calc)
            ∂lnϕ∂Py, voly = ∂lnϕ∂P(model, p, T, _y, Hϕx, phase=phasey, vol0=voly)
            ∂OF -= dot(∂lnϕ∂Py,w_calc)
        elseif !_pressure && second_order
            ∂lnϕ∂Tx, volx = ∂lnϕ∂T(model, p, T, _x, Hϕx, phase=phasex, vol0=volx)
            ∂OF = dot(∂lnϕ∂Tx,w_calc)
            ∂lnϕ∂Ty, voly = ∂lnϕ∂T(model, p, T, _y, Hϕx, phase=phasey, vol0=voly)
            ∂OF -= dot(∂lnϕ∂Ty,w_calc)
        elseif _pressure && !second_order
            if j == 1
                ∂OF = OF/sqrt(eps(p))
            else
                ∂OF = (OF - OF_old)/(p - p_old)
            end
        else
            if j == 1
                ∂OF = OF/sqrt(eps(T))
            else
                ∂OF = (OF - OF_old)/(T - T_old)
            end
        end
        if isnan(volx) || isnan(voly)
            break
        end

        if !_bubble_or_lle && second_order
            ∂OF = -∂OF
        end

        ∂step = OF / ∂OF
        if valid_iter && abs(∂step) < tol_pT || abs(OF) < tol_of
            converged = true
            break
        end

        if _pressure
            ∂step = clamp(∂step,-0.4*p,0.4*p)
            p_old = p
            p -= ∂step
        else
            ∂step = clamp(∂step,-0.05*T,0.05*T)
            T_old = T
            Tinv = 1/T + ∂step/(T*T)
            T = 1/Tinv
            #T -= ∂step
            update_temperature!(model,T)
        end

        propx =  _pressure ? p : T
        verbose && @info "$(__pad_val(j,13)) $(__pad_val(ss_count,12)) $(__pad_val(propx,16)) $(__pad_val(OF,16)) $(repr(K,context = :compact => true))"

        if !isfinite(∂step) #error, fail early, the NaN propagation is handled upstream
            converged = true
            break
        end
    end

    if !valid_iter
        w .= NaN
        lnK .= NaN
        volx,voly = w[1],w[1]
        if _pressure
            p = w[1]
        else
            T = w[1]
        end
    end
    if verbose && converged
        @info "algorithm converged in sucessive substitution phase."
    end
    return converged,(p,T,_x,_y,(volx,voly))
end

##general successive substitution method to solve bubble/dew problems via fugacity coefficients
##support for nonvolatiles/noncondensables

function _fug_OF_ss(modelx::EoSModel,modely::EoSModel,p,T,x,y,vol0,_view,data::FugData,cache)
    volx,voly = vol0
    converged = false

    method = data.method
    itmax_ss = data.itmax_ss
    itmax_newton = data.itmax_newton
    tol_xy = data.tol_xy
    tol_of = data.tol_of
    tol_pT = data.tol_pT
    second_order = data.second_order
    verbose = data.verbose
    _bubble_or_lle = FugEnum.is_bubble(method) || is_lle(method)
    _pressure = FugEnum.is_pressure(method)
    tol_stability = abs2(cbrt(tol_xy))

    lnK,K,w,w_old,w_calc,w_restart,_,Hϕx,Hϕy,u = cache

    OF = NaN*zero(eltype(lnK))

    if _bubble_or_lle
        w .= y
        u .= x
        _x,_y = u,w
    else
        w .= x
        u .= y
        _x,_y = w,u
    end

    if _pressure
    verbose && @info "___________________________________________________________________________________
      newton_iters  ss_iters     p                OF               K"    
    else
    verbose && @info "___________________________________________________________________________________
      newton_iters  ss_iters     T                OF               K"   
    end

    p_old,T_old = p,T
    valid_iter = false
    for j in 1:itmax_newton
        if isnan(volx) || isnan(voly)
            break
        end

        w_restart .= w_calc
        volx_restart,voly_restart = volx,voly
        valid_iter = true
        error = Inf*one(eltype(lnK))
        ss_count = 0
        for i in 1:itmax_ss
            error < tol_xy && break
            ss_count +=1
            lnϕx, volx = modified_lnϕ(modelx, p, T, _x, Hϕx, phase=:liquid, vol0=volx)
            lnϕy, voly = modified_lnϕ(modely, p, T, _y, Hϕy, phase=:vapour, vol0=voly)
            if isnan(volx) || isnan(voly)
                break
            end

            if _bubble_or_lle
                _lnϕx = view(lnϕx,_view)
                lnK .=_lnϕx .- lnϕy
            else
                _lnϕy = view(lnϕy,_view)
                lnK .= lnϕx .- _lnϕy
            end

            K .= exp.(lnK)

            w_old .=  w

            if _bubble_or_lle
                __x = view(_x,_view)
                w .= __x .* K
                w_calc .= w
            else
                __y = view(_y,_view)
                w .= __y ./ K
                w_calc .= w
            end

            w ./= sum(w)
            error = dnorm(w,w_old,Inf) #||x-x_old||∞

            if _bubble_or_lle
                tpd_x = view(_x,_view)
                stability = dnorm(_y,tpd_x,Inf)
            else
                tpd_y = view(_y,_view)
                stability = dnorm(tpd_y,_x,Inf)
            end

            if stability < tol_stability #the interation procedure went wrong. perform a T/P movement first
                w .= w_restart
                volx,voly = volx_restart,voly_restart
                w .= w_restart
                valid_iter = false
                volx,voly = volx_restart,voly_restart
                if _bubble_or_lle
                    __x = view(_x,_view)
                    K .= y ./ __x
                else
                    __y = view(_y,_view)
                    K .= __y ./ x
                end
                lnK .= log.(K)
                break
            end
        end

        OF_old = OF
        OF = sum(w_calc) - 1.0

        if _pressure && second_order
            ∂lnϕ∂Px, volx = ∂lnϕ∂P(modelx, p, T, _x, Hϕx, phase=:liquid, vol0=volx)
            ∂lnϕ∂Py, voly = ∂lnϕ∂P(modely, p, T, _y, Hϕy, phase=:vapour, vol0=voly)
            if _bubble_or_lle
                _∂lnϕ∂Px = view(∂lnϕ∂Px, _view)
                ∂OF = @sum(w[i]*(_∂lnϕ∂Px[i] - ∂lnϕ∂Py[i]))
            else
                _∂lnϕ∂Py = view(∂lnϕ∂Py,_view)
                ∂OF = @sum(w[i]*(∂lnϕ∂Px[i] - _∂lnϕ∂Py[i]))
            end
        elseif !_pressure && second_order
            ∂lnϕ∂Tx, volx = ∂lnϕ∂T(modelx, p, T, _x, Hϕx, phase=:liquid, vol0=volx)
            ∂lnϕ∂Ty, voly = ∂lnϕ∂T(modely, p, T, _y, Hϕy, phase=:vapour, vol0=voly)
            if _bubble_or_lle
                _∂lnϕ∂Tx = view(∂lnϕ∂Tx,_view)
                ∂OF = @sum(w_calc[i]*(_∂lnϕ∂Tx[i] - ∂lnϕ∂Ty[i]))
            else
                _∂lnϕ∂Ty = view(∂lnϕ∂Ty,_view)
                ∂OF = @sum(w_calc[i]*(∂lnϕ∂Tx[i] - _∂lnϕ∂Ty[i]))
            end
        elseif _pressure && !second_order
            if j == 1
                ∂OF = OF/sqrt(eps(p))
            else
                ∂OF = (OF - OF_old)/(p - p_old)
            end
        else
            if j == 1
                ∂OF = OF/sqrt(eps(T))
            else
                ∂OF = (OF - OF_old)/(T - T_old)
            end
        end

        if isnan(volx) || isnan(voly)
            break
        end

        if !_bubble_or_lle && second_order
            ∂OF = -∂OF
        end

        ∂step = OF / ∂OF

        if _pressure
            ∂step∂p = clamp(∂step,-0.4*p,0.4*p)
            p_old = p
            p -= ∂step∂p
        else
            ∂step∂T = clamp(∂step,-0.05*T,0.05*T)
            T_old = T
            T -= ∂step∂T
            _update_temperature_with_view!(modelx,modely,T,_view)
        end

        propx =  _pressure ? p : T
        verbose && @info "$(__pad_val(j,13)) $(__pad_val(ss_count,12)) $(__pad_val(propx,16)) $(__pad_val(OF,16)) $(repr(K,context = :compact => true))"

        if valid_iter && (abs(∂step) < tol_pT || abs(OF) < tol_of)
            converged = true
            break
        end

        if !isfinite(∂step) #error, fail early, the NaN propagation is handled upstream
            converged = true
            break
        end
    end

    if !valid_iter
        w .= NaN
        lnK .= NaN
        volx,voly = w[1],w[1]
        if _pressure
            p = w[1]
        else
            T = w[1]
        end
    end

    if verbose && converged
        @info "algorithm converged in sucessive substitution phase."
    end
    return converged,(p,T,_x,_y,(volx,voly))
end

##general multidimensional non linear system generator to solve bubble/dew problems via fugacity coefficients
function _fug_J_∂i∂j!(J,w,∂lnϕ∂nw)
    neq = length(w)
    for i in 1:neq
        J[i,i] = 1
    end
    J[neq+1, neq+1] = 0
    Jwlnϕw = @view J[1:neq, 1:neq]
    Jw = @view(J[neq+1, 1:neq])
    ∂lnϕ∂nw .= w .* ∂lnϕ∂nw
    w_∂lnϕ∂nw = transpose(∂lnϕ∂nw)
    Jwlnϕw .+= w_∂lnϕ∂nw
    Jw .= w
    return J
end

##general multidimensional non linear system generator to solve bubble/dew problems via fugacity coefficients
function _fug_OF_neq!(F,J,inc,model,prop,z,data::FugData,cache)
    method = data.method
    _bubble_or_lle = FugEnum.is_bubble(method) || is_lle(method)
    _pressure = FugEnum.is_pressure(method)
    phasex,phasey = FugEnum.phases(method)
    _,K,w,u,_,_,vol_cache,Hϕx = cache
    second_order = !isnothing(J)
    volx,voly = vol_cache[]
    lnK = @view inc[1:end-1]
    K .= exp.(lnK)
    u .= z

    if _bubble_or_lle
        w .= K .* u
        x,y = u,w
    else
        w .= u ./ K
        x,y = w,u
    end

    propx = exp(last(inc))
    propy = oftype(propx,prop)
    if _pressure
        p,T = propx,propy
    else
        p,T = propy,propx
        update_temperature!(model,T)
    end

    if second_order
        J .= 0.0
        J1 = @view J[1:(end-1), end]
        if _pressure
            lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, volx = ∂lnϕ∂n∂P(model, p, T, x, Hϕx, phase=phasex, vol0=volx)
            J1 .-= p .* ∂lnϕ∂Px
            !_bubble_or_lle && _fug_J_∂i∂j!(J,x,∂lnϕ∂nx)
            !isnothing(F) && (F[1:end-1] .= lnK .- lnϕx)
            lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, voly = ∂lnϕ∂n∂P(model, p, T, y, Hϕx, phase=phasey, vol0=voly)
            J1 .+= p .* ∂lnϕ∂Py
            _bubble_or_lle && _fug_J_∂i∂j!(J,y,∂lnϕ∂ny)
        else
            lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, ∂lnϕ∂Tx, volx = ∂lnϕ∂n∂P∂T(model, p, T, x, Hϕx, phase=phasex, vol0=volx)
            !_bubble_or_lle && _fug_J_∂i∂j!(J,x,∂lnϕ∂nx)
            J1 .-= T .* ∂lnϕ∂Tx
            !isnothing(F) && (F[1:end-1] .= lnK .- lnϕx)
            lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, ∂lnϕ∂Ty, voly = ∂lnϕ∂n∂P∂T(model, p, T, y, Hϕx, phase=phasey, vol0=voly)
            J1 .+= T .* ∂lnϕ∂Ty
            _bubble_or_lle && _fug_J_∂i∂j!(J,y,∂lnϕ∂ny)
        end
        if !isnothing(F)
            Feq = @view F[1:end-1]
            Feq .+= lnϕy
            F[end] = sum(y)  - sum(x)
        end
    else
        Feq = @view F[1:end-1]
        lnϕx, volx = modified_lnϕ(model, p, T, x, Hϕx, phase=phasex, vol0=volx)
        Feq.= lnK .- lnϕx
        lnϕy, voly = modified_lnϕ(model, p, T, y, Hϕx, phase=phasey, vol0=voly)
        Feq .+= lnϕy
        F[end] = sum(y)  - sum(x)
    end
    vol_cache[] = (volx,voly)
    return nothing
end

##general multidimensional non linear system generator to solve bubble/dew problems via fugacity coefficients
##support for noncondensables/nonvolatiles
function _fug_OF_neq!(F,J,inc,modelx,modely,prop,z,_view,data::FugData,cache)
    method = data.method
    _bubble_or_lle = FugEnum.is_bubble(method) || is_lle(method)
    _pressure = FugEnum.is_pressure(method)
    phasex,phasey = FugEnum.phases(method)
    _,K,w,w2,fw1,fw2,vol_cache,Hϕx,Hϕy,u = cache
    second_order = !isnothing(J)
    volx,voly = vol_cache[]
    lnK = @view inc[1:end-1]
    K .= exp.(lnK)
    u .= z

    if _bubble_or_lle
        w .= K .* @view(u[_view])
        x,y = u,w
    else
        w .= @view(u[_view]) ./ K
        x,y = w,u
    end

    propx = exp(last(inc))
    propy = oftype(propx,prop)
    if _pressure
        p,T = propx,propy
    else
        p,T = propy,propx
        _update_temperature_with_view!(modelx,modely,T,_view)
    end

    if second_order
        J .= 0.0
        if _pressure
            lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, volx = ∂lnϕ∂n∂P(modelx, p, T, x, Hϕx, phase=phasex, vol0=volx)
            lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, voly = ∂lnϕ∂n∂P(modely, p, T, y, Hϕy, phase=phasey, vol0=voly)
            if _bubble_or_lle
                J[1:(end-1), end] .= p .* (∂lnϕ∂Py .- @view(∂lnϕ∂Px[_view]))
            else
                J[1:(end-1), end] .= p .* (@view(∂lnϕ∂Py[_view]) .- ∂lnϕ∂Px)
            end
        else
            lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, ∂lnϕ∂Tx, volx = ∂lnϕ∂n∂P∂T(modelx, p, T, x, Hϕx, phase=phasex, vol0=volx)
            lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, ∂lnϕ∂Ty, voly = ∂lnϕ∂n∂P∂T(modely, p, T, y, Hϕy, phase=phasey, vol0=voly)
            if _bubble_or_lle
                J[1:(end-1), end] .= T .* (∂lnϕ∂Ty .- @view(∂lnϕ∂Tx[_view]))
            else
                J[1:(end-1), end] .= T .* (@view(∂lnϕ∂Ty[_view]) .- ∂lnϕ∂Tx)
            end
        end
        if _bubble_or_lle
            _fug_J_∂i∂j!(J,y,∂lnϕ∂ny)
        else
            _fug_J_∂i∂j!(J,x,∂lnϕ∂nx)
        end
        if F !== nothing
            if _bubble_or_lle
                lnϕview = @view lnϕx[_view]
                F[1:end-1] = lnK .+ lnϕy .- lnϕview
            else
                lnϕview = @view lnϕy[_view]
                F[1:end-1] = lnK .+ lnϕview .- lnϕx
            end
            F[end] = sum(y)  - sum(x)
        end
    else
        lnϕx, volx = modified_lnϕ(modelx, p, T, x, Hϕx, phase=phasex, vol0=volx)
        lnϕy, voly = modified_lnϕ(modely, p, T, y, Hϕy, phase=phasey, vol0=voly)
        if _bubble_or_lle
            lnϕview = @view lnϕx[_view]
            F[1:end-1] = lnK .+ lnϕy .- lnϕview
        else
            lnϕview = @view lnϕy[_view]
            F[1:end-1] = lnK .+ lnϕview .- lnϕx
        end
        F[end] = sum(y)  - sum(x)
    end
    vol_cache[] = (volx,voly)
    return nothing
end

function _fug_OF_neq(model,prop,z,data,cache)
    function f!(F,x)
        _fug_OF_neq!(F,nothing,x,model,prop,z,data,cache)
        F
    end

    function fj!(F,J,x)
        _fug_OF_neq!(F,J,x,model,prop,z,data,cache)
        F,J
    end

    function j!(J, α)
        _fug_OF_neq!(nothing,J,x,model,prop,z,data,cache)
        return ∇f
    end

    obj = NLSolvers.VectorObjective(f!,j!,fj!,nothing)
    return NEqProblem(obj, inplace = true)
end

function _fug_OF_neq(modelx,modely,prop,z,_view,data,cache)
    function f!(F,x)
        _fug_OF_neq!(F,nothing,x,modelx,modely,prop,z,_view,data,cache)
        F
    end

    function fj!(F,J,x)
        _fug_OF_neq!(F,J,x,modelx,modely,prop,z,_view,data,cache)
        F,J
    end

    function j!(J, α)
        _fug_OF_neq!(nothing,J,x,modelx,modely,prop,z,_view,data,cache)
        return ∇f
    end

    obj = NLSolvers.VectorObjective(f!,j!,fj!,nothing)
    return NEqProblem(obj, inplace = true)
end

function rhovec_bubbledew_obj!(F, J, x, prop, model, z, in_equilibria, cache, calc_type::Val{IS_TEMP}) where IS_TEMP
    if (F === nothing && J === nothing)
        TT = eltype(cache.aux)
        return TT(NaN),TT(NaN)
    end

    nc = length(model)
    nx = count(in_equilibria)
    ip = nx + 1
    rhoz = x[ip]
    R = Rgas(model)
    #x[ip+1] = τ = p/(R*T), with prop = p
    τ = IS_TEMP ? x[ip + 1] : zero(rhoz)
    T = IS_TEMP ? prop/(R*τ) : prop
    dTdτ = IS_TEMP ? -T/τ : zero(T)
    RT = R*T
    if IS_TEMP && J !== nothing && length(cache.aux) < nc + 1
        throw(ArgumentError("a ∂lnϕ_cache with Val{true} is needed to differentiate with respect to T"))
    end

    hess_fn = IS_TEMP ? __μprT_hess : __μpr_hess
    J === nothing || fill!(J, zero(eltype(J)))
    ρ = @view cache.aux[1:nc]

    #reference phase: ρ = rhoz*z
    @. ρ = rhoz*z
    if J === nothing
        μz, pz_res = __μpr(model, T, ρ, cache)
    else
        μz, pz_res, Hz, dpdTz = hess_fn(model, T, ρ, cache)
    end
    ptot_z = sum(ρ) + pz_res

    k = 0
    az = zero(ptot_z)
    for i in 1:nc
        if J !== nothing
            hz = zero(az)
            for m in 1:nc
                hz += Hz[i, m]*z[m]
            end
            az += z[i]*(1 + rhoz*hz)
        end
        in_equilibria[i] || continue
        k += 1
        F === nothing || (F[k] = log(ρ[i]) + μz[i])
        if J !== nothing
            J[k, ip] = 1/rhoz + hz
            IS_TEMP && (J[k, ip + 1] = Hz[i, nc + 1]*dTdτ)
        end
    end
    if J !== nothing
        J[ip, ip] = az
        if IS_TEMP
            J[ip, ip + 1] = dpdTz*dTdτ
            J[ip + 1, ip] = az
            J[ip + 1, ip + 1] = dpdTz*dTdτ - 1
        end
    end

    #incipient phase: decompress x into ρ
    k = 0
    ρt_w = zero(eltype(ρ))
    for i in 1:nc
        if in_equilibria[i]
            k += 1
            ρ[i] = x[k]
            ρt_w += x[k]
        else
            ρ[i] = zero(eltype(ρ))
        end
    end
    if J === nothing
        μw, pw_res = __μpr(model, T, ρ, cache)
    else
        μw, pw_res, Hw, dpdTw = hess_fn(model, T, ρ, cache)
    end
    ptot_w = ρt_w + pw_res

    k = 0
    for i in 1:nc
        in_equilibria[i] || continue
        k += 1
        F === nothing || (F[k] -= log(ρ[i]) + μw[i])
        if J !== nothing
            ai = one(ptot_w)
            for m in 1:nc
                ai += Hw[i, m]*ρ[m]
            end
            J[ip, k] = -ai
            l = 0
            for j in 1:nc
                in_equilibria[j] || continue
                l += 1
                J[k, l] = -Hw[i, j]
            end
            J[k, k] -= 1/ρ[i]
            IS_TEMP && (J[k, ip + 1] -= Hw[i, nc + 1]*dTdτ)
        end
    end
    IS_TEMP && J !== nothing && (J[ip, ip + 1] -= dpdTw*dTdτ)

    #pressure rows, scaled by 1/rhoz (dimensionless)
    G = ptot_z - ptot_w
    if F !== nothing
        F[ip] = G/rhoz
        IS_TEMP && (F[ip + 1] = (ptot_z - τ)/rhoz)
    end
    if J !== nothing
        invρz = 1/rhoz
        for j in axes(J, 2)
            J[ip, j] *= invρz
            IS_TEMP && (J[ip + 1, j] *= invρz)
        end
        J[ip, ip] -= G*invρz^2
        IS_TEMP && (J[ip + 1, ip] -= (ptot_z - τ)*invρz^2)
    end
    return ptot_z, ptot_w
end

#applies the Newton step dx to x, respecting the bounds. dx is overwritten with the step actually applied.
#returns true if the step was truncated by any bound.
function _bubbledew_bounded_step!(x, dx, prop, model, z, in_equilibria, wbuf, ub_z, θ, calc_type::Val{IS_TEMP}) where IS_TEMP
    nx = count(in_equilibria)
    ip = nx + 1
    clipped = false
    rhoz_old = x[ip]
    ρtw_old = zero(eltype(x))
    for k in 1:nx
        ρtw_old += x[k]
    end

    #lower bounds: ρw, ρz and τ = p/(RT) are positive
    if any(i -> x[i] + dx[i] <= 0, eachindex(x))
        α = positive_linesearch(x, dx)
        dx .*= α
        clipped = true
    end
    x .+= dx

    T = IS_TEMP ? prop/(Rgas(model)*x[ip + 1]) : prop

    #upper bound on ρz: bracket between the previous iterate and the bound
    ub = IS_TEMP ? 1/lb_volume(model, T, z) : ub_z
    if x[ip] > θ*ub
        rhoz_new = min((rhoz_old + ub)/2, θ*ub)
        dx[ip] = rhoz_new - rhoz_old
        x[ip] = rhoz_new
        clipped = true
    end

    #upper bound on sum(ρw): bracket the total density, keeping the trial composition.
    k = 0
    ρtw_new = zero(eltype(x))
    for i in eachindex(wbuf)
        if in_equilibria[i]
            k += 1
            wbuf[i] = x[k]
            ρtw_new += x[k]
        else
            wbuf[i] = zero(eltype(wbuf))
        end
    end
    φ = lb_volume(model, T, wbuf)           #lb_volume of ρw at V = 1
    if φ > θ
        ρub_w = ρtw_new/φ                   #lb_volume(model,T,c*z) = c*lb_volume(model,T,z)
        ρtw_target = min((ρtw_old + ρub_w)/2, θ*ρub_w)
        c = ρtw_target/ρtw_new
        for k in 1:nx
            xnew = c*x[k]
            dx[k] += xnew - x[k]
            x[k] = xnew
        end
        clipped = true
    end
    return clipped
end

function bubbledew_rhovec(model,prop,z,rhoz0,rhovecw,prop_obj0,cache,in_equilibria,method,calc_type::Val{IS_TEMP}) where IS_TEMP
    atol, rtol, max_iters = method.atol, method.rtol, method.max_iters
    nc = length(model)
    TT = eltype(cache.aux)
    nc_actual = count(in_equilibria)
    nx = nc_actual + 1 + IS_TEMP
    F = similar(rhovecw,TT,nx)
    J = similar(rhovecw,TT,(nx,nx))
    x = similar(F)
    dx = similar(F)
    xnc = @view x[1:nc_actual]
    xnc .= @view rhovecw[in_equilibria]
    x[nc_actual + 1] = rhoz0
    IS_TEMP && (x[nc_actual + 2] = prop/(Rgas(model)*prop_obj0))   #τ0 = p/(R*T0)
    piv = Vector{Int}(undef,nx)

    wbuf = cache.scratch
    ub_z = IS_TEMP ? zero(TT) : 1/lb_volume(model,prop,z)
    θ = TT(0.999)
    clipped = false

    pz,pw = zero(TT),zero(TT)
    F0 = zero(TT)
    converged = false
    for iter in 0:max_iters
        pz,pw = rhovec_bubbledew_obj!(F,J,x,prop,model,z,in_equilibria,cache,calc_type)
        Fnorm = norm(F,Inf)
        isfinite(Fnorm) || break
        iter == 0 && (F0 = Fnorm)

        small_F = Fnorm <= max(atol, rtol*F0)
        small_dx = iter > 0 && !clipped && all(i -> abs(dx[i]) <= rtol*abs(x[i]), eachindex(x))
        if small_F || small_dx
            converged = true
            break
        end
        iter == max_iters && break

        Fac = Solvers.unsafe_LU!(J,piv)
        ldiv!(Fac,F)                    #destroys J, stores J\F in F
        dx .= .-F
        clipped = _bubbledew_bounded_step!(x,dx,prop,model,z,in_equilibria,wbuf,ub_z,θ,calc_type)
    end

    if !converged
        nan = TT(NaN)
        resize!(F,nc)
        F .= nan
        return nan,nan,nan,nan,F
    end

    rhoz = x[nc_actual + 1]
    if IS_TEMP
        T = prop/(Rgas(model)*x[nc_actual + 2])
        p = TT(prop)
    else
        T = TT(prop)
        p = sqrt(pz*pw)*Rgas(model)*T
    end

    vz = sum(z)/rhoz
    rhow = sum(xnc)
    vw = 1/rhow

    resize!(F,nc)
    F .= 0
    Fnc = @view F[in_equilibria]
    w = F
    Fnc .= xnc
    w ./= rhow
    return p,T,vz,vw,w
end

