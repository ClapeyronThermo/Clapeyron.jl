##general successive substitution method to solve bubble/dew problems via fugacity coefficients
module FugEnum
    @enum BubbleDew begin
        BUBBLE_PRESSURE
        BUBBLE_TEMPERATURE
        DEW_PRESSURE
        DEW_TEMPERATURE
    end
    is_bubble(x::BubbleDew) = (x == BUBBLE_PRESSURE || x == BUBBLE_TEMPERATURE)
    is_dew(x::BubbleDew) = !is_bubble(x)
    is_temperature(x::BubbleDew) = (x == BUBBLE_TEMPERATURE  || x == DEW_TEMPERATURE)
    is_pressure(x::BubbleDew) = !is_temperature(x) 
end

function _fug_OF_ss(model::EoSModel,p,T,x,y,vol0,method;itmax_ss = 5, itmax_newton = 10,tol_pT = 1e-8,tol_xy = 1e-8,tol_of = 1e-8)
    volx,voly = vol0
    converged = false
    _bubble,_pressure = FugEnum.is_bubble(method),FugEnum.is_pressure(method)
    tol_stability = abs2(cbrt(tol_xy))
    #caches for ∂lnϕ∂n∂P∂T/∂lnϕ∂n∂P
    if _pressure
        Hϕx = ∂lnϕ_cache(model, p, T, x, Val{false}())
        Hϕy = ∂lnϕ_cache(model, p, T, y, Val{false}())
    else
        Hϕx = ∂lnϕ_cache(model, p, T, x, Val{true}())
        Hϕy = ∂lnϕ_cache(model, p, T, y, Val{true}())
    end

    n = length(model)
    lnK = similar(Hϕx[3],n)
    K = similar(lnK)
    w = similar(lnK)
    w_old = similar(lnK)
    w_calc = similar(lnK)
    w_restart = similar(lnK)

    if _bubble
        w .= y
        _x,_y = x,w
    else
        w .= x
        _x,_y = w,y
    end

    valid_iter = true
    for j in 1:itmax_newton
        w_restart .= w_calc
        volx_restart,voly_restart = volx,voly
        valid_iter = true
        if isnan(volx) || isnan(voly)
            break
        end
        for i in 1:itmax_ss
            lnϕx, volx = lnϕ(model, p, T, _x, Hϕx, vol0=volx, phase = :liquid)
            lnϕy, voly = lnϕ(model, p, T, _y, Hϕy, vol0=voly, phase = :vapour)

            if isnan(volx)
                lnϕx, volx = lnϕ(model, p, T, _x, Hϕx, phase = :liquid)
            end

            if isnan(voly)
                lnϕy, voly = lnϕ(model, 1.1p, T, _y, Hϕy, phase = :vapour)
            end
            if isnan(volx) || isnan(voly)
                break
            end

            lnK .= lnϕx .- lnϕy
            K .= exp.(lnK)
            w_old .=  w

            if _bubble
                w .= _x .* K
                w_calc .= w
            else
                w .= _y ./ K
                w_calc .= w
            end
            w ./= sum(w)

            error = dnorm(w,w_old,Inf) #||x-x_old||∞

            if error < tol_xy
                break
            end
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
                break
            end
        end
        if _pressure
            lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, volx = ∂lnϕ∂n∂P(model, p, T, _x, Hϕx, phase=:liquid, vol0=volx)
            lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, voly = ∂lnϕ∂n∂P(model, p, T, _y, Hϕy, phase=:vapour, vol0=voly)
            ∂OF =@sum(w_calc[i]*(∂lnϕ∂Px[i] - ∂lnϕ∂Py[i]))
        else
            lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, ∂lnϕ∂Tx, volx = ∂lnϕ∂n∂P∂T(model, p, T, _x, Hϕx, phase=:liquid, vol0=volx)
            lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, ∂lnϕ∂Ty, voly = ∂lnϕ∂n∂P∂T(model, p, T, _y, Hϕy, phase=:vapour, vol0=voly)
            ∂OF = @sum(w_calc[i]*(∂lnϕ∂Tx[i] - ∂lnϕ∂Ty[i]))
        end
        if isnan(volx) || isnan(voly)
            break
        end

        if !_bubble
            ∂OF = -∂OF
        end

        lnK .= lnϕx .- lnϕy
        K .= exp.(lnK)

        OF = sum(w_calc) - 1.0
        ∂step = OF / ∂OF
        if valid_iter && abs(∂step) < tol_pT || abs(OF) < tol_of
            converged = true
            break
        end



        if _pressure
            ∂step = clamp(∂step,-0.4*p,0.4*p)
            p -= ∂step
        else
            ∂step = clamp(∂step,-0.05*T,0.05*T)
            T -= ∂step
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

    return converged,(p,T,_x,_y,(volx,voly),lnK)
end

function _fug_OF_ss(modelx::EoSModel,modely::Nothing,p,T,x,y,vol0,method,_view;itmax_ss = 5, itmax_newton = 10,tol_pT = 1e-8,tol_xy = 1e-8,tol_of = 1e-8)
    return _fug_OF_ss(modelx,p,T,x,y,vol0,method;itmax_ss, itmax_newton,tol_pT,tol_xy,tol_of)
end

function _fug_OF_ss(modelx::Nothing,modely::EoSModel,p,T,x,y,vol0,method,_view;itmax_ss = 5, itmax_newton = 10,tol_pT = 1e-8,tol_xy = 1e-8,tol_of = 1e-8)
    return _fug_OF_ss(modely,p,T,x,y,vol0,method;itmax_ss, itmax_newton,tol_pT,tol_xy,tol_of)
end


##general successive substitution method to solve bubble/dew problems via fugacity coefficients
##support for nonvolatiles/noncondensables

function _fug_OF_ss(modelx::EoSModel,modely::EoSModel,p,T,x,y,vol0,method,_view;itmax_ss = 5, itmax_newton = 10,tol_pT = 1e-8,tol_xy = 1e-8,tol_of = 1e-8)
    volx,voly = vol0
    converged = false
    _bubble,_pressure = FugEnum.is_bubble(method),FugEnum.is_pressure(method)
    tol_stability = abs2(cbrt(tol_xy))
    if _bubble
        n = length(modely)
        zz = view(x,_view)
    else
        n = length(modelx)
        zz = view(y,_view)
    end


    #caches for ∂lnϕ∂n∂P∂T/∂lnϕ∂n∂P
    if _pressure
        Hϕx = ∂lnϕ_cache(modelx, p, T, x, Val{false}())
        Hϕy = ∂lnϕ_cache(modely, p, T, y, Val{false}())
    else
        Hϕx = ∂lnϕ_cache(modelx, p, T, x, Val{true}())
        Hϕy = ∂lnϕ_cache(modely, p, T, y, Val{true}())
    end

    lnK = similar(Hϕx[3],n)
    K = similar(lnK)
    w = similar(lnK)
    w_old = similar(lnK)
    w_calc = similar(lnK)
    w_restart = similar(lnK)

    if _bubble
        w .= y
        _x,_y = x,w
    else
        w .= x
        _x,_y = w,y
    end
    valid_iter = false
    for j in 1:itmax_newton
        if isnan(volx) || isnan(voly)
            break
        end
        w_restart .= w_calc
        volx_restart,voly_restart = volx,voly
        valid_iter = true

        for i in 1:itmax_ss

            lnϕx, volx = lnϕ(modelx, p, T, _x, Hϕx, phase=:liquid, vol0=volx)
            lnϕy, voly = lnϕ(modely, p, T, _y, Hϕy, phase=:vapour, vol0=voly)

            if isnan(volx)
                lnϕx, volx = lnϕ(modelx, p, T, _x, Hϕx, phase = :liquid)
            end

            if isnan(voly)
                lnϕy, voly = lnϕ(modely, 1.1p, T, _y, Hϕy, phase = :vapour)
            end

            if isnan(volx) || isnan(voly)
                break
            end

            if _bubble
                _lnϕx = view(lnϕx,_view)
                lnK .=_lnϕx .- lnϕy
            else
                _lnϕy = view(lnϕy,_view)
                lnK .= lnϕx .- _lnϕy
            end

            K .= exp.(lnK)
            w_old .=  w

            if _bubble
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
            if error < tol_xy
                break
            end

            if _bubble
                #tpd_lnϕx = view(lnϕx,_view)
                tpd_x = view(_x,_view)
                stability = dnorm(_y,tpd_x,Inf)
                #tpd = 1 + @sum(_y[i]*(lnϕy[i] + log(y[i]) - log(tpd_x[i]) - tpd_lnϕx[i] - 1))
            else
                #tpd_lnϕy = view(lnϕy,_view)
                tpd_y = view(_y,_view)
                stability = dnorm(tpd_y,_x,Inf)
                #tpd = 1 + @sum(_x[i]*(lnϕx[i] + log(x[i]) - log(tpd_y[i]) - tpd_lnϕy[i] - 1))
            end
            if stability < tol_stability #the interation procedure went wrong. perform a T/P movement first
                w .= w_restart
                volx,voly = volx_restart,voly_restart
                break
            end
        end

        if _pressure
            lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, volx = ∂lnϕ∂n∂P(modelx, p, T, _x, Hϕx, phase=:liquid, vol0=volx)
            lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, voly = ∂lnϕ∂n∂P(modely, p, T, _y, Hϕy, phase=:vapour, vol0=voly)
            if _bubble
                _∂lnϕ∂Px = view(∂lnϕ∂Px, _view)
                ∂OF = @sum(w[i]*(_∂lnϕ∂Px[i] - ∂lnϕ∂Py[i]))
            else
                _∂lnϕ∂Py = view(∂lnϕ∂Py,_view)
                ∂OF = @sum(w[i]*(∂lnϕ∂Px[i] - _∂lnϕ∂Py[i]))
            end
        else
            lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, ∂lnϕ∂Tx, volx = ∂lnϕ∂n∂P∂T(modelx, p, T, _x, Hϕx, phase=:liquid, vol0=volx)
            lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, ∂lnϕ∂Ty, voly = ∂lnϕ∂n∂P∂T(modely, p, T, _y, Hϕy, phase=:vapour, vol0=voly)
            if _bubble
                _∂lnϕ∂Tx = view(∂lnϕ∂Tx,_view)
                ∂OF = @sum(w_calc[i]*(_∂lnϕ∂Tx[i] - ∂lnϕ∂Ty[i]))
            else
                _∂lnϕ∂Ty = view(∂lnϕ∂Ty,_view)
                ∂OF = @sum(w_calc[i]*(∂lnϕ∂Tx[i] - _∂lnϕ∂Ty[i]))
            end
        end
        if isnan(volx) || isnan(voly)
            break
        end

        if !_bubble
            ∂OF = -∂OF
        end

        if _bubble
            _lnϕx = view(lnϕx,_view)
            lnK .=_lnϕx .- lnϕy
        else
            _lnϕy = view(lnϕy,_view)
            lnK .= lnϕx .- _lnϕy
        end
        K .= exp.(lnK)

        OF = sum(w_calc) - 1.0
        ∂step = OF / ∂OF

        if _pressure
            ∂step = clamp(∂step,-0.4*p,0.4*p)
            p -= ∂step
        else
            ∂step = clamp(∂step,-0.05*T,0.05*T)
            T -= ∂step
        end

        if valid_iter && abs(∂step) < tol_pT || abs(OF) < tol_of
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

    return converged,(p,T,_x,_y,(volx,voly),lnK)
end

##general multidimensional non linear system generator to solve bubble/dew problems via fugacity coefficients

function _select_xy(w,K,x,y,_bubble)
    if _bubble
        w .= K .* x
        return x, w
    else
        w .= y ./ K
        return w, y
    end
end

function _select_pT(inc,p,T,_pressure)
    if _pressure
        return exp(inc[end]),T
    else
        return p,exp(inc[end])
    end
end

function _fug_J(J,x,y,∂lnϕ∂nx,∂lnϕ∂ny,_bubble)
    for i in 1:size(J,1)
        J[i,i] = 1
    end
    J[end, end] = 0
    Jwlnϕw = @view J[1:(end-1), 1:(end-1)]
    Jw = @view(J[end, 1:(end-1)])
    if _bubble
        ∂lnϕ∂ny .= y .* ∂lnϕ∂ny
        w = y
        w_∂lnϕ∂nw = transpose(∂lnϕ∂ny)
    else
        ∂lnϕ∂nx .= x .* ∂lnϕ∂nx
        w_∂lnϕ∂nw = transpose(∂lnϕ∂nx)
        w = x
    end
    Jwlnϕw .+= w_∂lnϕ∂nw
    Jw .= w
    return J
end

function _fug_OF_neqsystem(model,_x, _y, _p, _T, vol_cache,method,_phase)
    
    _bubble,_pressure = FugEnum.is_bubble(method),FugEnum.is_pressure(method)
 
    #caches for ∂lnϕ∂n∂P∂T/∂lnϕ∂n∂P

    XX = something(_p,_T)
    _w = something(_x,_y)

    if _pressure
        Hϕx = ∂lnϕ_cache(model, XX, XX, _w, Val{false}())
        Hϕy = ∂lnϕ_cache(model, XX, XX, _w, Val{false}())
    else
        Hϕx = ∂lnϕ_cache(model, XX, XX, _w, Val{true}())
        Hϕy = ∂lnϕ_cache(model, XX, XX, _w, Val{true}())
    end

    K = similar(_w)
    w = similar(_w)

    function f!(F, inc)
        volx, voly = vol_cache
        lnK = @view inc[1:end-1]
        K .= exp.(lnK)
        p,T = _select_pT(inc,_p,_T,_pressure)
        x,y = _select_xy(w,K,_x,_y,_bubble)
        lnϕx, volx = lnϕ(model, p, T, x, Hϕx, phase=_phase[1], vol0=volx)
        lnϕy, voly = lnϕ(model, p, T, y, Hϕy, phase=_phase[2], vol0=voly)
        F[1:end-1] .= lnK .+ lnϕy .- lnϕx
        F[end] = sum(y)  - sum(x)
        vol_cache .= (volx, voly)
        return F
    end

    function fj!(F,J,inc)
        volx, voly = vol_cache
        lnK = @view inc[1:end-1]
        K .= exp.(lnK)
        p,T = _select_pT(inc,_p,_T,_pressure)
        x,y = _select_xy(w,K,_x,_y,_bubble)
        J .= 0.0
        if _pressure
            lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, volx = ∂lnϕ∂n∂P(model, p, T, x, Hϕx, phase=_phase[1], vol0=volx)
            lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, voly = ∂lnϕ∂n∂P(model, p, T, y, Hϕy, phase=_phase[2], vol0=voly)
            J[1:(end-1), end] .= p .* (∂lnϕ∂Py .- ∂lnϕ∂Px)
        else
            lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, ∂lnϕ∂Tx, volx = ∂lnϕ∂n∂P∂T(model, p, T, x, Hϕx, phase=_phase[1], vol0=volx)
            lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, ∂lnϕ∂Ty, voly = ∂lnϕ∂n∂P∂T(model, p, T, y, Hϕy, phase=_phase[2], vol0=voly)
            J[1:(end-1), end] .= T .* (∂lnϕ∂Ty .- ∂lnϕ∂Tx)
        end

        F[1:end-1] .= lnK .+ lnϕy .- lnϕx
        F[end] = sum(y)  - sum(x)
        _fug_J(J,x,y,∂lnϕ∂nx,∂lnϕ∂ny,_bubble)
        vol_cache .= (volx, voly)
        return F,J
    end

    function j!(J,inc)
        volx, voly = vol_cache
        lnK = @view inc[1:end-1]
        K .= exp.(lnK)
        p,T = _select_pT(inc,_p,_T,_pressure)
        x,y = _select_xy(w,K,_x,_y,_bubble)
        J .= 0.0
        if _pressure
            lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, volx = ∂lnϕ∂n∂P(model, p, T, x, phase=_phase[1], vol0=volx)
            lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, voly = ∂lnϕ∂n∂P(model, p, T, y, phase=_phase[2], vol0=voly)
            J[1:(end-1), end] .= p .* (∂lnϕ∂Py .- ∂lnϕ∂Px)
        else
            lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, ∂lnϕ∂Tx, volx = ∂lnϕ∂n∂P∂T(model, p, T, x, phase=_phase[1], vol0=volx)
            lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, ∂lnϕ∂Ty, voly = ∂lnϕ∂n∂P∂T(model, p, T, y, phase=_phase[2], vol0=voly)
            J[1:(end-1), end] .= T .* (∂lnϕ∂Ty .- ∂lnϕ∂Tx)
        end
        _fug_J(J,x,y,∂lnϕ∂nx,∂lnϕ∂ny,_bubble)
        vol_cache .= (volx, voly)
        return J
    end

    return Solvers.NLSolvers.VectorObjective(f!,j!,fj!,nothing) |> Solvers.NLSolvers.NEqProblem
end

##general multidimensional non linear system generator to solve bubble/dew problems via fugacity coefficients
##support for noncondensables/nonvolatiles

#dispatch to simplified function
function _fug_OF_neqsystem(modelx::EoSModel,modely::Nothing,_x, _y, _p, _T, vol_cache,method,_phase,_view)
    return  _fug_OF_neqsystem(modelx,_x, _y, _p, _T, vol_cache,method,_phase)
end

function _fug_OF_neqsystem(modelx::Nothing,modely::EoSModel,_x, _y, _p, _T, vol_cache,method,_phase,_view)
    return  _fug_OF_neqsystem(modely,_x, _y, _p, _T, vol_cache,method,_phase)
end


#support for views
function _select_xy(w,K,x,y,_bubble,_view)
    if _bubble
        xv = @view x[_view]
        w .= K .* xv
        return x, w
    else
        yv = @view y[_view]
        w .= yv ./ K
        return w , y
    end
end

#main function
function _fug_OF_neqsystem(modelx::EoSModel,modely::EoSModel,_x, _y, _p, _T, vol_cache,method,_phase,_view)
    _bubble,_pressure = FugEnum.is_bubble(method),FugEnum.is_pressure(method)
    if _bubble
        w = similar(_y)
        wx = x
        wy = w
    else
        w = similar(_x)
        wx = w
        wy = y
    end

    K = similar(w)
    XX = something(_p,_T)
    if _pressure
        Hϕx = ∂lnϕ_cache(model, XX, XX, wx, Val{false}())
        Hϕy = ∂lnϕ_cache(model, XX, XX, wy, Val{false}())
    else
        Hϕx = ∂lnϕ_cache(model, XX, XX, wx, Val{true}())
        Hϕy = ∂lnϕ_cache(model, XX, XX, wy, Val{true}())
    end

    function f!(F, inc)
        volx, voly = vol_cache
        lnK = @view inc[1:end-1]
        K .= exp.(lnK)
        p,T = _select_pT(inc,_p,_T,_pressure)
        x,y = _select_xy(w,K,_x,_y,_bubble,_view)

        lnϕx, volx = lnϕ(modelx, p, T, x, Hϕx, phase=_phase[1], vol0=volx)
        lnϕy, voly = lnϕ(modely, p, T, y, Hϕy, phase=_phase[2], vol0=voly)

        if _bubble
            lnϕview = @view lnϕx[_view]
            F[1:end-1] .= lnK .+ lnϕy .- lnϕview
        else
            lnϕview = @view lnϕy[_view]
            F[1:end-1] .= lnK .+ lnϕview .- lnϕx
        end
        F[end] = sum(y)  - sum(x)
        vol_cache .= (volx, voly)
        return F
    end

    function fj!(F,J,inc)
        volx, voly = vol_cache
        lnK = @view inc[1:end-1]
        K .= exp.(lnK)
        p,T = _select_pT(inc,_p,_T,_pressure)
        x,y = _select_xy(w,K,_x,_y,_bubble,_view)
        J .= 0.0
        if _pressure
            lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, volx = ∂lnϕ∂n∂P(modelx, p, T, x, Hϕx, phase=_phase[1], vol0=volx)
            lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, voly = ∂lnϕ∂n∂P(modely, p, T, y, Hϕy, phase=_phase[2], vol0=voly)
            if _bubble
                J[1:(end-1), end] .= p .* (∂lnϕ∂Py .- @view(∂lnϕ∂Px[_view]))
            else
                J[1:(end-1), end] .= p .* (@view(∂lnϕ∂Py[_view]) .- ∂lnϕ∂Px)
            end
        else
            lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, ∂lnϕ∂Tx, volx = ∂lnϕ∂n∂P∂T(modelx, p, T, x, Hϕx, phase=_phase[1], vol0=volx)
            lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, ∂lnϕ∂Ty, voly = ∂lnϕ∂n∂P∂T(modely, p, T, y, Hϕy, phase=_phase[2], vol0=voly)
            if _bubble
                J[1:(end-1), end] .= T .* (∂lnϕ∂Ty .- @view(∂lnϕ∂Tx[_view]))
            else
                J[1:(end-1), end] .= T .* (@view(∂lnϕ∂Ty[_view]) .- ∂lnϕ∂Tx)
            end
        end

        if _bubble
            lnϕview = @view lnϕx[_view]
            F[1:end-1] = lnK .+ lnϕy .- lnϕview
        else
            lnϕview = @view lnϕy[_view]
            F[1:end-1] = lnK .+ lnϕview .- lnϕx
        end
        F[end] = sum(y)  - sum(x)
        _fug_J(J,x,y,∂lnϕ∂nx,∂lnϕ∂ny,_bubble)
        vol_cache .= (volx, voly)
        return F,J
    end

    function j!(J,inc)
        volx, voly = vol_cache
        lnK = @view inc[1:end-1]
        K .= exp.(lnK)
        p,T = _select_pT(inc,_p,_T,_pressure)
        x,y = _select_xy(K,_x,_y,_bubble,_view)
        J .= 0.0
        if _pressure
            lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, volx = ∂lnϕ∂n∂P(modelx, p, T, x, Hϕx, phase=_phase[1], vol0=volx)
            lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, voly = ∂lnϕ∂n∂P(modely, p, T, y, Hϕy, phase=_phase[2], vol0=voly)
            if _bubble
                J[1:(end-1), end] .= p .* (∂lnϕ∂Py .- @view(∂lnϕ∂Px[_view]))
            else
                J[1:(end-1), end] .= p .* (@view(∂lnϕ∂Py[_view]) .- ∂lnϕ∂Px)
            end
        else
            lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, ∂lnϕ∂Tx, volx = ∂lnϕ∂n∂P∂T(modelx, p, T, x, Hϕx, phase=_phase[1], vol0=volx)
            lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, ∂lnϕ∂Ty, voly = ∂lnϕ∂n∂P∂T(modely, p, T, y, Hϕy, phase=_phase[2], vol0=voly)
            if _bubble
                J[1:(end-1), end] .= T .* (∂lnϕ∂Ty .- @view(∂lnϕ∂Tx[_view]))
            else
                J[1:(end-1), end] .= T .* (@view(∂lnϕ∂Ty[_view]) .- ∂lnϕ∂Tx)
            end
        end

        _fug_J(J,x,y,∂lnϕ∂nx,∂lnϕ∂ny,_bubble)
        vol_cache .= (volx, voly)
        return J
    end

    return Solvers.NLSolvers.VectorObjective(f!,j!,fj!,nothing) |> Solvers.NLSolvers.NEqProblem
end
