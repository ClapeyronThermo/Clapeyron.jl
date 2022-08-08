##general successive substitution method to solve bubble/dew problems via fugacity coefficients

function _fug_OF_ss(model::EoSModel,p,T,x,y,vol0,_bubble,_pressure;itmax_ss = 5, itmax_newton = 10,tol_pT = 1e-8,tol_xy = 1e-8,tol_of = 1e-8)
    volx,voly = vol0
    converged = false

    lnϕx, volx = lnϕ(model, p, T, x, phase=:liquid, vol0=volx)
    lnϕy, voly = lnϕ(model, p, T, y, phase=:vapor, vol0=voly)

    n = length(model)
    lnK = zeros(n)
    K = zeros(n)
    w = zeros(n)
    w_old = zeros(n)
    w_calc = zeros(n)

    if _bubble
        w .= y
        _x,_y = x,w
    else
        w .= x
        _x,_y = w,y
    end

    for j in 1:itmax_newton
        lnϕx, volx = lnϕ(model, p, T, _x, phase=:liquid, vol0=volx)
        lnϕy, voly = lnϕ(model, p, T, _y, phase=:vapor, vol0=voly)

        for i in 1:itmax_ss
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

            lnϕx, volx = lnϕ(model, p, T, _x, phase=:liquid, vol0=volx)
            lnϕy, voly = lnϕ(model, p, T, _y, phase=:vapor, vol0=voly)

        end

        if _pressure
            lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, volx = ∂lnϕ∂n∂P(model, p, T, _x, phase=:liquid, vol0=volx)
            lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, voly = ∂lnϕ∂n∂P(model, p, T, _y, phase=:vapor, vol0=voly)
            ∂OF = sum(w.*(∂lnϕ∂Px .- ∂lnϕ∂Py))
        else
            lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, ∂lnϕ∂Tx, volx = ∂lnϕ∂n∂P∂T(model, p, T, _x, phase=:liquid, vol0=volx)
            lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, ∂lnϕ∂Ty, voly = ∂lnϕ∂n∂P∂T(model, p, T, _y, phase=:vapor, vol0=voly)
            ∂OF = sum(w.*(∂lnϕ∂Tx .- ∂lnϕ∂Ty))
        end

        if !_bubble
            ∂OF = -∂OF
        end

        lnK .= lnϕx .- lnϕy
        K .= exp.(lnK)

        OF = sum(w_calc) - 1.0
        ∂step = OF / ∂OF

        if abs(∂step) < tol_pT || abs(OF) < tol_of
            converged = true
            break
        end

        if !isfinite(∂step) #error, fail early, the NaN propagation is handled upstream
            converged = true
            break
        end

        if _pressure
            ∂step = min(∂step,0.4*p)
            p -= ∂step
        else
            ∂step = min(∂step,0.2*T)
            T -= ∂step
        end
    end

    return converged,(p,T,_x,_y,(volx,voly),lnK)
end

function _fug_OF_ss(modelx::EoSModel,modely::Nothing,p,T,x,y,vol0,_bubble,_pressure,_view;itmax_ss = 5, itmax_newton = 10,tol_pT = 1e-8,tol_xy = 1e-8,tol_of = 1e-8)
    return _fug_OF_ss(modelx,p,T,x,y,vol0,_bubble,_pressure;itmax_ss, itmax_newton,tol_pT,tol_xy,tol_of)
end

function _fug_OF_ss(modelx::Nothing,modely::EoSModel,p,T,x,y,vol0,_bubble,_pressure,_view;itmax_ss = 5, itmax_newton = 10,tol_pT = 1e-8,tol_xy = 1e-8,tol_of = 1e-8)
    return _fug_OF_ss(modely,p,T,x,y,vol0,_bubble,_pressure;itmax_ss, itmax_newton,tol_pT,tol_xy,tol_of)
end


##general successive substitution method to solve bubble/dew problems via fugacity coefficients
##support for nonvolatiles/noncondensables

function _fug_OF_ss(modelx::EoSModel,modely::EoSModel,p,T,x,y,vol0,_bubble,_pressure,_view;itmax_ss = 5, itmax_newton = 10,tol_pT = 1e-8,tol_xy = 1e-8,tol_of = 1e-8)
    volx,voly = vol0
    converged = false

    lnϕx, volx = lnϕ(modelx, p, T, x, phase=:liquid, vol0=volx)
    lnϕy, voly = lnϕ(modely, p, T, y, phase=:vapor, vol0=voly)
    if _bubble
        n = length(modely)
    else
        n = length(modelx)
    end

    lnK = zeros(n)
    K = zeros(n)
    w = zeros(n)
    w_old = zeros(n)
    w_calc = zeros(n)

    if _bubble
        w .= y
        _x,_y = x,w
    else
        w .= x
        _x,_y = w,y
    end
    for j in 1:itmax_newton
        lnϕx, volx = lnϕ(modelx, p, T, _x, phase=:liquid, vol0=volx)
        lnϕy, voly = lnϕ(modely, p, T, _y, phase=:vapor, vol0=voly)
        for i in 1:itmax_ss
            if _bubble
                lnK .= lnϕx[_view] .- lnϕy
            else
                lnK .= lnϕx .- lnϕy[_view]
            end
            
            K .= exp.(lnK)
            w_old .=  w

            if _bubble
                w .= _x[_view] .* K
                w_calc .= w
            else
                w .= _y[_view] ./ K
                w_calc .= w
            end
            w ./= sum(w)
            error = dnorm(w,w_old,Inf) #||x-x_old||∞
            if error < tol_xy
                break
            end
            lnϕx, volx = lnϕ(modelx, p, T, _x, phase=:liquid, vol0=volx)
            lnϕy, voly = lnϕ(modely, p, T, _y, phase=:vapor, vol0=voly)

        end
       
        if _pressure
            lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, volx = ∂lnϕ∂n∂P(modelx, p, T, _x, phase=:liquid, vol0=volx)
            lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, voly = ∂lnϕ∂n∂P(modely, p, T, _y, phase=:vapor, vol0=voly)
            if _bubble
                ∂OF = sum(w.*(∂lnϕ∂Px[_view] .- ∂lnϕ∂Py))
            else
                ∂OF = sum(w.*(∂lnϕ∂Px .- ∂lnϕ∂Py[_view]))
            end
        else
            lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, ∂lnϕ∂Tx, volx = ∂lnϕ∂n∂P∂T(modelx, p, T, _x, phase=:liquid, vol0=volx)
            lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, ∂lnϕ∂Ty, voly = ∂lnϕ∂n∂P∂T(modely, p, T, _y, phase=:vapor, vol0=voly)
            if _bubble
                ∂OF = sum(w.*(∂lnϕ∂Tx[_view] .- ∂lnϕ∂Ty))
            else
                ∂OF = sum(w.*(∂lnϕ∂Tx .- ∂lnϕ∂Ty[_view]))
            end
        end

        if !_bubble
            ∂OF = -∂OF
        end

        if _bubble
            lnK .= lnϕx[_view] .- lnϕy
        else
            lnK .= lnϕx .- lnϕy[_view]
        end
        K .= exp.(lnK)

        OF = sum(w_calc) - 1.0
        ∂step = OF / ∂OF

        if _pressure
            #∂step = min(∂step,0.4*p)
            if ∂step > p
                ∂step = 0.4*p
            end
            p -= ∂step
        else
            if ∂step > T
                ∂step = 0.2*T
            end
            T -= ∂step
        end
        if abs(∂step) < tol_pT || abs(OF) < tol_of
            converged = true
            break
        end
        
        if !isfinite(∂step) #error, fail early, the NaN propagation is handled upstream
            converged = true
            break
        end
    end

    return converged,(p,T,_x,_y,(volx,voly),lnK)
end

##general multidimensional non linear system generator to solve bubble/dew problems via fugacity coefficients

function _select_xy(K,x,y,_bubble)
    if _bubble
        return x, K .* x
    else
        return y ./ K, y
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
    if _bubble
        J[diagind(J)] .= 1.
        J[1:(end-1), 1:(end-1)] += (y .* ∂lnϕ∂ny)'
        J[end, 1:(end-1)] = y
        J[end, end] = 0.
    else
        J[diagind(J)] .= 1.
        J[1:(end-1), 1:(end-1)] += (x .* ∂lnϕ∂nx)'
        J[end, 1:(end-1)] = x
        J[end, end] = 0.
    end
    return J
end

function _fug_OF_neqsystem(model,_x, _y, _p, _T, vol_cache,_bubble,_pressure,_phase)
    function f!(F, inc)
        volx, voly = vol_cache
        lnK = inc[1:end-1]
        K = exp.(lnK)
        p,T = _select_pT(inc,_p,_T,_pressure)
        x,y = _select_xy(K,_x,_y,_bubble)
        lnϕx, volx = lnϕ(model, p, T, x, phase=_phase[1], vol0=volx)
        lnϕy, voly = lnϕ(model, p, T, y, phase=_phase[2], vol0=voly)

        F[1:end-1] = lnK .+ lnϕy .- lnϕx
        F[end] = sum(y)  - sum(x)

        vol_cache .= (volx, voly)
        return F
    end

    function fj!(F,J,inc)
        volx, voly = vol_cache
        lnK = inc[1:end-1]
        K = exp.(lnK)
        p,T = _select_pT(inc,_p,_T,_pressure)
        x,y = _select_xy(K,_x,_y,_bubble)
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

        F[1:end-1] = lnK .+ lnϕy .- lnϕx
        F[end] = sum(y)  - sum(x)
        _fug_J(J,x,y,∂lnϕ∂nx,∂lnϕ∂ny,_bubble)
        vol_cache .= (volx, voly)
        return F,J
    end

    function j!(J,inc)
        volx, voly = vol_cache
        lnK = inc[1:end-1]
        K = exp.(lnK)
        p,T = _select_pT(inc,_p,_T,_pressure)
        x,y = _select_xy(K,_x,_y,_bubble)
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

    function jv!(inc)
        return nothing
    end

    return Solvers.NLSolvers.VectorObjective(f!,j!,fj!,jv!) |> Solvers.NLSolvers.NEqProblem
end

##general multidimensional non linear system generator to solve bubble/dew problems via fugacity coefficients
##support for noncondensables/nonvolatiles

#dispatch to simplified function
function _fug_OF_neqsystem(modelx::EoSModel,modely::Nothing,_x, _y, _p, _T, vol_cache,_bubble,_pressure,_phase,_view)
    return  _fug_OF_neqsystem(modelx,_x, _y, _p, _T, vol_cache,_bubble,_pressure,_phase)
end

function _fug_OF_neqsystem(modelx::Nothing,modely::EoSModel,_x, _y, _p, _T, vol_cache,_bubble,_pressure,_phase,_view)
    return  _fug_OF_neqsystem(modely,_x, _y, _p, _T, vol_cache,_bubble,_pressure,_phase)
end


#support for views
function _select_xy(K,x,y,_bubble,_view)
    if _bubble
        return x, K .* x[_view]
    else
        return y[_view] ./ K , y
    end
end

#main function
function _fug_OF_neqsystem(modelx::EoSModel,modely::EoSModel,_x, _y, _p, _T, vol_cache,_bubble,_pressure,_phase,_view)
    function f!(F, inc)
        volx, voly = vol_cache
        lnK = inc[1:end-1]
        K = exp.(lnK)
        p,T = _select_pT(inc,_p,_T,_pressure)
        x,y = _select_xy(K,_x,_y,_bubble,_view)

        lnϕx, volx = lnϕ(modelx, p, T, x, phase=_phase[1], vol0=volx)
        lnϕy, voly = lnϕ(modely, p, T, y, phase=_phase[2], vol0=voly)
        
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
        lnK = inc[1:end-1]
        K = exp.(lnK)
        p,T = _select_pT(inc,_p,_T,_pressure)
        x,y = _select_xy(K,_x,_y,_bubble,_view)  
        J .= 0.0
        if _pressure
            lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, volx = ∂lnϕ∂n∂P(modelx, p, T, x, phase=_phase[1], vol0=volx)
            lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, voly = ∂lnϕ∂n∂P(modely, p, T, y, phase=_phase[2], vol0=voly)
            if _bubble
                J[1:(end-1), end] .= p .* (∂lnϕ∂Py .- @view(∂lnϕ∂Px[_view]))
            else
                J[1:(end-1), end] .= p .* (@view(∂lnϕ∂Py[_view]) .- ∂lnϕ∂Px)
            end
        else
            lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, ∂lnϕ∂Tx, volx = ∂lnϕ∂n∂P∂T(modelx, p, T, x, phase=_phase[1], vol0=volx)
            lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, ∂lnϕ∂Ty, voly = ∂lnϕ∂n∂P∂T(modely, p, T, y, phase=_phase[2], vol0=voly)
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
        lnK = inc[1:end-1]
        K = exp.(lnK)
        p,T = _select_pT(inc,_p,_T,_pressure)
        x,y = _select_xy(K,_x,_y,_bubble,_view)
        J .= 0.0
        if _pressure
            lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, volx = ∂lnϕ∂n∂P(modelx, p, T, x, phase=_phase[1], vol0=volx)
            lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, voly = ∂lnϕ∂n∂P(modely, p, T, y, phase=_phase[2], vol0=voly)
            if _bubble
                J[1:(end-1), end] .= p .* (∂lnϕ∂Py .- @view(∂lnϕ∂Px[_view]))
            else
                J[1:(end-1), end] .= p .* (@view(∂lnϕ∂Py[_view]) .- ∂lnϕ∂Px)
            end
        else
            lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, ∂lnϕ∂Tx, volx = ∂lnϕ∂n∂P∂T(modelx, p, T, x, phase=_phase[1], vol0=volx)
            lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, ∂lnϕ∂Ty, voly = ∂lnϕ∂n∂P∂T(modely, p, T, y, phase=_phase[2], vol0=voly)
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

    function jv!(inc)
        return nothing
    end

    return Solvers.NLSolvers.VectorObjective(f!,j!,fj!,jv!) |> Solvers.NLSolvers.NEqProblem
end