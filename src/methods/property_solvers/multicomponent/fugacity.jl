function _select_xy(K,x::Nothing,y)
    return y./ K,y
end

function _select_xy(K,x,y::Nothing)
    return x,K .* x
end

function _select_pT(inc,p,T::Nothing)
    return p,exp(inc[end])
end

function _select_pT(inc,p::Nothing,T)
    return exp(inc[end]),T
end

function _fug_J(J,_x,_y::Nothing,x,y,∂lnϕ∂nx,∂lnϕ∂ny)
    J[diagind(J)] .= 1.
    J[1:(end-1), 1:(end-1)] += (y .* ∂lnϕ∂ny)'
    J[end, 1:(end-1)] = y
    J[end, end] = 0.
    return J
end

function _fug_J(J,_x::Nothing,_y,x,y,∂lnϕ∂nx,∂lnϕ∂ny)
    J[diagind(J)] .= 1.
    J[1:(end-1), 1:(end-1)] += (x .* ∂lnϕ∂nx)'
    J[end, 1:(end-1)] = x
    J[end, end] = 0.
    return J
end
function _tasd()
function f!(F, inc)
    volx, voly = vol_cache[:]
    lnK = inc[1:end-1]
    lnp = inc[end]
    p = exp(lnp)

    y = exp.(lnK) .* x

    lnϕx, volx = lnϕ(model, p, T, x, phase=:liquid, vol0=volx)
    lnϕy, voly = lnϕ(model, p, T, y, phase=:vapor, vol0=voly)

    F[1:end-1] = lnK .+ lnϕy .- lnϕx
    F[end] = sum(y .- x)

    vol_cache[:] .= (volx, voly)
    return F
end

function fj!(F,J,inc)
    volx, voly = vol_cache[:]
    lnK = inc[1:end-1]
    lnp = inc[end]
    K = exp.(lnK)
    p = exp(lnp)

    y = K .* x
    J .= 0
    lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, volx = ∂lnϕ∂n∂P(model, p, T, x, phase=:liquid, vol0=volx)
    lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, voly = ∂lnϕ∂n∂P(model, p, T, y, phase=:vapor, vol0=voly)

    F[1:end-1] = lnK .+ lnϕy .- lnϕx
    F[end] = sum(y .- x)

    J[diagind(J)] .= 1.
    J[1:(end-1), 1:(end-1)] += (y .* ∂lnϕ∂ny)'
    J[end, 1:(end-1)] = K .* x

    J[1:(end-1), end] = p * (∂lnϕ∂Py .- ∂lnϕ∂Px)
    J[end, end] = 0.
    vol_cache[:] .= (volx, voly)

    return F,J


end

function j!(J,inc)
    volx, voly = vol_cache[:]
    lnK = inc[1:end-1]
    lnp = inc[end]
    K = exp.(lnK)
    p = exp(lnp)

    y = K .* x
    J .= 0

    lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, volx = ∂lnϕ∂n∂P(model, p, T, x, phase=:liquid, vol0=volx)
    lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, voly = ∂lnϕ∂n∂P(model, p, T, y, phase=:vapor, vol0=voly)

    J[diagind(J)] .= 1.
    J[1:(end-1), 1:(end-1)] += (y .* ∂lnϕ∂ny)'
    J[end, 1:(end-1)] = K .* x

    J[1:(end-1), end] = p * (∂lnϕ∂Py .- ∂lnϕ∂Px)
    J[end, end] = 0.
    vol_cache[:] .= (volx, voly)
    return J

end

function jv!(inc)
    return nothing
end

return Solvers.NLSolvers.VectorObjective(f!,j!,fj!,jv!) |> Solvers.NLSolvers.NEqProblem

end

function generic_OF_fug(model,_x, _y, _p, _T, vol_cache,_phase)
    function f!(F, inc)
        volx, voly = vol_cache[:]
        lnK = inc[1:end-1]
        K = exp.(lnK)
        p,T = _select_pT(inc,_p,_T)
        x,y = _select_xy(K,_x,_y)

        lnϕx, volx = lnϕ(model, p, T, x, phase=_phase[1], vol0=volx)
        lnϕy, voly = lnϕ(model, p, T, y, phase=_phase[2], vol0=voly)

        F[1:end-1] = lnK .+ lnϕy .- lnϕx
        F[end] = sum(y)  - sum(x)

        vol_cache[:] .= (volx, voly)
        return F
    end

    function fj!(F,J,inc)
        volx, voly = vol_cache[:]
        lnK = inc[1:end-1]
        K = exp.(lnK)
        p,T = _select_pT(inc,_p,_T)
        x,y = _select_xy(K,_x,_y)
        J .= 0.0
        if _T === nothing
            lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, ∂lnϕ∂Tx, volx = ∂lnϕ∂n∂P∂T(model, p, T, x, phase=_phase[1], vol0=volx)
            lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, ∂lnϕ∂Ty, voly = ∂lnϕ∂n∂P∂T(model, p, T, y, phase=_phase[2], vol0=voly)
            J[1:(end-1), end] .= T .* (∂lnϕ∂Ty .- ∂lnϕ∂Tx)
        else
            lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, volx = ∂lnϕ∂n∂P(model, p, T, x, phase=_phase[1], vol0=volx)
            lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, voly = ∂lnϕ∂n∂P(model, p, T, y, phase=_phase[2], vol0=voly)
            J[1:(end-1), end] .= p .* (∂lnϕ∂Py .- ∂lnϕ∂Px)
        end

        F[1:end-1] = lnK .+ lnϕy .- lnϕx
        F[end] = sum(y)  - sum(x)
        _fug_J(J,_x,_y,x,y,∂lnϕ∂nx,∂lnϕ∂ny)
        vol_cache[:] .= (volx, voly)
        return F,J
    end

    function j!(J,inc)
        volx, voly = vol_cache[:]
        lnK = inc[1:end-1]
        K = exp.(lnK)
        p,T = _select_pT(inc,_p,_T)
        x,y = _select_xy(K,_x,_y)
        J .= 0.0
        if _T === nothing
            lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, ∂lnϕ∂Tx, volx = ∂lnϕ∂n∂P∂T(model, p, T, x, phase=_phase[1], vol0=volx)
            lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, ∂lnϕ∂Ty, voly = ∂lnϕ∂n∂P∂T(model, p, T, y, phase=_phase[2], vol0=voly)
            J[1:(end-1), end] .= T .* (∂lnϕ∂Ty .- ∂lnϕ∂Tx)
        else
            lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, volx = ∂lnϕ∂n∂P(model, p, T, x, phase=_phase[1], vol0=volx)
            lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, voly = ∂lnϕ∂n∂P(model, p, T, y, phase=_phase[2], vol0=voly)
            J[1:(end-1), end] .= p .* (∂lnϕ∂Py .- ∂lnϕ∂Px)
        end
        _fug_J(J,_x,_y,x,y,∂lnϕ∂nx,∂lnϕ∂ny)
        vol_cache[:] .= (volx, voly)
        return J
    end

    function jv!(inc)
        return nothing
    end

    return Solvers.NLSolvers.VectorObjective(f!,j!,fj!,jv!) |> Solvers.NLSolvers.NEqProblem
end

function generic_OF_fug(model_long::EoSModel,model_short::EoSModel,_x_long, _x_short, _p, _T, vol_cache,_phase,_view)
    
    function f!(F, inc)
        volx_long, volx_short = vol_cache[:]
        lnK = inc[1:end-1]
        K = exp.(lnK)
        p,T = _select_pT(inc,_p,_T)
        x_long,x_short,y = _select_xy(K,_x_long,_x_short)
        lnϕx_long, volx_long = lnϕ(model_long, p, T, x_long, phase=_phase[1], vol0=volx_long)
        lnϕx_short, volx_short = lnϕ(model_short, p, T, x_short, phase=_phase[2], vol0=volx_short)
        lnϕx_long = lnϕx_long[_view]
        F[1:end-1] = lnK .+ lnϕx_short .- lnϕx_long
        F[end] = sum(y)  - sum(x)
        vol_cache[:] .= (volx_long, volx_short)
        return F
    end

    function fj!(F,J,inc)
        volx_long, volx_short = vol_cache[:]
        lnK = inc[1:end-1]
        K = exp.(lnK)
        p,T = _select_pT(inc,_p,_T)
        x_long,x_short,y = _select_xy(K,_x_long,_x_short)

        if _T === nothing
            lnϕx_long, ∂lnϕ∂nx_long, ∂lnϕ∂Px_long, ∂lnϕ∂Tx_long, volx_long = ∂lnϕ∂n∂P∂T(model, p, T, x_long, phase=_phase[1], vol0=volx_long)
            lnϕx_short, ∂lnϕ∂nx_short, ∂lnϕ∂Px_short, ∂lnϕ∂Tx_short, volx_short = ∂lnϕ∂n∂P∂T(model, p, T, x_short, phase=_phase[2], vol0=volx_short)
            J[1:(end-1), end] .= T .* (∂lnϕ∂Tx_short .- ∂lnϕ∂Tx_long[_view])
        else
            lnϕx_long, ∂lnϕ∂nx_long, ∂lnϕ∂Px_long, volx_long = ∂lnϕ∂n∂P(model, p, T, x_long, phase=_phase[1], vol0=volx_long)
            lnϕx_short, ∂lnϕ∂nx_short, ∂lnϕ∂Px_short, volx_short = ∂lnϕ∂n∂P(model, p, T, x_short, phase=_phase[2], vol0=volx_short)
            J[1:(end-1), end] .= p .* (∂lnϕ∂Px_short .- ∂lnϕ∂Px_long[_view])
        end

        lnϕx_long = lnϕx_long[_view]
        F[1:end-1] = lnK .+ lnϕx_short .- lnϕx_long
        F[end] = sum(y)  - sum(x)
        _fug_J(J,_x_long,_x_short,x_long,x_short,∂lnϕ∂nx_long,∂lnϕ∂nx_short)
        vol_cache[:] .= (volx_long, volx_short)
        return F,J
    end

    function j!(J,inc)
        volx_long, volx_short = vol_cache[:]
        lnK = inc[1:end-1]
        K = exp.(lnK)
        p,T = _select_pT(inc,_p,_T)
        x_long,x_short,y = _select_xy(K,_x_long,_x_short)

        if _p === nothing
            lnϕx_long, ∂lnϕ∂nx_long, ∂lnϕ∂Px_long, ∂lnϕ∂Tx_long, volx_long = ∂lnϕ∂n∂P∂T(model, p, T, x_long, phase=_phase[1], vol0=volx_long)
            lnϕx_short, ∂lnϕ∂nx_short, ∂lnϕ∂Px_short, ∂lnϕ∂Tx_short, volx_short = ∂lnϕ∂n∂P∂T(model, p, T, x_short, phase=_phase[2], vol0=volx_short)
            J[1:(end-1), end] .= T .* (∂lnϕ∂Tx_short .- ∂lnϕ∂Tx_long[_view])
        else
            lnϕx_long, ∂lnϕ∂nx_long, ∂lnϕ∂Px_long, volx_long = ∂lnϕ∂n∂P(model, p, T, x_long, phase=_phase[1], vol0=volx_long)
            lnϕx_short, ∂lnϕ∂nx_short, ∂lnϕ∂Px_short, volx_short = ∂lnϕ∂n∂P(model, p, T, x_short, phase=_phase[2], vol0=volx_short)
            J[1:(end-1), end] .= p .* (∂lnϕ∂Px_short .- ∂lnϕ∂Px_long[_view])
        end
        _fug_J(J,_x_long,_x_short,x_long,x_short,∂lnϕ∂nx_long,∂lnϕ∂nx_short)

        vol_cache[:] .= (volx_long, volx_short)
        return J
    end

    function jv!(inc)
        return nothing
    end

    return Solvers.NLSolvers.VectorObjective(f!,j!,fj!,jv!) |> Solvers.NLSolvers.NEqProblem
end


function __x(model::EoSModel,p,T,x1,x2,_phases;itmax_ss = 5, itmax_newton = 10)
    
    for j in 1:itmax_newton
        for i in 1:itmax_ss

            lnK = lnϕx .- lnϕy
            K = exp.(lnK)

            y_old = 1. * y
            y_calc = x .* K

            # non-volalite
            y_calc[non_volatile] .= 0.

            y = y_calc / sum(y_calc)
            error = sum(abs2, y_old - y)
            # println(i, y, error)
            if error < tol_y
                break
            end

            lnϕx, volx = lnϕ(model, p, T, x, phase=:liquid, vol0=volx)
            lnϕy, voly = lnϕ(model, p, T, y, phase=:vapor, vol0=voly)
        end

        lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, ∂lnϕ∂Tx, volx = ∂lnϕ∂n∂P∂T(model, p, T, x, phase=:liquid, vol0=volx)
        lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, ∂lnϕ∂Ty, voly = ∂lnϕ∂n∂P∂T(model, p, T, y, phase=:vapor, vol0=voly)
        lnK = lnϕx .- lnϕy
        K = exp.(lnK)

        OF = sum(y_calc) - 1.
        dOFdT = sum((x.*K.*(∂lnϕ∂Tx .- ∂lnϕ∂Ty))[volatile])
        dT = OF / dOFdT
        # to avoid negative temperatures
        if dT > T
            dT = 0.2*T
        end

        T -= dT    
    end
end