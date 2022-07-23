"""
OF_bubblepy!(model::EoSModel, x, T, vol_cache)

Objective function to compute bubble pressure using a multidimensional
system of equations via fugacity coefficients.

Inputs:
model: equation of state model
x: liquid phase composition
T: temperature ['K']
vol_cache: array used to update the phases' volumes

Returns: NLSolvers.NEqProblem
"""
function OF_bubblepy!(model, x, T, vol_cache)
# Objetive function to solve bubble point using multidimensional-Newton's method
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

"""
 bubble_pressure_fug(model::EoSModel, T, x, y0, p0; vol0=(nothing,nothing),
                         itmax_newton = 10, itmax_ss = 5, tol_y = 1e-8,
                         tol_p = 1e-8, tol_of = 1e-8)

Function to compute bubble pressure via fugacity coefficients. First it uses
successive substitution to update the phase composition and a outer newtown
loop to update the pressure. If no convergence is reached after itmax_newton
iterations, the system is solved using a multidimensional non-linear
systems of equations.

Inputs:
model: equation of state model
T: bubble temperature ['K']
x: liquid phase composition
y0: initial guess for the vapor phase composition
p0: initial guess for the bubble pressure ['Pa']
vol0: optional, initial guesses for the liquid and vapor phase volumes
itmax_newton: optional, number of iterations to update the pressure using newton's method
itmax_ss: optional, number of iterations to update the liquid phase composition using successive substitution
tol_x: optional, tolerance to stop successive substitution cycle
tol_p: optional, tolerance to stop newton cycle
tol_of: optional, tolerance to check if the objective function is zero.

Returns:
p: bubble pressure
volx: saturared liquid volume
voly: saturared vapor volume
y: saturated vapor composition
"""
function bubble_pressure_fug(model::EoSModel, T, x, y0, p0; vol0=(nothing,nothing),
                         itmax_newton = 10, itmax_ss = 5, tol_y = 1e-8,
                         tol_p = 1e-8, tol_of = 1e-8)

    # Setting the initial guesses for volumes
    vol0 === nothing && (vol0 = (nothing,nothing))
    volx, voly = vol0

    p = 1. * p0
    y = 1. * y0

    nc = length(model)

    # to access this values outside the for loop
    lnϕx = zeros(nc)
    lnϕy = zeros(nc)
    OF = 1.

    for j in 1:itmax_newton

    lnϕx, volx = lnϕ(model, p, T, x, phase=:liquid, vol0=volx)
    lnϕy, voly = lnϕ(model, p, T, y, phase=:vapor, vol0=voly)

    y_calc = 1. * y


    for i in 1:itmax_ss

        lnK = lnϕx .- lnϕy
        K = exp.(lnK)

        y_old = 1. * y
        y_calc = x .* K
        y = y_calc / sum(y_calc)
        error = sum(abs2, y_old - y)
        # println(i, y, error)
        if error < tol_y
            break
        end

        lnϕx, volx = lnϕ(model, p, T, x, phase=:liquid, vol0=volx)
        lnϕy, voly = lnϕ(model, p, T, y, phase=:vapor, vol0=voly)

    end

    lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, volx = ∂lnϕ∂n∂P(model, p, T, x, phase=:liquid, vol0=volx)
    lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, voly = ∂lnϕ∂n∂P(model, p, T, y, phase=:vapor, vol0=voly)
    lnK = lnϕx .- lnϕy
    K = exp.(lnK)

    OF = sum(y_calc) - 1.
    dOFdP = sum(x.*K.*(∂lnϕ∂Px .- ∂lnϕ∂Py))
    dp = OF / dOFdP
    # to avoid negative pressures
    if dp > p
        dp = 0.4*p
    end

    p -= dp

    # println(j, " ", OF, " ", p, " ", dp, " ", y)

    if abs(dp) < tol_p
        break
    end


    end

    if abs(OF) > tol_of
        lnK = lnϕx .- lnϕy
        inc0 = vcat(lnK, log(p))
        vol_cache = [volx, voly]
        problem = OF_bubblepy!(model, x, T, vol_cache)
        sol = Solvers.nlsolve(problem, inc0, Solvers.LineSearch(Solvers.Newton()))
        inc = Solvers.x_sol(sol)
        lnK = inc[1:(end-1)]
        lnp = inc[end]

        y = exp.(lnK) .* x
        p = exp.(lnp)
        volx, voly = vol_cache[:]
        # println("Second order method ", p, " ", y)
    end

    return p, volx, voly, y

end


################ Bubble temperature solver

"""
OF_bubbleTy!(model::EoSModel, y, p, vol_cache)

Objective function to compute bubble temperature using a multidimensional
system of equations via fugacity coefficients.

Inputs:
model: equation of state model
y: vapor phase composition
p: pressure ['Pa']
vol_cache: array used to update the phases' volumes


Returns: NLSolvers.NEqProblem
"""
function OF_bubbleTy!(model, x, p, vol_cache)
# Objetive function to solve bubble point using multidimensional-Newton's method
    function f!(F, inc)
        volx, voly = vol_cache[:]
        lnK = inc[1:end-1]
        lnT = inc[end]
        K = exp.(lnK)
        T = exp(lnT)

        y = K .* x

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
        lnT = inc[end]
        K = exp.(lnK)
        T = exp(lnT)

        y = exp.(lnK) .* x

        lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, ∂lnϕ∂Tx, volx = ∂lnϕ∂n∂P∂T(model, p, T, x, phase=:liquid, vol0=volx)
        lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, ∂lnϕ∂Ty, voly = ∂lnϕ∂n∂P∂T(model, p, T, y, phase=:vapor, vol0=voly)


        F[1:end-1] = lnK .+ lnϕy .- lnϕx
        F[end] = sum(y .- x)

        J[diagind(J)] .= 1.
        J[1:(end-1), 1:(end-1)] += (y .* ∂lnϕ∂ny)'
        J[end, 1:(end-1)] = K.* x

        J[1:(end-1), end] = T * (∂lnϕ∂Ty .- ∂lnϕ∂Tx)
        J[end, end] = 0.
        vol_cache[:] .= (volx, voly)
        return F,J


    end

    function j!(J,inc)
        volx, voly = vol_cache[:]
        lnK = inc[1:end-1]
        lnT = inc[end]
        K = exp.(lnK)
        T = exp(lnT)

        y = exp.(lnK) .* x

        lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, ∂lnϕ∂Tx, volx = ∂lnϕ∂n∂P∂T(model, p, T, x, phase=:liquid, vol0=volx)
        lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, ∂lnϕ∂Ty, voly = ∂lnϕ∂n∂P∂T(model, p, T, y, phase=:vapor, vol0=voly)

        J[diagind(J)] .= 1.
        J[1:(end-1), 1:(end-1)] += (y .* ∂lnϕ∂ny)'
        J[end, 1:(end-1)] = K .* x

        J[1:(end-1), end] = T * (∂lnϕ∂Ty .- ∂lnϕ∂Tx)
        J[end, end] = 0.
        vol_cache[:] .= (volx, voly)
        return J

    end

    function jv!(inc)
        return nothing
    end

    return Solvers.NLSolvers.VectorObjective(f!,j!,fj!,jv!) |> Solvers.NLSolvers.NEqProblem
end

"""
bubble_temperature_fug(model::EoSModel, p, x, y0, T0; vol0=(nothing,nothing),
                         itmax_newton = 10, itmax_ss = 5, tol_y = 1e-8,
                         tol_T = 1e-8, tol_of = 1e-8)

Function to compute bubble temperature via fugacity coefficients. First it uses
successive substitution to update the phase composition and a outer newtown
loop to update the temperature. If no convergence is reached after
itmax_newton iterations, the system is solved using a multidimensional
non-linear systems of equations.

Inputs:
model: equation of state model
P: pressure ['Pa']
x: liquid phase composition
y: initial guess for the vapor phase composition
T0: initial guess for the bubble temperature ['K']
vol0: optional, initial guesses for the liquid and vapor phase volumes
itmax_newton: optional, number of iterations to update the temperature using newton's method
itmax_ss: optional, number of iterations to update the liquid phase composition using successive substitution
tol_x: optional, tolerance to stop successive substitution cycle
tol_T: optional, tolerance to stop newton cycle
tol_of: optional, tolerance to check if the objective function is zero.

Returns:
T: bubble temperature
volx: saturared liquid volume
voly: saturared vapor volume
y: saturated vapor composition
"""
function bubble_temperature_fug(model::EoSModel, p, x, y0, T0; vol0=(nothing,nothing),
                         itmax_newton = 10, itmax_ss = 5, tol_y = 1e-8,
                         tol_T = 1e-8, tol_of = 1e-8)

    # Setting the initial guesses for volumes
    vol0 === nothing && (vol0 = (nothing,nothing))
    volx, voly = vol0

    T = 1. * T0
    y = 1. * y0

    nc = length(model)
    # to access this values outside the for loop
    lnϕx = zeros(nc)
    lnϕy = zeros(nc)
    OF = 1.

    for j in 1:itmax_newton

    lnϕx, volx = lnϕ(model, p, T, x, phase=:liquid, vol0=volx)
    lnϕy, voly = lnϕ(model, p, T, y, phase=:vapor, vol0=voly)

    y_calc = 1. * y

    for i in 1:itmax_ss

        lnK = lnϕx .- lnϕy
        K = exp.(lnK)

        y_old = 1. * y
        y_calc = x .* K
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
    dOFdT = sum(x.*K.*(∂lnϕ∂Tx .- ∂lnϕ∂Ty))
    dT = OF / dOFdT
    # to avoid negative temperatures
    if dT > T
        dT = 0.2*T
    end

    T -= dT

    # println(j, " ", OF, " ", T, " ", dT, " ", y)

    if abs(dT) < tol_T
        break
    end


    end

    # println(T, " ", volx, " ", voly, " ", y)


    if abs(OF) > tol_of
        lnK = lnϕx .- lnϕy
        inc0 = vcat(lnK, log(T))
        vol_cache = [volx, voly]
        problem = OF_bubbleTy!(model, x, p, vol_cache)
        sol = Solvers.nlsolve(problem, inc0, Solvers.LineSearch(Solvers.Newton()))
        inc = Solvers.x_sol(sol)
        lnK = inc[1:(end-1)]
        lnT = inc[end]

        y = exp.(lnK) .* x
        T = exp(lnT)
        volx, voly = vol_cache[:]
        # println("Second order method ", T, " ", y, " ", volx, " ", voly)
    end

    return T, volx, voly, y
end
