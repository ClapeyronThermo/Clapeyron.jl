using LinearAlgebra

########## Dew pressure calculation
"""
OF_dewPx_condensable!(model::EoSModel, y, T, x, vol_cache, condensable)

Objective function to compute dew pressure using a multidimensional
system of equations via fugacity coefficients.

Inputs:
model: equation of state model
y: vapor phase composition
T: temperature ['K']
x: liquid composition array
vol_cache: array used to update the phases' volumes
condensable: array, bools of condensable molecules

Returns: NLSolvers.NEqProblem
"""
function OF_dewPx_condensable!(model, y, T, x, vol_cache, condensable)

    function f!(F, inc)
        volx, voly = vol_cache[:]
        lnK = inc[1:end-1]
        lnp = inc[end]
        K = exp.(lnK)
        p = exp(lnp)

        x[condensable] .= y[condensable] ./ exp.(lnK)

        lnϕx, volx = lnϕ(model, p, T, x, phase=:liquid, vol0=volx)
        lnϕy, voly = lnϕ(model, p, T, y, phase=:vapor, vol0=voly)

        F[1:end-1] = lnK .+ lnϕy[condensable] .- lnϕx[condensable]
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

        x[condensable] .= y[condensable] ./ exp.(lnK)

        lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, ∂lnϕ∂Tx, volx = ∂lnϕ∂n∂P∂T(model, p, T, x, phase=:liquid, vol0=volx)
        lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, ∂lnϕ∂Ty, voly = ∂lnϕ∂n∂P∂T(model, p, T, y, phase=:vapor, vol0=voly)


        F[1:end-1] = lnK .+ lnϕy[condensable] .- lnϕx[condensable]
        F[end] = sum(y .- x)

        J[diagind(J)] .= 1.
        J[1:(end-1), 1:(end-1)] += ((x .* ∂lnϕ∂nx)')[condensable, condensable]
        J[end, 1:(end-1)] = y[condensable] ./ exp.(lnK)
        J[1:(end-1), end] = p * (∂lnϕ∂Py .- ∂lnϕ∂Px)[condensable]
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

        x[condensable] .= y[condensable] ./ exp.(lnK)

        lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, ∂lnϕ∂Tx, volx = ∂lnϕ∂n∂P∂T(model, p, T, x, phase=:liquid, vol0=volx)
        lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, ∂lnϕ∂Ty, voly = ∂lnϕ∂n∂P∂T(model, p, T, y, phase=:vapor, vol0=voly)

        J[diagind(J)] .= 1.
        J[1:(end-1), 1:(end-1)] += ((x .* ∂lnϕ∂nx)')[condensable, condensable]
        J[end, 1:(end-1)] = y[condensable] ./ exp.(lnK)
        J[1:(end-1), end] = p * (∂lnϕ∂Py .- ∂lnϕ∂Px)[condensable]
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
function dew_pressure_fug_condensable(model::EoSModel, T, y, x0, p0; vol0=(nothing,nothing),
                             itmax_newton = 10, itmax_ss = 5, tol_x = 1e-8,
                             tol_p = 1e-8, tol_of = 1e-8,
                             non_condensable_list=[])

Function to compute dew pressure via fugacity coefficients. First it uses
successive substitution to update the phase composition and a outer newtown
loop to update the pressure. If no convergence is reached after itmax_newton
iterations, the system is solved using a multidimensional non-linear
systems of equations.

Inputs:
model: equation of state model
T: dew temperature ['K']
y: vapor phase composition
x0: initial guess for the liquid phase composition
p0: initial guess for the dew pressure ['Pa']
vol0: optional, initial guesses for the liquid and vapor phase volumes
itmax_newton: optional, number of iterations to update the pressure using newton's method
itmax_ss: optional, number of iterations to update the liquid phase composition using successive substitution
tol_x: optional, tolerance to stop successive substitution cycle
tol_p: optional, tolerance to stop newton cycle
tol_of: optional, tolerance to check if the objective function is zero.
non_condensable_list: optional, array with index of molecules that are non-condensable

Returns:
p: dew pressure
volx: saturared liquid volume
voly: saturared vapor volume
x: saturated liquid composition
"""
function dew_pressure_fug_condensable(model::EoSModel, T, y, x0, p0; vol0=(nothing,nothing),
                             itmax_newton = 10, itmax_ss = 5, tol_x = 1e-8,
                             tol_p = 1e-8, tol_of = 1e-8,
                             non_condensable_list=[])

     # Setting the initial guesses for volumes
     vol0 === nothing && (vol0 = (nothing,nothing))
     volx, voly = vol0

     p = 1. * p0
     x = 1. * x0

     nc = length(model)
     # to access this values outside the for loop
     lnϕx = zeros(nc)
     lnϕy = zeros(nc)
     OF = 1.

     # constructing non-condesables list
     non_condensable = Bool.(zeros(nc))
     non_condensable[non_condensable_list] .= true
     condensable = .!non_condensable

     for j in 1:itmax_newton

     lnϕx, volx = lnϕ(model, p, T, x, phase=:liquid, vol0=nothing)
     lnϕy, voly = lnϕ(model, p, T, y, phase=:vapor, vol0=nothing)

     x_calc = 1. * x


     for i in 1:itmax_ss

         lnK = lnϕx .- lnϕy
         K = exp.(lnK)

         x_old = 1. * x
         x_calc = y ./ K
         x_calc[non_condensable] .= 0.
         x = x_calc / sum(x_calc)
         error = sum(abs2, x_old - x)
         # println(i, x, error)
         if error < tol_x
             break
         end

         lnϕx, volx = lnϕ(model, p, T, x, phase=:liquid, vol0=volx)
         lnϕy, voly = lnϕ(model, p, T, y, phase=:vapor, vol0=voly)

     end

     lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, volx = ∂lnϕ∂n∂P(model, p, T, x, phase=:liquid, vol0=volx)
     lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, voly = ∂lnϕ∂n∂P(model, p, T, y, phase=:vapor, vol0=voly)
     lnK = lnϕx .- lnϕy
     K = exp.(lnK)

     OF = sum(x_calc) - 1.
     dOFdp = sum(((- y./K) .* (∂lnϕ∂Px .- ∂lnϕ∂Py))[condensable])
     dp = OF / dOFdp
     if dp > p
         dp = 0.4*p
     end

     p -= dp

     if abs(dp) < tol_p
         break
     end

     end

     if abs(OF) > tol_of
         lnK = lnϕx .- lnϕy
         inc0 = vcat(lnK[condensable], log(p))
         vol_cache = [volx, voly]
         x[non_condensable] .= 0.
         problem = OF_dewPx_condensable!(model, y, T, x, vol_cache, condensable)
         sol = Solvers.nlsolve(problem, inc0, Solvers.LineSearch(Solvers.Newton()))
         inc = Solvers.x_sol(sol)
         lnK = inc[1:(end-1)]
         lnp = inc[end]

         K = exp.(lnK)
         p = exp(lnp)
         x[condensable] .= y[condensable] ./ exp.(lnK)
         volx, voly = vol_cache[:]
         # println("Second order method ", p, " ", x, " ", volx, " ", voly)
     end

     return p, volx, voly, x
end


################# Dew temperature calculation

"""
OF_dewTx_condensable!(model, y, p, x, vol_cache, condensable)

Objective function to compute dew temperature using a multidimensional
system of equations via fugacity coefficients.

Inputs:
model: equation of state model
y: vapor phase composition
P: pressure ['Pa']
x: liquid composition array
vol_cache: array used to update the phases' volumes
condensable: array, bools of condensable molecules


Returns: NLSolvers.NEqProblem
"""
function OF_dewTx_condensable!(model, y, p, x, vol_cache, condensable)

    function f!(F, inc)
        volx, voly = vol_cache[:]
        lnK = inc[1:end-1]
        lnT = inc[end]
        K = exp.(lnK)
        T = exp(lnT)

        x[condensable] .= y[condensable] ./ exp.(lnK)

        lnϕx, volx = lnϕ(model, p, T, x, phase=:liquid, vol0=volx)
        lnϕy, voly = lnϕ(model, p, T, y, phase=:vapor, vol0=voly)

        F[1:end-1] = lnK .+ lnϕy[condensable] .- lnϕx[condensable]
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

        x[condensable] .= y[condensable] ./ exp.(lnK)

        lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, ∂lnϕ∂Tx, volx = ∂lnϕ∂n∂P∂T(model, p, T, x, phase=:liquid, vol0=volx)
        lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, ∂lnϕ∂Ty, voly = ∂lnϕ∂n∂P∂T(model, p, T, y, phase=:vapor, vol0=voly)


        F[1:end-1] = lnK .+ lnϕy[condensable] .- lnϕx[condensable]
        F[end] = sum(y .- x)

        J[diagind(J)] .= 1.
        J[1:(end-1), 1:(end-1)] += ((x .* ∂lnϕ∂nx)')[condensable, condensable]
        J[end, 1:(end-1)] = y[condensable] ./ exp.(lnK)

        J[1:(end-1), end] = T * (∂lnϕ∂Ty .- ∂lnϕ∂Tx)[condensable]
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

        x[condensable] .= y[condensable] ./ exp.(lnK)

        lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, ∂lnϕ∂Tx, volx = ∂lnϕ∂n∂P∂T(model, p, T, x, phase=:liquid, vol0=volx)
        lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, ∂lnϕ∂Ty, voly = ∂lnϕ∂n∂P∂T(model, p, T, y, phase=:vapor, vol0=voly)

        J[diagind(J)] .= 1.
        J[1:(end-1), 1:(end-1)] += ((x .* ∂lnϕ∂nx)')[condensable, condensable]
        J[end, 1:(end-1)] = y[condensable] ./ exp.(lnK)

        J[1:(end-1), end] = T * (∂lnϕ∂Ty .- ∂lnϕ∂Tx)[condensable]
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
dew_temperature_fug_condensable(model::EoSModel, p, y, x0, T0; vol0=(nothing,nothing),
                             itmax_newton = 10, itmax_ss = 5, tol_x = 1e-8,
                             tol_T = 1e-8, tol_of = 1e-8,
                             non_condensable_list=[])

Function to compute dew temperature via fugacity coefficients. First it uses
successive substitution to update the phase composition and a outer newtown
loop to update the temperature. If no convergence is reached after
itmax_newton iterations, the system is solved using a multidimensional
non-linear systems of equations.

Inputs:
model: equation of state model
P: pressure ['Pa']
y: vapor phase composition
x0: initial guess for the liquid phase composition
T0: initial guess for the dew temperature ['K']
vol0: optional, initial guesses for the liquid and vapor phase volumes
itmax_newton: optional, number of iterations to update the temperature using newton's method
itmax_ss: optional, number of iterations to update the liquid phase composition using successive substitution
tol_x: optional, tolerance to stop successive substitution cycle
tol_T: optional, tolerance to stop newton cycle
tol_of: optional, tolerance to check if the objective function is zero.
non_condensable_list: optional, array with index of molecules that are non-condensable

Returns:
T: dew temperature
volx: saturared liquid volume
voly: saturared vapor volume
x: saturated liquid composition
"""
function dew_temperature_fug_condensable(model::EoSModel, p, y, x0, T0; vol0=(nothing,nothing),
                             itmax_newton = 10, itmax_ss = 5, tol_x = 1e-8,
                             tol_T = 1e-8, tol_of = 1e-8,
                             non_condensable_list=[])

     # Setting the initial guesses for volumes
     vol0 === nothing && (vol0 = (nothing,nothing))
     volx, voly = vol0


     T = 1. * T0
     x = 1. * x0

     nc = length(model)
     # to access this values outside the for loop
     lnϕx = zeros(nc)
     lnϕy = zeros(nc)
     OF = 1.

     # constructing non-condesables list
     non_condensable = Bool.(zeros(nc))
     non_condensable[non_condensable_list] .= true
     condensable = .!non_condensable

     for j in 1:itmax_newton

     lnϕx, volx = lnϕ(model, p, T, x, phase=:liquid, vol0=nothing)
     lnϕy, voly = lnϕ(model, p, T, y, phase=:vapor, vol0=nothing)

     x_calc = 1. * x


     for i in 1:itmax_ss

         lnK = lnϕx .- lnϕy
         K = exp.(lnK)

         x_old = 1. * x
         x_calc = y ./ K

         # non - condensable
         x_calc[non_condensable] .= 0.

         x = x_calc / sum(x_calc)
         error = sum(abs2, x_old - x)
         # println(i, x, error)
         if error < tol_x
             break
         end

         lnϕx, volx = lnϕ(model, p, T, x, phase=:liquid, vol0=volx)
         lnϕy, voly = lnϕ(model, p, T, y, phase=:vapor, vol0=voly)

     end

     lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, ∂lnϕ∂Tx, volx = ∂lnϕ∂n∂P∂T(model, p, T, x, phase=:liquid, vol0=volx)
     lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, ∂lnϕ∂Ty, voly = ∂lnϕ∂n∂P∂T(model, p, T, y, phase=:vapor, vol0=voly)
     lnK = lnϕx .- lnϕy
     K = exp.(lnK)

     OF = sum(x_calc) - 1.
     dOFdT = sum(((- y./K) .*(∂lnϕ∂Tx .- ∂lnϕ∂Ty))[condensable])
     dT = OF / dOFdT
     if dT > T
         dT = 0.2*T
     end

     T -= dT

     # println(j, " ", OF, " ", T, " ", dT, " ", x)

     if abs(dT) < tol_T
         break
     end

     end

     # println(T, " ", volx, " ", voly, " ", x)


    if abs(OF) > tol_of
        lnK = lnϕx .- lnϕy
        inc0 = vcat(lnK[condensable], log(T))
        vol_cache = [volx, voly]
        x[non_condensable] .= 0.
        problem = OF_dewTx_condensable!(model, y, p, x, vol_cache, condensable)
        sol = Solvers.nlsolve(problem, inc0, Solvers.LineSearch(Solvers.Newton()))
        inc = Solvers.x_sol(sol)
        lnK = inc[1:(end-1)]
        lnT = inc[end]

        # x = y./ exp.(lnK)
        x[condensable] .= y[condensable] ./ exp.(lnK)
        T = exp(lnT)
        volx, voly = vol_cache[:]
        # println("Second order method ", T, " ", x, " ", volx, " ", voly)
    end
    return T, volx, voly, x
end
