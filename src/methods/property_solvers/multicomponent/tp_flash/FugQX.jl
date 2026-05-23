struct FugQX{T} <: FlashMethod
    data::FugEnum.BubbleDew
    vol0::Union{Nothing,Tuple{T,T}}
    prop0::Union{Nothing,T}
    w0::Union{Nothing,Vector{T}}
    non_in_w::Union{Nothing,Vector{String}}
    f_limit::Float64
    atol::Float64
    rtol::Float64
    max_iters::Int
    itmax_newton::Int
    itmax_ss::Int
    tol_w::Float64
    tol_pT::Float64
    tol_of::Float64
    second_order::Bool
    verbose::Bool
end

function Solvers.primalval(method::FugQX{T}) where T
    if T == Nothing
        return Solvers.primalval_struct(method,T)
    else
        return Solvers.primalval_struct(method,Solvers.primal_eltype(T))
    end
end

function FugQX(data;vol0 = nothing,
                    prop0 = nothing,
                    w0 = nothing,
                    nonvolatiles = nothing,
                    f_limit = 0.0,
                    atol = 1e-8,
                    rtol = 1e-12,
                    max_iters = 10^4,
                    itmax_newton = 10,
                    itmax_ss = 5,
                    tol_w = 1e-8,
                    tol_pT = 1e-8,
                    tol_of = 1e-8,
                    second_order = true,
                    verbose = false)

    if prop0 == w0 == vol0 == nothing
        return FugQX{Nothing}(data,vol0,prop0,w0,non_in_w,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_w,tol_pT,tol_of,second_order,verbose)
    elseif (prop0 == w0 == nothing) && !isnothing(vol0)
        vl,vv = promote(vol0[1],vol0[2])
        return FugQX{typeof(vl)}(data,vol0,prop0,w0,non_in_w,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_w,tol_pT,tol_of,second_order,verbose)
    elseif (vol0 == w0 == nothing) && !isnothing(prop0)
        prop0 = float(prop0)
        return FugQX{typeof(prop0)}(data,vol0,prop0,w0,non_in_w,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_w,tol_pT,tol_of,second_order,verbose)
    elseif (prop0 == vol0 == nothing) && !isnothing(w0)
        T = eltype(w0)
        return FugQX{T}(data,vol0,prop0,w0,non_in_w,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_w,tol_pT,tol_of,second_order,verbose)
    elseif !isnothing(vol0) && !isnothing(prop0) && !isnothing(w0)
        vl,vv,prop0,_ = promote(vol0[1],vol0[2],prop0,first(w0))
        T = eltype(vl)
        w0 = convert(Vector{T},w0)
        return FugQX{T}(data,vol0,prop0,w0,non_in_w,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_w,tol_pT,tol_of,second_order,verbose)
    elseif !isnothing(vol0) && !isnothing(w0)
        vl,vv,_ = promote(vol0[1],vol0[2],first(w0))
        T = eltype(vl)
        w0 = convert(Vector{T},w0)
        return FugQX{T}(data,vol0,prop0,w0,non_in_w,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_w,tol_pT,tol_of,second_order,verbose)
    elseif  !isnothing(prop0) && !isnothing(w0)
        prop0,_ = promote(prop0,first(w0))
        T = eltype(prop0)
        w0 = convert(Vector{T},w0)
        return FugQX{T}(data,vol0,prop0,w0,non_in_w,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_w,tol_pT,tol_of,second_order,verbose)
    else
        invalid_bd_error(data)
    end
end

function bdt_flash_impl(model::EoSModel, T, z, method::FugQX)
    data = method.data

    in_equilibria = comps_in_equilibria(component_list(model),method.non_in_w)
    verbose = get_verbosity(method)

    if is_bubble(data)
        p0,vz0,vw0,w00 = bubble_pressure_init(model,T,z,method.vol0,method.prop0,method.w0,in_equilibria,verbose)
    elseif is_lle(data)
        p0,vz0,vw0,w00 = LLE_pressure_init(model,T,z,method.vol0,method.prop0,method.w0,in_equilibria,verbose)
    else
        p0,vw0,vz0,w00 = dew_pressure_init(model,T,z,method.vol0,method.prop0,method.w0,in_equilibria,verbose)

    end
    zz = similar(w00)
    zz .= z
    model_w,_ = index_reduction(model,in_equilibria)
    w0 = w00[in_equilibria]
    neq = count(in_equilibria)
    cache = Clapeyron.fug_bubbledew_cache(model,model_w,T,T,z,z,Val{false}())

    method_data = FugData(data,
                    method.itmax_ss,
                    method.itmax_newton,
                    method.tol_pT,
                    method.tol_w,
                    method.tol_of,
                    method.second_order,
                    false)

    cache = fug_bubbledew_cache(model,model_y,p0,T,zz,w0,Val{false}())

    if all(in_equilibria)
        if is_dew(data)
            converged,res_ss = _fug_OF_ss(model,p0,T,w00,zz,vol0,method_data,cache)
        else
            converged,res_ss = _fug_OF_ss(model,p0,T,zz,w00,vol0,method_data,cache)
        end

    else
        if is_dew(data)
            converged,res_ss = _fug_OF_ss(model_w,model,p0,T,w00,zz,vol0,in_equilibria,method_data,cache)
        else
            converged,res_ss = _fug_OF_ss(model,model_w,p0,T,zz,w00,vol0,in_equilibria,method_data,cache)
        end

    end

    p_ss,T_ss,x_ss,y_ss,vol_ss = res_ss
    volx_ss,voly_ss = vol_ss
    phasex,phasey = FugEnum.phases(data)

    if converged || isnan(volx_ss) || isnan(voly_ss)
        if iszero(volx_ss) && model isa PTFlashWrapper
            vx = volume(model,p_ss,T,x_ss,phase = phasex)
            volx_ss = oftype(volx_ss,x_ss)
        end

        if iszero(voly_ss) && model isa PTFlashWrapper
            vy = volume(model,p_ss,T,y_ss,phase = phasey)
            voly_ss = oftype(voly_ss,vy)
        end

        if is_bubble(data) || is_lle(data)
            w_ss = index_expansion(y_ss,in_equilibria)
        else
            w_ss = index_expansion(x_ss,in_equilibria)
        end
        return p_ss,volx_ss,voly_ss,w_ss
    end

    lnK,K,w,w_old,w_calc,w_restart,vol_cache,Hϕx,Hϕy = cache
    inc0 = vcat(lnK, log(p))
    vol_cache[] = (volx_ss,voly_ss)
    opts = NLSolvers.NEqOptions(method)
    if all(in_equilibria)
        problem = _fug_OF_neq(model,T,zz,data,cache)
        sol = Solvers.nlsolve(problem, inc0, Solvers.LineSearch(Solvers.Newton2(inc0)),opts)
        inc = Solvers.x_sol(sol)
        !__check_convergence(sol) && (inc.= NaN)
    else
        problem = _fug_OF_neq(model,model_y,T,zz,in_equilibria,data,cache)
        sol = Solvers.nlsolve(problem, inc0, Solvers.LineSearch(Solvers.Newton2(inc0)),opts)
        inc = Solvers.x_sol(sol)
        !__check_convergence(sol) && (inc.= NaN)
    end

    lnK = @view inc[1:(end-1)]
    volx,voly = vol_cache[]
    p = exp(inc[end])
    if is_dew(data)
        w_r .= @view(z[in_equilibria]) ./  exp.(lnK)
    else
        w_r .= exp.(lnK) .* @view(z[in_equilibria])
    end
    w = index_expansion(w_r,in_equilibria)
    return (p, volx, voly, w)
end

function bdp_flash_impl(model::EoSModel, p, z, method::FugQX)
    data = method.data

    in_equilibria = comps_in_equilibria(component_list(model),method.non_in_w)
    verbose = get_verbosity(method)

    if is_bubble(data)
        T0,vz0,vw0,w00 = bubble_temperature_init(model,p,z,method.vol0,method.prop0,method.w0,in_equilibria,verbose)
    elseif is_lle(data)
        T0,vz0,vw0,w00 = LLE_temperature_init(model,p,z,method.vol0,method.prop0,method.w0,in_equilibria,verbose)
    else
        T0,vw0,vz0,w00 = dew_temperature_init(model,p,z,method.vol0,method.prop0,method.w0,in_equilibria,verbose)

    end
    zz = similar(w00)
    zz .= z
    model_w,_ = index_reduction(model,in_equilibria)
    w0 = w00[in_equilibria]
    neq = count(in_equilibria)
    cache = Clapeyron.fug_bubbledew_cache(model,model_w,p,p,z,z,Val{false}())

    method_data = FugData(data,
                    method.itmax_ss,
                    method.itmax_newton,
                    method.tol_pT,
                    method.tol_w,
                    method.tol_of,
                    method.second_order,
                    false)

    cache = fug_bubbledew_cache(model,model_y,p0,T,zz,w0,Val{false}())

    if all(in_equilibria)
        if is_dew(data)
            converged,res_ss = _fug_OF_ss(model,p,T0,w00,zz,vol0,method_data,cache)
        else
            converged,res_ss = _fug_OF_ss(model,p,T0,zz,w00,vol0,method_data,cache)
        end

    else
        if is_dew(data)
            converged,res_ss = _fug_OF_ss(model_w,model,p,T0,w00,zz,vol0,in_equilibria,method_data,cache)
        else
            converged,res_ss = _fug_OF_ss(model,model_w,p,T0,zz,w00,vol0,in_equilibria,method_data,cache)
        end

    end

    p_ss,T_ss,x_ss,y_ss,vol_ss = res_ss
    volx_ss,voly_ss = vol_ss
    phasex,phasey = FugEnum.phases(data)

    if converged || isnan(volx_ss) || isnan(voly_ss)
        if iszero(volx_ss) && model isa PTFlashWrapper
            vx = volume(model,p,T_ss,x_ss,phase = phasex)
            volx_ss = oftype(volx_ss,x_ss)
        end

        if iszero(voly_ss) && model isa PTFlashWrapper
            vy = volume(model,p,T_ss,y_ss,phase = phasey)
            voly_ss = oftype(voly_ss,vy)
        end

        if is_bubble(data) || is_lle(data)
            w_ss = index_expansion(y_ss,in_equilibria)
        else
            w_ss = index_expansion(x_ss,in_equilibria)
        end
        return p_ss,volx_ss,voly_ss,w_ss
    end

    lnK,K,w,w_old,w_calc,w_restart,vol_cache,Hϕx,Hϕy = cache
    inc0 = vcat(lnK, log(T))
    vol_cache[] = (volx_ss,voly_ss)
    opts = NLSolvers.NEqOptions(method)
    if all(in_equilibria)
        problem = _fug_OF_neq(model,T,zz,data,cache)
        sol = Solvers.nlsolve(problem, inc0, Solvers.LineSearch(Solvers.Newton2(inc0)),opts)
        inc = Solvers.x_sol(sol)
        !__check_convergence(sol) && (inc.= NaN)
    else
        problem = _fug_OF_neq(model,model_y,T,zz,in_equilibria,data,cache)
        sol = Solvers.nlsolve(problem, inc0, Solvers.LineSearch(Solvers.Newton2(inc0)),opts)
        inc = Solvers.x_sol(sol)
        !__check_convergence(sol) && (inc.= NaN)
    end

    lnK = @view inc[1:(end-1)]
    volx,voly = vol_cache[]
    T = exp(inc[end])
    volx,voly = vol_cache[]
    if is_dew(data)
        w_r .= @view(z[in_equilibria]) ./  exp.(lnK)
    else
        w_r .= exp.(lnK) .* @view(z[in_equilibria])
    end
    w = index_expansion(w_r,in_equilibria)
    return (p, volx, voly, w)
end

function bubble_pressure_impl(model::EoSModel, T, x, method::FugQX)
    data = method.data
    @assert data == FugEnum.BUBBLE_PRESSURE
    if has_a_res(model)
        return bdt_flash_impl(model,T,x,method)
    else
        cached_model = __tpflash_cache_model(model,NaN,T,x,:vle)
        return bdt_flash_impl(model,T,x,method)
    end
end

function dew_pressure_impl(model::EoSModel, T, y, method::FugQX)
    data = method.data
    @assert data == FugEnum.DEW_PRESSURE
    if has_a_res(model)
        return bdt_flash_impl(model,T,y,method)
    else
        cached_model = __tpflash_cache_model(model,NaN,T,y,:vle)
        return bdt_flash_impl(model,T,y,method)
    end
end

function bubble_temperature_impl(model::EoSModel, p, x, method::FugQX)
    data = method.data
    @assert data == FugEnum.BUBBLE_TEMPERATURE
    if has_a_res(model)
        return bdp_flash_impl(model,p,x,method)
    else
        cached_model = __tpflash_cache_model(model,p,NaN,x,:vle)
        return bdp_flash_impl(cached_model,p,x,method)
    end
end

function dew_temperature_impl(model::EoSModel, p, y, method::FugQX)
    data = method.data
    @assert data == FugEnum.DEW_TEMPERATURE
    if has_a_res(model)
        return bdp_flash_impl(model,p,y,method)
    else
        cached_model = __tpflash_cache_model(model,p,NaN,y,:vle)
        return bdp_flash_impl(cached_model,p,y,method)
    end
end

function LLE_pressure_impl(model::EoSModel, T, z, method::FugQX)
    data = method.data
    @assert data == FugEnum.LLE_PRESSURE
    if has_a_res(model)
        return bdt_flash_impl(model,T,z,method)
    else
        cached_model = __tpflash_cache_model(model,NaN,T,z,:lle)
        return bdt_flash_impl(cached_model,T,z,method)
    end
end

function LLE_temperature_impl(model::EoSModel, p, z, method::FugQX)
    data = method.data
    @assert data == FugEnum.LLE_TEMPERATURE
    if has_a_res(model)
        return bdp_flash_impl(model,p,z,method)
    else
        cached_model = __tpflash_cache_model(model,p,NaN,z,:vle)
        return bdp_flash_impl(cached_model,p,z,method)
    end
end

"""
    FugBubblePressure(kwargs...)

Function to compute [`bubble_pressure`](@ref) via fugacity coefficients. First it uses
successive substitution to update the phase composition and an outer Newton's
loop to update the pressure. If no convergence is reached after `itmax_newton`
iterations, the system is solved using a multidimensional non-linear
system of equations.

Inputs:
- `y0 = nothing`: optional, initial guess for the vapor phase composition
- `p0 = nothing`: optional, initial guess for the bubble pressure `[Pa]`
- `vol0 = nothing`: optional, initial guesses for the liquid and vapor phase volumes `[m³]`
- `itmax_newton = 10`: optional, number of iterations to update the pressure using Newton's method
- `itmax_ss = 5`: optional, number of iterations to update the liquid phase composition using successive substitution
- `tol_x = 1e-8`: optional, tolerance to stop successive substitution cycle
- `tol_p = 1e-8`: optional, tolerance to stop Newton's cycle
- `tol_of = 1e-8`: optional, tolerance to check if the objective function is zero.
- `nonvolatiles = nothing`: optional, Vector of strings containing non volatile compounds. Those will be set to zero on the vapour phase.
- `second_order`: optional, decide if the algorithm uses second order information when updating the guess estimates. Second order methods are slower, but more reliable.
- `verbose = false`: optional, if set to `true`, the method will display additional information in the REPL.
"""
function FugBubblePressure(;vol0 = nothing,
                                p0 = nothing,
                                y0 = nothing,
                                nonvolatiles = nothing,
                                f_limit = 0.0,
                                atol = 1e-8,
                                rtol = 1e-12,
                                max_iters = 10^4,
                                itmax_newton = 10,
                                itmax_ss = 5,
                                tol_y = 1e-8,
                                tol_p = 1e-8,
                                tol_of = 1e-8,
                                second_order = true,
                                verbose = false)
    w0 = y0
    prop0 = p0
    tol_pT = tol_p
    tol_w = tol_y
    non_in_w = nonvolatiles
    return FugQX(FugEnum.BUBBLE_PRESSURE;vol0,prop0,w0,non_in_w,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_w,tol_pT,tol_of,second_order,verbose)
end

"""
    FugBubbleTemperature(kwargs...)

Function to compute [`bubble_temperature`](@ref)  via fugacity coefficients.
First it uses successive substitution to update the phase composition and an outer Newton's (or secant) loop to update the temperature.
If no convergence is reached after `itmax_newton` iterations, the system is solved using a multidimensional non-linear system of equations.


Inputs:
- `y = nothing`: optional, initial guess for the vapor phase composition.
- `T0 = nothing`: optional, initial guess for the bubble temperature `[K]`.
- `vol0 = nothing`: optional, initial guesses for the liquid and vapor phase volumes `[m³]`
- `itmax_newton = 10`: optional, number of iterations to update the pressure using Newton's (or secant) method
- `itmax_ss = 5`: optional, number of iterations to update the liquid phase composition using successive substitution
- `tol_x = 1e-8`: optional, tolerance to stop successive substitution cycle
- `tol_T = 1e-8`: optional, tolerance to stop Newton's cycle
- `tol_of = 1e-8`: optional, tolerance to check if the objective function is zero.
- `nonvolatiles`: optional, Vector of strings containing non volatile compounds. Those will be set to zero on the vapour phase.
- `second_order`: optional, decide if the algorithm uses second order information when updating the guess estimates. Second order methods are slower, but more reliable.
- `verbose = false`: optional, if set to `true`, the method will display additional information in the REPL.
"""
function FugBubbleTemperature(;vol0 = nothing,
                                T0 = nothing,
                                y0 = nothing,
                                nonvolatiles = nothing,
                                f_limit = 0.0,
                                atol = 1e-8,
                                rtol = 1e-12,
                                max_iters = 10^4,
                                itmax_newton = 10,
                                itmax_ss = 5,
                                tol_y = 1e-8,
                                tol_T = 1e-8,
                                tol_of = 1e-8,
                                second_order = true,
                                verbose = false)

    w0 = y0
    prop0 = T0
    tol_pT = tol_T
    tol_w = tol_y
    non_in_w = nonvolatiles
    return FugQX(FugEnum.BUBBLE_TEMPERATURE;vol0,prop0,w0,non_in_w,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_w,tol_pT,tol_of,second_order,verbose)
end

"""
    FugDewPressure(kwargs...)

Method to compute [`dew_pressure`](@ref) via fugacity coefficients. First it uses
successive substitution to update the phase composition and an outer Newton's
loop to update the pressure. If no convergence is reached after `itmax_newton`
iterations, the system is solved using a multidimensional non-linear
system of equations.

Inputs:
- `x0 = nothing`: optional, initial guess for the liquid phase composition
- `p0 = nothing`: optional, initial guess for the dew pressure `[Pa]`
- `vol0 = nothing`: optional, initial guesses for the liquid and vapor phase volumes `[m³]`
- `itmax_newton = 10`: optional, number of iterations to update the pressure using Newton's method
- `itmax_ss = 5`: optional, number of iterations to update the liquid phase composition using successive substitution
- `tol_x = 1e-8`: optional, tolerance to stop successive substitution cycle
- `tol_p = 1e-8`: optional, tolerance to stop Newton's cycle
- `tol_of = 1e-8`: optional, tolerance to check if the objective function is zero.
- `noncondensables = nothing`: optional, Vector of strings containing non condensable compounds. Those will be set to zero on the liquid phase.
- `second_order`: optional, decide if the algorithm uses second order information when updating the guess estimates. Second order methods are slower, but more reliable.
- `verbose = false`: optional, if set to `true`, the method will display additional information in the REPL.
"""
function FugDewPressure(;vol0 = nothing,
                                p0 = nothing,
                                x0 = nothing,
                                noncondensables = nothing,
                                f_limit = 0.0,
                                atol = 1e-8,
                                rtol = 1e-12,
                                max_iters = 10^4,
                                itmax_newton = 10,
                                itmax_ss = 5,
                                tol_x = 1e-8,
                                tol_p = 1e-8,
                                tol_of = 1e-8,
                                second_order = true,
                                verbose = false)

    w0 = x0
    prop0 = p0
    tol_pT = tol_p
    tol_w = tol_x
    non_in_w = noncondensables
    return FugQX(FugEnum.DEW_PRESSURE;vol0,prop0,w0,non_in_w,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_w,tol_pT,tol_of,second_order,verbose)
end

"""
    FugDewTemperature(kwargs...)

Method to compute [`dew_temperature`](@ref) via fugacity coefficients. First it uses
successive substitution to update the phase composition and an outer Newton's
loop to update the temperature. If no convergence is reached after
`itmax_newton` iterations, the system is solved using a multidimensional
non-linear system of equations.

Inputs:
- `x0 = nothing`: optional, initial guess for the liquid phase composition
- `T0 = nothing`: optional, initial guess for the dew temperature `[K]`
- `vol0 = nothing`: optional, initial guesses for the liquid and vapor phase volumes `[m³]`
- `itmax_newton = 10`: optional, number of iterations to update the temperature using Newton's method
- `itmax_ss = 5`: optional, number of iterations to update the liquid phase composition using successive substitution
- `tol_x = 1e-8`: optional, tolerance to stop successive substitution cycle
- `tol_T = 1e-8`: optional, tolerance to stop Newton's cycle
- `tol_of = 1e-8`: optional, tolerance to check if the objective function is zero.
- `noncondensables = nothing`: optional, Vector of strings containing non condensable compounds. Those will be set to zero on the liquid phase.
- `second_order`: optional, decide if the algorithm uses second order information when updating the guess estimates. Second order methods are slower, but more reliable.
- `verbose = false`: optional, if set to `true`, the method will display additional information in the REPL.
"""
function FugDewTemperature(;vol0 = nothing,
                            T0 = nothing,
                            x0 = nothing,
                            noncondensables = nothing,
                            f_limit = 0.0,
                            atol = 1e-8,
                            rtol = 1e-12,
                            max_iters = 10^4,
                            itmax_newton = 10,
                            itmax_ss = 5,
                            tol_x = 1e-8,
                            tol_T = 1e-8,
                            tol_of = 1e-8,
                            second_order = true,
                            verbose = false)

    w0 = x0
    prop0 = T0
    tol_pT = tol_T
    tol_w = tol_x
    non_in_w = noncondensables
    return FugQX(FugEnum.DEW_TEMPERATURE;vol0,prop0,w0,non_in_w,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_w,tol_pT,tol_of,second_order,verbose)
end

"""
    FugLLEPressure(kwargs...)

Function to compute [`LLE_pressure`](@ref) via fugacity coefficients. First it uses
successive substitution to update the phase composition and an outer Newton's
loop to update the pressure. If no convergence is reached after `itmax_newton`
iterations, the system is solved using a multidimensional non-linear
system of equations.

Inputs:
- `w0 = nothing`: optional, initial guess for the incipent liquid phase composition.
- `p0 = nothing`: optional, initial guess for the LLE pressure `[Pa]`
- `vol0 = nothing`: optional, initial guesses for the bulk liquid and incipient phase volumes `[m³]`
- `itmax_newton = 10`: optional, number of iterations to update the pressure using Newton's method
- `itmax_ss = 5`: optional, number of iterations to update the bulk liquid phase composition using successive substitution
- `tol_x = 1e-8`: optional, tolerance to stop successive substitution cycle
- `tol_p = 1e-8`: optional, tolerance to stop Newton's cycle
- `tol_of = 1e-8`: optional, tolerance to check if the objective function is zero.
- `non_in_w = nothing`: optional, Vector of strings containing compounds that will be excluded from the incipient phase.
- `second_order`: optional, decide if the algorithm uses second order information when updating the guess estimates. Second order methods are slower, but more reliable.
- `verbose = false`: optional, if set to `true`, the method will display additional information in the REPL.
"""
function FugLLEPressure(;vol0 = nothing,
                                p0 = nothing,
                                w0 = nothing,
                                non_in_w = nothing,
                                f_limit = 0.0,
                                atol = 1e-8,
                                rtol = 1e-12,
                                max_iters = 10^4,
                                itmax_newton = 10,
                                itmax_ss = 5,
                                tol_w = 1e-8,
                                tol_p = 1e-8,
                                tol_of = 1e-8,
                                second_order = true,
                                verbose = false)
    prop0 = p0
    tol_pT = tol_p
    return FugQX(FugEnum.LLE_PRESSURE;vol0,prop0,w0,non_in_w,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_w,tol_pT,tol_of,second_order,verbose)
end

"""
    FugLLETemperature(kwargs...)

Function to compute [`LLE_temperature`](@ref) via fugacity coefficients.
First it uses successive substitution to update the phase composition and an outer Newton's (or secant) loop to update the temperature.
If no convergence is reached after `itmax_newton` iterations, the system is solved using a multidimensional non-linear system of equations.


Inputs:
- `w0 = nothing`: optional, initial guess for the incipent liquid phase composition.
- `T0 = nothing`: optional, initial guess for the LLE temperature `[K]`.
- `vol0 = nothing`: optional, initial guesses for the bulk liquid and incipient phase volumes `[m³]`
- `itmax_newton = 10`: optional, number of iterations to update the pressure using Newton's (or secant) method
- `itmax_ss = 5`: optional, number of iterations to update the bulk liquid phase composition using successive substitution
- `tol_x = 1e-8`: optional, tolerance to stop successive substitution cycle
- `tol_T = 1e-8`: optional, tolerance to stop Newton's cycle
- `tol_of = 1e-8`: optional, tolerance to check if the objective function is zero.
- `non_in_w = nothing`: optional, Vector of strings containing compounds that will be excluded from the incipient phase.
- `second_order`: optional, decide if the algorithm uses second order information when updating the guess estimates. Second order methods are slower, but more reliable.
- `verbose = false`: optional, if set to `true`, the method will display additional information in the REPL.
"""
function FugLLETemperature(;vol0 = nothing,
                                T0 = nothing,
                                w0 = nothing,
                                non_in_w = nothing,
                                f_limit = 0.0,
                                atol = 1e-8,
                                rtol = 1e-12,
                                max_iters = 10^4,
                                itmax_newton = 10,
                                itmax_ss = 5,
                                tol_w = 1e-8,
                                tol_T = 1e-8,
                                tol_of = 1e-8,
                                second_order = true,
                                verbose = false)

    prop0 = T0
    tol_pT = tol_T
    return FugQX(FugEnum.LLE_TEMPERATURE;vol0,prop0,w0,non_in_w,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_w,tol_pT,tol_of,second_order,verbose)
end

export FugBubblePressure, FugBubbleTemperature
export FugDewPressure, FugDewTemperature
export FugLLEPressure, FugLLETemperature