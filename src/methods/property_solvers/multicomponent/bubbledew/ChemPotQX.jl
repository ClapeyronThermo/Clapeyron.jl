function _mu_OF_neq!(F,J,inc,model,prop,z,method,cache)
    _bubble = FugEnum.is_bubble(method)
    _pressure = FugEnum.is_pressure(method)
    phasex,phasey = FugEnum.phases(method)
    _,K,w,u,_,_,p_cache,Hϕx = cache
    second_order = !isnothing(J)
    neq = length(inc) - 2 - Int(!_pressure)
    lnK = @view inc[1:neq]
    K .= exp.(lnK)
    u .= z
    vx = exp(inc[neq+1])
    vy = exp(inc[neq+2])
    if _bubble
        w .= K .* u
        x,y = u,w
    else
        w .= u ./ K
        x,y = w,u
    end

    propx = exp(last(inc))
    propy = oftype(propx,prop)
    if _pressure
        p,T = NaN*propy,propy
    else
        p,T = propy,propx
    end
    ps = p_scale(model,z)
    RT = Rgas(model)*T
    if second_order
        J .= 0.0
        Jxv = @view J[1:neq, neq+1]
        Jyv = @view J[1:neq, neq+2]
        if _pressure
            lnfx, ∂lnf∂nx, ∂lnfvx, ∂P∂nx, px, ∂P∂vx = ∂lnf∂n∂V(model, vx, T, x, Hϕx)
            !_bubble && _fug_J_∂i∂j!(J,x,∂lnf∂nx)
            Jxv .-= vx .* ∂lnfvx
            !isnothing(F) && (F[1:neq] .= lnK .- lnfx)
            !_bubble && (J[neq+2, 1:neq] .= 0.0 .- ∂P∂nx .* x ./ ps)

            lnfy, ∂lnf∂ny, ∂lnf∂vy, ∂P∂ny, py, ∂P∂vy = ∂lnf∂n∂V(model, vy, T, y, Hϕx)
            Jyv .+= vy .* ∂lnf∂vy
            _bubble && _fug_J_∂i∂j!(J,y,∂lnf∂ny)
            _bubble && (J[neq+2, 1:neq] .= ∂P∂ny .* y ./ ps)

        else
            lnfx, ∂lnf∂nx, ∂lnfvx, ∂lnf∂Tx, ∂P∂nx, px, ∂P∂vx, ∂P∂Tx = ∂lnf∂n∂V∂T(model, vx, T, x, Hϕx)
            JxT = @view J[1:neq, neq+3]
            Jxv .-= vx .* ∂lnfvx
            JxT .-= T .* ∂lnf∂Tx
            !isnothing(F) && (F[1:neq] .= lnK .- lnfx)
            !_bubble && _fug_J_∂i∂j!(J,x,∂lnf∂nx)
            !_bubble && (J[neq+2, 1:neq] .= 0.0 .- ∂P∂nx .* x ./ ps)
            lnfy, ∂lnf∂ny, ∂lnf∂vy, ∂lnf∂Ty, ∂P∂ny, py, ∂P∂vy, ∂P∂Ty = ∂lnf∂n∂V∂T(model, vy, T, y, Hϕx)
            Jyv = @view J[1:neq, neq+2]
            Jyv .+= vy .* ∂lnf∂vy
            JxT .+= T .* ∂lnf∂Ty
            _bubble && _fug_J_∂i∂j!(J,y,∂lnf∂ny)
            _bubble && (J[neq+2, 1:neq] .= ∂P∂ny .* y ./ ps)
            J[neq+2,neq+3] = T*(∂P∂Ty - ∂P∂Tx)/ps
            J[neq+3,1:neq] .= ∂P∂ny .* y ./ ps
            J[neq+3,neq+2] = ∂P∂vy  * vy / ps
            J[neq+3,neq+3] = T*∂P∂Ty/ps
        end
        J[neq + 2,neq + 1] = -∂P∂vx * vx / ps
        J[neq + 2,neq + 2] = ∂P∂vy * vy / ps
        #_fug_J_∂i∂j!(J,x,y,∂lnϕ∂nx,∂lnϕ∂ny,_bubble)
        if !isnothing(F)
            Fneq = @view F[1:neq]
            Fneq .+= lnfy
            F[neq + 1] = sum(y)  - sum(x)
            F[neq + 2] = (py - px)/ps
            if !_pressure
                F[neq + 3] = (py - p)/ps
            end
        end
    else
        Fneq = @view F[1:neq]
        lnfx, px = lnf(model, vx, T, x, Hϕx)
        Fneq .= lnK .- lnfx
        lnfy, py = lnf(model, vy, T, y, Hϕx)
        Fneq .+= lnfy
        F[neq + 1] = sum(y)  - sum(x)
        F[neq + 2] = (py - px)/ps
        if !_pressure
            F[neq + 3] = (py - p)/ps
        end
    end
    p_cache[] = (px,py)
    return nothing
end

function copy_view!(out,in,_view)
    k = 0
    for i in eachindex(in)
        if _view[i]
            k += 1
            out[k] = in[i]
        end
    end
    #out .= @view in[_view]
    out
end

function _mu_OF_neq!(F,J,inc,modelx,modely,prop,z,_view::V,method,cache) where V
    _bubble = FugEnum.is_bubble(method)
    _pressure = FugEnum.is_pressure(method)
    phasex,phasey = FugEnum.phases(method)
    _,K,w,w2,fw1,fw2,p_cache,Hϕx,Hϕy,u = cache
    second_order = !isnothing(J)
    neq = length(_view)
    lnK = @view inc[1:neq]
    K .= exp.(lnK)
    u .= z

    vx = exp(inc[neq+1])
    vy = exp(inc[neq+2])
    w .= 0
    copy_view!(w,u,_view)
    if _bubble
        w .*= K
        x,y = u,w
    else
        w ./= K
        x,y = w,u
    end

    propx = exp(last(inc))
    propy = oftype(propx,prop)
    if _pressure
        p,T = NaN*propy,propy
    else
        p,T = propy,propx
    end

    ps = p_scale(modelx,z)
    RT = Rgas(modelx)*T

    if second_order
        J .= 0.0
        Jxv = @view J[1:neq, neq+1]
        Jyv = @view J[1:neq, neq+2]
        Jdpdn = @view J[neq+2, 1:neq]
        if _pressure
            # Phase X derivatives
            lnfx, ∂lnf∂nx, ∂lnf∂vx, ∂P∂nx, px, ∂P∂vx = ∂lnf∂n∂V(modelx, vx, T, x, Hϕx)
            lnfy, ∂lnf∂ny, ∂lnf∂vy, ∂P∂ny, py, ∂P∂vy = ∂lnf∂n∂V(modely, vy, T, y, Hϕy)
            
            if _bubble
                _fug_J_∂i∂j!(J, y, ∂lnf∂ny)
                _∂lnf∂vx = copy_view!(fw1,∂lnf∂vx,_view)
                Jxv .-= vx .* _∂lnf∂vx
                Jyv .+= vy .* ∂lnf∂vy
                Jdpdn .= ∂P∂ny .* y ./ ps
            else
                _fug_J_∂i∂j!(J, x, ∂lnf∂nx)
                _∂lnf∂vy = copy_view!(fw1,∂lnf∂vy,_view)
                Jxv .-= vx .* ∂lnf∂vx
                Jyv .+= vy .* _∂lnf∂vy
                Jdpdn .= (0.0 .- ∂P∂nx) .* x ./ ps
            end

        else
            # Temperature case - with temperature derivatives
            lnfx, ∂lnf∂nx, ∂lnf∂vx, ∂lnf∂Tx, ∂P∂nx, px, ∂P∂vx, ∂P∂Tx = ∂lnf∂n∂V∂T(modelx, vx, T, x, Hϕx)
            lnfy, ∂lnf∂ny, ∂lnf∂vy, ∂lnf∂Ty, ∂P∂ny, py, ∂P∂vy, ∂P∂Ty = ∂lnf∂n∂V∂T(modely, vy, T, y, Hϕy)
            JxT = @view J[1:neq, neq+3]

            J[neq+2, neq+3] = T * (∂P∂Ty - ∂P∂Tx) / ps
            J[neq+3, neq+2] .= ∂P∂vy * vy / ps
            J[neq+3, neq+3] .= T * ∂P∂Ty / ps
            if _bubble
                _fug_J_∂i∂j!(J, y, ∂lnf∂ny)
                _∂lnf∂vx = copy_view!(fw1,∂lnf∂vx,_view)
                _∂lnf∂Tx = copy_view!(fw2,∂lnf∂Tx,_view)
                Jxv .-= vx .* _∂lnf∂vx
                Jyv .+= vy .* ∂lnf∂vy
                JxT .+= T .* (∂lnf∂Ty .- _∂lnf∂Tx)
                J[neq+2, 1:neq] .= ∂P∂ny .* y ./ ps
                J[neq+3, 1:neq] .= ∂P∂ny .* y ./ ps
            else
                _fug_J_∂i∂j!(J, x, ∂lnf∂nx)
                _∂lnf∂vy = copy_view!(fw1,∂lnf∂vy,_view)
                _∂lnf∂Ty = copy_view!(fw2,∂lnf∂Ty,_view)
                Jxv .-= vx .* ∂lnf∂vx
                Jyv .+= vy .* _∂lnf∂vy
                JxT .+= T .* (_∂lnf∂Ty .- ∂lnf∂Tx)
                J[neq+2, 1:neq] .= (0.0 .- ∂P∂nx) .* x ./ ps
                _y = copy_view!(fw1,y,_view)
                _∂P∂ny = copy_view!(fw2,∂P∂ny,_view)
                J[neq+3, 1:neq] .= _∂P∂ny .* _y ./ ps
            end
        end
        J[neq + 2, neq + 1] = -∂P∂vx * vx / ps
        J[neq + 2, neq + 2] = ∂P∂vy * vy / ps
        if F !== nothing
            if _bubble
                lnfview = copy_view!(fw1,lnfx,_view)
                F[1:neq] .= lnK .+ lnfy .- lnfview
            else
                lnfview = copy_view!(fw1,lnfy,_view)
                F[1:neq] .= lnK .+ lnfview .- lnfx
            end

            F[neq + 1] = sum(y) - sum(x)
            F[neq + 2] = (py - px) / ps

            if !_pressure
                F[neq + 3] = (py - p) / ps
            end
        end
    else
        # First-order only (no Jacobian)
        lnfx, px = lnf(modelx, vx, T, x, Hϕx)
        lnfy, py = lnf(modely, vy, T, y, Hϕy)

        if _bubble
            lnfview = copy_view!(fw1,lnfx,_view)
            F[1:neq] = lnK .+ lnfy .- lnfview
        else
            lnfview = copy_view!(fw1,lnfy,_view)
            F[1:neq] = lnK .+ lnfview .- lnfx
        end

        F[neq + 1] = sum(y) - sum(x)
        F[neq + 2] = (py - px) / ps

        if !_pressure
            F[neq + 3] = (py - p) / ps
        end
    end

    p_cache[] = (px, py)
    return nothing
end

function _mu_OF_neq(model,prop,z,data,cache)
    function f!(F,x)
        _mu_OF_neq!(F,nothing,x,model,prop,z,data,cache)
        F
    end

    function fj!(F,J,x)
        _mu_OF_neq!(F,J,x,model,prop,z,data,cache)
        F,J
    end

    function j!(J, α)
        fx = _mu_OF_neq!(nothing,J,x,model,prop,z,data,cache)
        return ∇f
    end

    obj = NLSolvers.VectorObjective(f!,j!,fj!,nothing)
    prob = NEqProblem(obj, inplace = true)
end

function _mu_OF_neq(modelx,modely,prop,z,_view,data,cache)
    function f!(F,x)
        _mu_OF_neq!(F,nothing,x,modelx,modely,prop,z,_view,data,cache)
        F
    end

    function fj!(F,J,x)
        _mu_OF_neq!(F,J,x,modelx,modely,prop,z,_view,data,cache)
        F,J
    end

    function j!(J, α)
        fx = _mu_OF_neq!(nothing,J,x,modelx,modely,prop,z,_view,data,cache)
        return ∇f
    end

    obj = NLSolvers.VectorObjective(f!,j!,fj!,nothing)
    prob = NEqProblem(obj, inplace = true)
end

## Bubble pressure solver
struct ChemPotQX{T} <: FlashMethod
    data::FugEnum.BubbleDew
    vol0::Union{Nothing,Tuple{T,T}}
    prop0::Union{Nothing,T}
    w0::Union{Nothing,Vector{T}}
    non_in_w::Union{Nothing,Vector{String}}
    f_limit::Float64
    atol::Float64
    rtol::Float64
    max_iters::Int
    verbose::Bool
end

function ChemPotQX(data;vol0 = nothing,
    prop0 = nothing,
    w0 = nothing,
    non_in_w = nothing,
    f_limit = 0.0,
    atol = 1e-8,
    rtol = 1e-12,
    max_iters = 10^3,
    verbose = false)

    if prop0 == w0 == vol0 == nothing
        return ChemPotQX{Nothing}(data,vol0,prop0,w0,non_in_w,f_limit,atol,rtol,max_iters,verbose)
    elseif (prop0 == w0 == nothing) && !isnothing(vol0)
        vl,vv = promote(vol0[1],vol0[2])
        return ChemPotQX{typeof(vl)}(data,vol0,prop0,w0,non_in_w,f_limit,atol,rtol,max_iters,verbose)
    elseif (vol0 == w0 == nothing) && !isnothing(prop0)
        prop0 = float(prop0)
        return ChemPotQX{typeof(prop0)}(data,vol0,prop0,w0,non_in_w,f_limit,atol,rtol,max_iters,verbose)
    elseif (prop0 == vol0 == nothing) && !isnothing(w0)
        T = eltype(w0)
        return ChemPotQX{T}(data,vol0,prop0,w0,non_in_w,f_limit,atol,rtol,max_iters,verbose)
    elseif !isnothing(vol0) && !isnothing(prop0) && !isnothing(w0)
        vl,vv,prop0,_ = promote(vol0[1],vol0[2],prop0,first(w0))
        T = eltype(vl)
        w0 = convert(Vector{T},w0)
        return ChemPotQX{T}(data,vol0,prop0,w0,non_in_w,f_limit,atol,rtol,max_iters,verbose)
    elseif !isnothing(vol0) && !isnothing(w0)
        vl,vv,_ = promote(vol0[1],vol0[2],first(w0))
        T = eltype(vl)
        w0 = convert(Vector{T},w0)
        return ChemPotQX{T}(data,vol0,prop0,w0,non_in_w,f_limit,atol,rtol,max_iters,verbose)
    elseif  !isnothing(prop0) && !isnothing(w0)
        prop0,_ = promote(prop0,first(w0))
        T = eltype(prop0)
        w0 = convert(Vector{T},w0)
        return ChemPotQX{T}(data,vol0,prop0,w0,non_in_w,f_limit,atol,rtol,max_iters,verbose)
    else
        invalid_bd_error(data)
    end
end

function index_reduction(method::ChemPotQX,idx_r)
    if !isnothing(method.w0)
        method_r = deepcopy(method)
        w0_new = method.w0[idx_r]
        resize!(method_r.w0,length(w0_new))
        method_r.w0 .= w0_new
        return method_r
    end
    return method
end

function Solvers.primalval(method::ChemPotQX{T}) where T
    if T == Nothing
        return Solvers.primalval_struct(method,T)
    else
        return Solvers.primalval_struct(method,Solvers.primal_eltype(T))
    end
end

#like qt, but only 0 or 1
function bdt_flash_impl(model::EoSModel, T, z, method::ChemPotQX)
    data = method.data

    in_equilibria = comps_in_equilibria(component_list(model),method.non_in_w)
    verbose = get_verbosity(method)

    if FugEnum.is_bubble(data)
        p0,vz0,vw0,w00 = bubble_pressure_init(model,T,z,method.vol0,method.prop0,method.w0,in_equilibria,verbose)
    elseif is_lle(data)
        p0,vz0,vw0,w00 = LLE_pressure_init(model,T,z,method.vol0,method.prop0,method.w0,in_equilibria,verbose)
    else
        p0,vw0,vz0,w00 = dew_pressure_init(model,T,z,method.vol0,method.prop0,method.w0,in_equilibria,verbose)
    end

    model_w,_ = index_reduction(model,in_equilibria)
    w0 = w00[in_equilibria]
    neq = count(in_equilibria)
    cache = Clapeyron.fug_bubbledew_cache(model,model_w,T,T,z,z,Val{false}())
    w_r,_,_,_,_,_,p_cache,_,_ = cache
    inc0 = similar(w_r,neq+2)
    inc0_view = @view inc0[1:neq]
    if FugEnum.is_bubble(data) || is_lle(data)
        inc0_view .= log.(w0 ./ @view(z[in_equilibria]))
        inc0[neq+1] = log(vz0)
        inc0[neq+2] = log(vw0)
    else
        inc0_view .= log.(@view(z[in_equilibria]) ./ w0)
        inc0[neq+2] = log(vz0)
        inc0[neq+1] = log(vw0)
    end

    opts = NLSolvers.NEqOptions(method)
    if all(in_equilibria)
        problem = _mu_OF_neq(model,T,z,data,cache)
        sol = Solvers.nlsolve(problem, inc0, Solvers.LineSearch(Solvers.Newton2(inc0),Static(1.0)),opts)
        inc = Solvers.x_sol(sol)
        !all(<(sol.options.f_abstol),sol.info.best_residual) && (inc .= NaN)
    else
        problem = _mu_OF_neq(model,model_w,T,z,in_equilibria,data,cache)
        sol_vol = Solvers.nlsolve(problem, inc0, Solvers.LineSearch(Solvers.Newton2(inc0),Static(1.0)),opts)
        inc = Solvers.x_sol(sol_vol)
        !all(<(sol_vol.options.f_abstol),sol_vol.info.best_residual) && (inc .= NaN)
    end
    lnK = @view inc[1:neq]
    pz,pw = p_cache[]
    p = 0.5*(pw + pw)
    if FugEnum.is_dew(data)
        w_r .= @view(z[in_equilibria]) ./  exp.(lnK)
        vz = exp(inc[neq+2])
        vw = exp(inc[neq+1])
    else
        w_r .= exp.(lnK) .* @view(z[in_equilibria])
        vz = exp(inc[neq+1])
        vw = exp(inc[neq+2])
    end
    w = index_expansion(w_r,in_equilibria)
    if FugEnum.is_bubble(data) || is_lle(data) #q = 0
        return (p, vz, vw, w)
    else #q = 1
        return (p, vw, vz, w)
    end
end

function bdp_flash_impl(model::EoSModel, p, z, method::ChemPotQX)
    data = method.data
    in_equilibria = comps_in_equilibria(component_list(model),method.non_in_w)
    verbose = get_verbosity(method)
    if FugEnum.is_bubble(data)
        T0,vz0,vw0,w00 = bubble_temperature_init(model,p,z,method.vol0,method.prop0,method.w0,in_equilibria,verbose)
    elseif is_lle(data)
        T0,vz0,vw0,w00 = LLE_temperature_init(model,p,z,method.vol0,method.prop0,method.w0,in_equilibria,verbose)
    else
        T0,vw0,vz0,w00 = dew_temperature_init(model,p,z,method.vol0,method.prop0,method.w0,in_equilibria,verbose)
    end

    model_w,_ = index_reduction(model,in_equilibria)
    w0 = w00[in_equilibria]
    neq = count(in_equilibria)
    cache = Clapeyron.fug_bubbledew_cache(model,model_w,p,p,z,z,Val{true}())
    w_r,_,_,_,_,_,_,_,_ = cache
    inc0 = similar(w_r,neq+3)
    inc0_view = @view inc0[1:neq]
    if FugEnum.is_bubble(data) || is_lle(data)
        inc0_view .= log.(w0 ./ @view(z[in_equilibria]))
        inc0[neq+1] = log(vz0)
        inc0[neq+2] = log(vw0)
    else
        inc0_view .= log.(@view(z[in_equilibria]) ./ w0)
        inc0[neq+2] = log(vz0)
        inc0[neq+1] = log(vw0)
    end
    inc0[neq+3] = log(T0)
    opts = NLSolvers.NEqOptions(method)
    if all(in_equilibria)
        problem = _mu_OF_neq(model,p,z,data,cache)
        sol = Solvers.nlsolve(problem, inc0, Solvers.LineSearch(Solvers.Newton2(inc0),Static(1.0)),opts)
        inc = Solvers.x_sol(sol)
        !all(<(sol.options.f_abstol),sol.info.best_residual) && (inc .= NaN)
    else
        problem = _mu_OF_neq(model,model_w,p,z,in_equilibria,data,cache)
        sol_vol = Solvers.nlsolve(problem, inc0, Solvers.LineSearch(Solvers.Newton2(inc0),Static(1.0)),opts)
        inc = Solvers.x_sol(sol_vol)
        !all(<(sol_vol.options.f_abstol),sol_vol.info.best_residual) && (inc .= NaN)
    end
    lnK = @view inc[1:neq]
    T = exp(inc[neq+3])

    if FugEnum.is_dew(data)
        w_r .= @view(z[in_equilibria]) ./  exp.(lnK)
        vz = exp(inc[neq+2])
        vw = exp(inc[neq+1])
    else
        w_r .= exp.(lnK) .* @view(z[in_equilibria])
        vz = exp(inc[neq+1])
        vw = exp(inc[neq+2])
    end
    w = index_expansion(w_r,in_equilibria)
    if FugEnum.is_bubble(data) || is_lle(data) #q = 0
        return (T, vz, vw, w)
    else #q = 1
        return (T, vw, vz, w)
    end
end

function bubble_pressure_impl(model::EoSModel, T, x, method::ChemPotQX)
    data = method.data
    @assert data == FugEnum.BUBBLE_PRESSURE
    return bdt_flash_impl(model,T,x,method)
end

function dew_pressure_impl(model::EoSModel, T, y, method::ChemPotQX)
    data = method.data
    @assert data == FugEnum.DEW_PRESSURE
    return bdt_flash_impl(model,T,y,method)
end

function bubble_temperature_impl(model::EoSModel, p, x, method::ChemPotQX)
    data = method.data
    @assert data == FugEnum.BUBBLE_TEMPERATURE
    return bdp_flash_impl(model,p,x,method)
end

function dew_temperature_impl(model::EoSModel, p, y, method::ChemPotQX)
    data = method.data
    @assert data == FugEnum.DEW_TEMPERATURE
    return bdp_flash_impl(model,p,y,method)
end

function LLE_pressure_impl(model::EoSModel, T, z, method::ChemPotQX)
    data = method.data
    @assert data == FugEnum.LLE_PRESSURE
    return bdt_flash_impl(model,T,z,method)
end

function LLE_temperature_impl(model::EoSModel, p, z, method::ChemPotQX)
    data = method.data
    @assert data == FugEnum.LLE_TEMPERATURE
    return bdp_flash_impl(model,p,z,method)
end

"""
    ChemPotBubblePressure(kwargs...)

Function to compute [`bubble_pressure`](@ref) via chemical potentials.
It directly solves the equality of chemical potentials system of equations.

Inputs:
- `y0 = nothing`: optional, initial guess for the vapor phase composition
- `p0 = nothing`: optional, initial guess for the bubble pressure `[Pa]`
- `vol0 = nothing`: optional, initial guesses for the liquid and vapor phase volumes `[m³]`
- `atol = 1e-8`: optional, absolute tolerance of the non linear system of equations
- `rtol = 1e-12`: optional, relative tolerance of the non linear system of equations
- `max_iters = 1000`: optional, maximum number of iterations
- `nonvolatiles = nothing`: optional, Vector of strings containing non volatile compounds. Those will be set to zero on the vapour phase.
- `verbose = false`: optional, if set to `true`, the method will display additional information in the REPL.
"""
function ChemPotBubblePressure(;vol0 = nothing,
                                p0 = nothing,
                                y0 = nothing,
                                nonvolatiles = nothing,
                                f_limit = 0.0,
                                atol = 1e-8,
                                rtol = 1e-12,
                                max_iters = 10^3,
                                verbose = false)
    prop0 = p0
    w0 = y0
    non_in_w = nonvolatiles
    return ChemPotQX(FugEnum.BUBBLE_PRESSURE;vol0,prop0,w0,non_in_w,f_limit,atol,rtol,max_iters,verbose)
end

"""
    ChemPotBubbleTemperature(kwargs...)

Function to compute [`bubble_temperature`](@ref) via chemical potentials.
It directly solves the equality of chemical potentials system of equations.

Inputs:
- `y0 = nothing`: optional, initial guess for the vapor phase composition.
- `T0 = nothing`: optional, initial guess for the bubble temperature `[K]`.
- `vol0 = nothing`: optional, initial guesses for the liquid and vapor phase volumes `[m³]`
- `atol = 1e-8`: optional, absolute tolerance of the non linear system of equations
- `rtol = 1e-12`: optional, relative tolerance of the non linear system of equations
- `max_iters = 1000`: optional, maximum number of iterations
- `nonvolatiles = nothing`: optional, Vector of strings containing non volatile compounds. Those will be set to zero on the vapour phase.
- `verbose = false`: optional, if set to `true`, the method will display additional information in the REPL.
"""
function ChemPotBubbleTemperature(;vol0 = nothing,
                                    T0 = nothing,
                                    y0 = nothing,
                                    nonvolatiles = nothing,
                                    f_limit = 0.0,
                                    atol = 1e-8,
                                    rtol = 1e-12,
                                    max_iters = 10^3,
                                    verbose = false)

    data = FugEnum.BUBBLE_TEMPERATURE
    prop0 = T0
    w0 = y0
    non_in_w = nonvolatiles
    return ChemPotQX(FugEnum.BUBBLE_TEMPERATURE;vol0,prop0,w0,non_in_w,f_limit,atol,rtol,max_iters,verbose)
end

"""
    ChemPotDewPressure(kwargs...)

Function to compute [`dew_pressure`](@ref) via chemical potentials.
It directly solves the equality of chemical potentials system of equations.

Inputs:
- `x0 = nothing`: optional, initial guess for the liquid phase composition
- `p0 = nothing`: optional, initial guess for the dew pressure `[Pa]`
- `vol0 = nothing`: optional, initial guesses for the liquid and vapor phase volumes `[m³]`
- `atol = 1e-8`: optional, absolute tolerance of the non linear system of equations
- `rtol = 1e-12`: optional, relative tolerance of the non linear system of equations
- `max_iters = 1000`: optional, maximum number of iterations
- `noncondensables = nothing`: optional, Vector of strings containing non condensable compounds. Those will be set to zero on the liquid phase.
- `verbose = false`: optional, if set to `true`, the method will display additional information in the REPL.
"""
function ChemPotDewPressure(;vol0 = nothing,
                                p0 = nothing,
                                x0 = nothing,
                                noncondensables = nothing,
                                f_limit = 0.0,
                                atol = 1e-8,
                                rtol = 1e-12,
                                max_iters = 10^3,
                                verbose = false)

    prop0 = p0
    w0 = x0
    non_in_w = noncondensables
    return ChemPotQX(FugEnum.DEW_PRESSURE;vol0,prop0,w0,non_in_w,f_limit,atol,rtol,max_iters,verbose)
end

"""
    ChemPotDewTemperature(kwargs...)

Function to compute [`dew_temperature`](@ref) via chemical potentials.
It directly solves the equality of chemical potentials system of equations.

Inputs:
- `x0 = nothing`: optional, initial guess for the liquid phase composition
- `T0  = nothing`: optional, initial guess for the dew temperature `[K]`
- `vol0 = nothing`: optional, initial guesses for the liquid and vapor phase volumes `[m³]`
- `atol = 1e-8`: optional, absolute tolerance of the non linear system of equations
- `rtol = 1e-12`: optional, relative tolerance of the non linear system of equations
- `max_iters = 1000`: optional, maximum number of iterations
- `noncondensables = nothing`: optional, Vector of strings containing non condensable compounds. Those will be set to zero on the liquid phase.
- `verbose = false`: optional, if set to `true`, the method will display additional information in the REPL.
"""
function ChemPotDewTemperature(;vol0 = nothing,
    T0 = nothing,
    x0 = nothing,
    noncondensables = nothing,
    f_limit = 0.0,
    atol = 1e-8,
    rtol = 1e-12,
    max_iters = 10^3,
    verbose = false)

    prop0 = T0
    w0 = x0
    non_in_w = noncondensables
    return ChemPotQX(FugEnum.DEW_TEMPERATURE;vol0,prop0,w0,non_in_w,f_limit,atol,rtol,max_iters,verbose)
end

"""
    ChemPotLLEPressure(kwargs...)

Function to compute [`LLE_pressure`](@ref) via chemical potentials.
It directly solves the equality of chemical potentials system of equations.

Inputs:
- `w0 = nothing`: optional, initial guess for the incipent liquid phase composition.
- `p0 = nothing`: optional, initial guess for the LLE pressure `[Pa]`
- `vol0 = nothing`: optional, initial guesses for the bulk liquid and incipient phase volumes `[m³]`
- `atol = 1e-8`: optional, absolute tolerance of the non linear system of equations
- `rtol = 1e-12`: optional, relative tolerance of the non linear system of equations
- `max_iters = 1000`: optional, maximum number of iterations
- `non_in_w = nothing`: optional, Vector of strings containing compounds that will be excluded from the incipient phase.
- `verbose = false`: optional, if set to `true`, the method will display additional information in the REPL.
"""
function ChemPotLLEPressure(;vol0 = nothing,
                                p0 = nothing,
                                w0 = nothing,
                                non_in_w = nothing,
                                f_limit = 0.0,
                                atol = 1e-8,
                                rtol = 1e-12,
                                max_iters = 10^3,
                                verbose = false)

    prop0 = p0
    return ChemPotQX(FugEnum.LLE_PRESSURE;vol0,prop0,w0,non_in_w,f_limit,atol,rtol,max_iters,verbose)
end

"""
    ChemPotLLETemperature(kwargs...)

Function to compute [`LLE_temperature`](@ref) via chemical potentials.
It directly solves the equality of chemical potentials system of equations.

Inputs:
- `w0 = nothing`: optional, initial guess for the incipent liquid phase composition.
- `T0 = nothing`: optional, initial guess for the LLE temperature `[K]`.
- `vol0 = nothing`: optional, initial guesses for the bulk liquid and incipient phase volumes `[m³]`
- `atol = 1e-8`: optional, absolute tolerance of the non linear system of equations
- `rtol = 1e-12`: optional, relative tolerance of the non linear system of equations
- `max_iters = 1000`: optional, maximum number of iterations
- `non_in_w = nothing`: optional, Vector of strings containing compounds that will be excluded from the incipient phase.
- `verbose = false`: optional, if set to `true`, the method will display additional information in the REPL.
"""
function ChemPotLLETemperature(;vol0 = nothing,
                                T0 = nothing,
                                w0 = nothing,
                                non_in_w = nothing,
                                f_limit = 0.0,
                                atol = 1e-8,
                                rtol = 1e-12,
                                max_iters = 10^3,
                                verbose = false)

    prop0 = T0
    return ChemPotQX(FugEnum.LLE_TEMPERATURE;vol0,prop0,w0,non_in_w,f_limit,atol,rtol,max_iters,verbose)
end

#default initializers

function init_preferred_method(method::typeof(bubble_pressure),model::EoSModel,kwargs)
    return ChemPotBubblePressure(;kwargs...)
end

function init_preferred_method(method::typeof(bubble_temperature),model::EoSModel,kwargs)
    return ChemPotBubbleTemperature(;kwargs...)
end

function init_preferred_method(method::typeof(dew_pressure),model::EoSModel,kwargs)
    return ChemPotDewPressure(;kwargs...)
end

function init_preferred_method(method::typeof(dew_temperature),model::EoSModel,kwargs)
    return ChemPotDewTemperature(;kwargs...)
end

function init_preferred_method(method::typeof(LLE_pressure),model::EoSModel,kwargs)
    return ChemPotLLEPressure(;kwargs...)
end

function init_preferred_method(method::typeof(LLE_temperature),model::EoSModel,kwargs)
    return ChemPotLLETemperature(;kwargs...)
end


export ChemPotBubblePressure, ChemPotBubbleTemperature
export ChemPotDewPressure, ChemPotDewTemperature
export ChemPotLLEPressure, ChemPotLLETemperature