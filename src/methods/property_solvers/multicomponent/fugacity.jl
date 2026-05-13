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
    is_dew(x::BubbleDew) = !is_bubble(x)
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
    _bubble = FugEnum.is_bubble(method)
    _pressure = FugEnum.is_pressure(method)
    converged = false
    tol_stability = abs2(cbrt(tol_xy))
    #caches for ∂lnϕ∂n∂P∂T/∂lnϕ∂n∂P
    lnK,K,w,w_old,w_calc,w_restart,_,Hϕx = cache

    OF = NaN*zero(eltype(lnK))
    if _bubble
        w .= y
        _x,_y = x,w
    else
        w .= x
        _x,_y = w,y
    end

    valid_iter = true
    T_old,p_old = T,p
    for j in 1:itmax_newton
        w_restart .= w_calc
        volx_restart,voly_restart = volx,voly
        valid_iter = true
        if isnan(volx) || isnan(voly)
            break
        end
        error = Inf*one(eltype(lnK))
        for i in 1:itmax_ss
            error < tol_xy && break
            lnϕx, volx = modified_lnϕ(model, p, T, _x, Hϕx, vol0=volx, phase = phasex)
            if isnan(volx)
                lnϕx, volx = lnϕ(model, 1.1p, T, _x, Hϕx, phase = phasex)
            end
            lnK .= lnϕx

            lnϕy, voly = modified_lnϕ(model, p, T, _y, Hϕx, vol0=voly, phase = phasey)
            if isnan(voly)
                lnϕy, voly = lnϕ(model, p, T, _y, Hϕx, phase = phasey)
            end
            lnK .-= lnϕy

            if isnan(volx) || isnan(voly)
                break
            end

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

        if !_bubble && second_order
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

    _bubble,_pressure = FugEnum.is_bubble(method),FugEnum.is_pressure(method)
    tol_stability = abs2(cbrt(tol_xy))
    if _bubble
        n = length(modely)
    else
        n = length(modelx)
    end


    lnK,K,w,w_old,w_calc,w_restart,_,Hϕx,Hϕy,u = cache

    OF = NaN*zero(eltype(lnK))

    if _bubble
        w .= y
        u .= x
        _x,_y = u,w
    else
        w .= x
        u .= y
        _x,_y = w,u
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
        for i in 1:itmax_ss
            error < tol_xy && break
            lnϕx, volx = modified_lnϕ(modelx, p, T, _x, Hϕx, phase=:liquid, vol0=volx)
            lnϕy, voly = modified_lnϕ(modely, p, T, _y, Hϕy, phase=:vapour, vol0=voly)
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

            if _bubble
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
                if _bubble
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
            if _bubble
                _∂lnϕ∂Px = view(∂lnϕ∂Px, _view)
                ∂OF = @sum(w[i]*(_∂lnϕ∂Px[i] - ∂lnϕ∂Py[i]))
            else
                _∂lnϕ∂Py = view(∂lnϕ∂Py,_view)
                ∂OF = @sum(w[i]*(∂lnϕ∂Px[i] - _∂lnϕ∂Py[i]))
            end
        elseif !_pressure && second_order
            ∂lnϕ∂Tx, volx = ∂lnϕ∂T(modelx, p, T, _x, Hϕx, phase=:liquid, vol0=volx)
            ∂lnϕ∂Ty, voly = ∂lnϕ∂T(modely, p, T, _y, Hϕy, phase=:vapour, vol0=voly)
            if _bubble
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

        if !_bubble && second_order
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
    _bubble = FugEnum.is_bubble(method)
    _pressure = FugEnum.is_pressure(method)
    phasex,phasey = FugEnum.phases(method)
    _,K,w,u,_,_,vol_cache,Hϕx = cache
    second_order = !isnothing(J)
    volx,voly = vol_cache[]
    lnK = @view inc[1:end-1]
    K .= exp.(lnK)
    u .= z

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
            !_bubble && _fug_J_∂i∂j!(J,x,∂lnϕ∂nx)
            !isnothing(F) && (F[1:end-1] .= lnK .- lnϕx)
            lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, voly = ∂lnϕ∂n∂P(model, p, T, y, Hϕx, phase=phasey, vol0=voly)
            J1 .+= p .* ∂lnϕ∂Py
            _bubble && _fug_J_∂i∂j!(J,y,∂lnϕ∂ny)
        else
            lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, ∂lnϕ∂Tx, volx = ∂lnϕ∂n∂P∂T(model, p, T, x, Hϕx, phase=phasex, vol0=volx)
            !_bubble && _fug_J_∂i∂j!(J,x,∂lnϕ∂nx)
            J1 .-= T .* ∂lnϕ∂Tx
            !isnothing(F) && (F[1:end-1] .= lnK .- lnϕx)
            lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, ∂lnϕ∂Ty, voly = ∂lnϕ∂n∂P∂T(model, p, T, y, Hϕx, phase=phasey, vol0=voly)
            J1 .+= T .* ∂lnϕ∂Ty
            _bubble && _fug_J_∂i∂j!(J,y,∂lnϕ∂ny)
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
    _bubble = FugEnum.is_bubble(method)
    _pressure = FugEnum.is_pressure(method)
    phasex,phasey = FugEnum.phases(method)
    _,K,w,w2,fw1,fw2,vol_cache,Hϕx,Hϕy,u = cache
    second_order = !isnothing(J)
    volx,voly = vol_cache[]
    lnK = @view inc[1:end-1]
    K .= exp.(lnK)
    u .= z

    if _bubble
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
            if _bubble
                J[1:(end-1), end] .= p .* (∂lnϕ∂Py .- @view(∂lnϕ∂Px[_view]))
            else
                J[1:(end-1), end] .= p .* (@view(∂lnϕ∂Py[_view]) .- ∂lnϕ∂Px)
            end
        else
            lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, ∂lnϕ∂Tx, volx = ∂lnϕ∂n∂P∂T(modelx, p, T, x, Hϕx, phase=phasex, vol0=volx)
            lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, ∂lnϕ∂Ty, voly = ∂lnϕ∂n∂P∂T(modely, p, T, y, Hϕy, phase=phasey, vol0=voly)
            if _bubble
                J[1:(end-1), end] .= T .* (∂lnϕ∂Ty .- @view(∂lnϕ∂Tx[_view]))
            else
                J[1:(end-1), end] .= T .* (@view(∂lnϕ∂Ty[_view]) .- ∂lnϕ∂Tx)
            end
        end
        if _bubble
            _fug_J_∂i∂j!(J,y,∂lnϕ∂ny)
        else
            _fug_J_∂i∂j!(J,x,∂lnϕ∂nx)
        end
        if F !== nothing
            if _bubble
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
        if _bubble
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
        fx = _fug_OF_neq!(nothing,J,x,model,prop,z,data,cache)
        return ∇f
    end

    obj = NLSolvers.VectorObjective(f!,j!,fj!,nothing)
    prob = NEqProblem(obj, inplace = true)
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
        fx = _fug_OF_neq!(nothing,J,x,modelx,modely,prop,z,_view,data,cache)
        return ∇f
    end

    obj = NLSolvers.VectorObjective(f!,j!,fj!,nothing)
    prob = NEqProblem(obj, inplace = true)
end
