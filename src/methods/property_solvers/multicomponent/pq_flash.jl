function pq_flash_x0(model,p,Q,z,T00)

    #Algorithm 1: An initialization subroutine for the PQ-flash inside-out algorithm

    #modified to use some Clapeyron initialization functions.
    pure = split_model(model)

    #1. T guess. original: 0.8*dot(z,Tc). here: 0.5*dot(z,Ti) where Ti = extended_saturation_temperature(model,p)
    #   guess α = 0.5
    
    
    #2. Solve 5) and 11) for T,α  assuming ideal thermophysical properties and Raoult’s Law
    T,α = pq_flash_x0_step_2(model,p,Q,z,T0,α0)
    
    #3,4: Tref = T - 1, T̃ = T + 1
    Tᵣ = T - 1
    T̃ = T + 1
    Tₘ = T #midpoint

    #5: calculate x,y (we use thermopack aproximation instead of raoult, that uses real propery models)
    Kₘ = suggest_K(model,p,T,z,pure)
    y = rr_flash_vapor(Kₘ,z,α)
    x = rr_flash_liquid(Kₘ,z,α)

    #6,7 calculate K(Tr), K(T̃), using real property models
    Kᵣ = suggest_K(model,p,Tᵣ,z,pure)
    K̃ = suggest_K(model,p,T̃,z,pure)
    
    logK̃ = log.(K̃)
    logKᵣ = log.(Kᵣ)
    logKₘ = log.(Kₘ)
    #finite differences
    dKᵣdT = logKₘ .- logKᵣ
    dK̃dT = logK̃ .- logKₘ
    ti = similar(Kᵣ)
    #8. calculate kb and k̃b using 15) - 17)
    ti .= y .* dKᵣdT ./(1 .+ α.*(Kᵣ .-1))
    ti ./= sum(ti)
    kb = exp(dot(ti,logKᵣ))

    ti .= y .* dK̃dT ./(1 .+ α.*(K̃ .-1))
    ti ./= sum(ti)
    k̃b = exp(dot(ti,logK̃))

    #9. set kb0 <- k̃b
    kb₀ = k̃b

    #10. calculate u,A,B from 14) 19), 20)
    u = log.(Kₘ ./ kb)
    B = log(k̃b/kb)/(1/T̃ - 1/T)
    A = log(kb) - B*(1/T - 1/Tᵣ)

    #11. calculate ΔhV,ΔhL
    ΔhV = enthalpy(model,p,T,y,phase = :v)
    ΔhL = enthalpy(model,p,T,x,phase = :l)

    #12. calculate ΔhVᵣ,ΔhLᵣ
    ΔhVᵣ = enthalpy(model,p,Tᵣ,y,phase = :v)
    ΔhLᵣ = enthalpy(model,p,Tᵣ,x,phase = :l)

    #13. set C = ΔhVᵣ
    C = ΔhVᵣ

    #14. calculate D from 31)
    D = (ΔhV - C)/(T - Tᵣ)

    #15. set E = ΔhLᵣ
    E = ΔhLᵣ

    #16. calculate F from 32)
    F = (ΔhL - E)/(T - Tᵣ)

    𝕧 = u
    append!(𝕧,(A,B,C,D,E,F))
    return 𝕧,T̃,Tᵣ,kb₀,kb
end

function pq_flash_x0_step_2(model,p,Q,z,pure,T00 = nothing)
    if T00 != nothing
        #solve just for alpha
        K = suggest_K(model,p,T00,z,pure)
        α = rachfordrice(K,z)
        T = T0*one(α)
        T,α
    end


    sat = extended_saturation_temperature.(pure,p,crit_retry = false)
    _crit = __crit_pure.(sat,pure)
    fix_sat_ti!(sat,pure,_crit,p)
    T0 = sum(z[i]*sat_t[i][1] for i in 1:length(z))

    sat0 = extended_saturation_temperature.(pure,T0,_crit)
    dPdTsat = __dlnPdTinvsat.(pure,sat0,_crit,T0,true,false)
    
    #=

    P0:
    dlnpdTinv,logp0,Tcinv = __dlnPdTinvsat(pure,sat,crit,T,volatile,false)
    lnp = logp0 + dlnpdTinv*(1/T - Tcinv)
    =#
    T0 = sum(z[i]*sat_t[i][1] for i in 1:length(z))
    Kx = fill(T0,length(model))
    xx = fill(T0,length(model))
    yx = fill(T0,length(model))
    Tsearch = first.(sat)
    push!(Tsearch,T0)


    vf0 = volume(model,p,T0,z)
    hf0 = hf = VT_enthalpy(model,vf0,T0,z)
    cpf = VT_isobaric_heat_capacity(model,vf0,T0,z)
    function Q_balance(Tx)
        Kx .= __sat_p_aprox(dPdTsat,Tx)
        αx = rachfordrice(Kx,z)
        x = Clapeyron.rr_flash_liquid!(xx,Kx,z,αx)
        y = Clapeyron.rr_flash_vapor!(yx,Kx,z,αx)
        
        hf = hf0 + cpf*(Tx - T0)
        #todo: replace those with an aproximation.
        hl = enthalpy(model,p,Tx,x,phase = :l)
        hv = enthalpy(model,p,Tx,y,phase = :v)
        isnan(hv) && iszero(αx) && (hv = zero(hv))
        isnan(hl) && isone(αx) && (hl = zero(hl))
        L,V,F = (1-αx),αx,one(αx)
        return abs(V*hv + L*hl - Q - F*hf)
    end
    
    _,Tmin = findmin(Q_balance,Tsearch)
    #todo: refine this
    Kx .= __sat_p_aprox(dPdTsat,Tmin)
    α = rachfordrice(Kx,z)
    return Tmin,α
end

function __sat_p_aprox(dPdTsat,T)
    dlnpdTinv,logp0,Tcinv = dPdTsat
    lnp = logp0 + dlnpdTinv*(1/T - Tcinv)
    return exp(lnp)
end

function pq_flash(model,p,Q,z,T0 = nothing)

end  