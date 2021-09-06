function vϕ(model::EoSModel,P,T,z;phase=:unknown)
    v = volume(model,P,T,z,phase = phase)
    vx = volume(model,P,T,z)

    n = sum(z)
    p = pressure(model,v,T,z)
    @show p
    lnZ = log(p*v/(n*R̄*T))
    fun(z) = a_res(model,v,T,z)*n
    ϕ = ForwardDiff.gradient(fun,z)
    return v,ϕ .- lnZ
end

function flash_tpd(xi, lnϕ,x0,lnϕ0)
    f(xiᵢ,lnϕᵢ,x0ᵢ,lnϕ0ᵢ) = xiᵢ*(log(xiᵢ) + lnϕᵢ-log(x0ᵢ)- lnϕ0ᵢ)
    return mapreduce(f,+,xi, lnϕ,x0,lnϕ0)
end

function x0_PT_flash(model::EoSModel,T,P,x0)
    #check each T with T_scale, if treshold is over, replace Pi with inf
 
 _0 = zero(P+T+first(x0))
 pure = split_model(model)
 crit = crit_pure.(pure)
 
 T_c = [tup[1] for tup in crit]
 V_c = [tup[3] for tup in crit]
 nan = _0/_0 
 sat_nan = (nan,nan,nan)
 replaceP = ifelse.(T_c .< T,true,false)

 sat = [if !replaceP[i] sat_pure(pure[i],T) else sat_nan end for i in 1:length(pure)]
 
 P_sat = [tup[1] for tup in sat]
 Ki = zeros(eltype(x0),length(x0))
 for i in 1:length(x0)
     if !replaceP[i]
         Ki[i] = P_sat[i][1]/P
     else 
         Ki[i] = pressure(pure[i],V_c[i],T)/P
     end
 end

return Ki
end 

function PT_flash(model::EoSModel,P,T,z,K0=nothing)
    
    _0 = zero(P+T+first(z))
    _1 = one(_0)
    if K0 === nothing
        K0 = x0_PT_flash(model,T,P,z)
    end

    β = rr_vle_vapor_fraction(K0,z)
    v0 = volume(model,P,T,z)
    if β <= 0  #bubble point assumption, initial phase supposed 
        β = _0
        β0 = _0
        xil = copy(z)
        xiv = rr_flash_vapor(K0,z,_0)

    elseif β >= 0 #dew point assumption
        β = _1
        β0 = _1
        xil = rr_flash_liquid(K0,z,_1)
        xiv = copy(z)
    else #two phase assumption        
        xil = rr_flash_liquid(K0,z,β)
        xiv = rr_flash_vapor(K0,z,β)
    end

    logk = similar(xil)
    k = similar(xil)
    v0,logϕ0= vϕ(model,P,T,z) 
    vv,logϕv = vϕ(model,P,T,xiv,phase=:v) 
    vl,logϕl = vϕ(model,P,T,xil,phase=:l) 
    stable_phase = false
    @show K0
    #Rachford Rice sucessive substitution, 5 iters to improve the result a bit
    for i = 1:5
        logk = logϕl - logϕv
        k .= exp.(logk)
        β = rr_vle_vapor_fraction(k,z)
        @show k,β
        if !(_0 <= β <= _1)
            stable_phase = true
            break
        end 
        xil = rr_flash_liquid(k,z,β)
        xiv = rr_flash_vapor(k,z,β)
        vv,logϕv = vϕ(model,P,T,xiv,phase=:v) 
        
        vl,logϕl = vϕ(model,P,T,xil,phase=:l) 
        if !isfinite(vl+vv)
            stable_phase = true
            break
        end
    end
    #TPD analysis
    if stable_phase == false
        tpdl = flash_tpd(xil, logϕl,z,logϕ0)
        tpdv = flash_tpd(xiv, logϕv,z,logϕ0)
        ΔG = (1 - β) * tpdl + β * tpdv
        if all(>(0),(ΔG,tpdv,tpdl))
            stable_phase = true # more detailed stability phase analisis
        end
    end
    yy = copy(xiv)
    prepend!(yy,log10.([vl,vv]))
    @show yy
    #if TPD analysis worked, pass to optimization.
    #if not, try to check if a phase could be made stable.
    return xil, xiv
end

 
