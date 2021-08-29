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
    g_0 = rr_flash_eval(K0,z,_0)
    g_1 = rr_flash_eval(K0,z,_1)
    
    if g_0 <= 0  #bubble point assumption
        println("bubble point")
        β = _0
        β0 = _0
        xil = copy(z0)
        xiv = rr_flash_vapor(K0,z,_0)
    elseif g_1 >= 0 #dew point assumption
        println("dew point")
        β = _1
        β0 = _1
        xil = rr_flash_liquid(K0,z,_1)
        xiv = copy(z)
    else #two phase assumption
        println("flash point")
        β = rr_vle_vapor_fraction(K0,z)
        
        xil = rr_flash_liquid(K0,z,β)
        xiv = rr_flash_vapor(K0,z,β)
    end
    logϕl = similar(xil)
    logϕv = similar(xiv)

    vv = volume(model,P,T,xiv,phase=:v)
    vl = volume(model,P,T,xil,phase=:l)
    @show vl
    @show vv
    @show β
    
    return xil, xiv
end

 
