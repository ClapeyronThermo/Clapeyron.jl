
abstract type AbstractTPInterpolation <: SatPureAproximation end


"""
    invTlogPInterpolation(max_length::Int = 256)

performs interpolation of available values of saturation equilibria.

"""
struct invTlogPInterpolation 
    max_length::Int
end

function sizehint(cachedmodel::CACHED_SAT_PURE_APROX{invTlogPInterpolation})
    return cachedmodel.saturation_pressure_aprox.max_length
end

const TPInterpolation = invTlogPInterpolation

TPInterpolation() = TPInterpolation(256) #max length of 200 data points


function x0_sat_pure(model::CACHED_SAT_PURE_APROX{invTlogPInterpolation},T)
    if !haskey(model.cache,:saturation_pressure)
        return x0_sat_pure(model.model,T)
    end
    
    T_vec,P_vec,Vl_vec,Vv_vec = model.cache[:saturation_pressure]
    isone(length(T_vec)) && return x0_sat_pure(model.model,T) #use underlying x0_sat_pure
    Tmin,Tc = first(T_vec),last(T_vec)
    nan = zero(T)/zero(T)
    T > Tc && return [nan,nan] #error, because T>Tc
    T < Tmin && return x0_sat_pure(model.model,T) #use underlying x0_sat_pure
    pos = searchsorted(T_vec,T)
    
    if length(pos) == 1 #the value is already there
        log10vl = log10(only(Vl_vec[pos]))
        log10vv = log10(only(Vv_vec[pos]))
        return [log10vl,log10vv]
    end
    adyacent_pos = last(pos):first(pos)
    T1,T2 = T_vec[adyacent_pos]
    P1,P2 = P_vec[adyacent_pos]

    #invTlogP: log(Pi)  ≈ A-B/Ti
    #log(P1) - log(P2) = -B(1/T1 - 1/T2)
    #log(P1/P2)/(1/T1 - 1/T2)
    B = -log(P1/P2)/(1/T1 - 1/T2)
    A = log(P1) + B/T1
    
    P0 = exp(A - B/T)
    vv0 = volume_virial(model.model,P0,T)

    vl1,vl2 = Vl_vec[adyacent_pos]
    #rackett approximation
    #log(V1) = C + D(1-Tr1)^(2/7) = C+D*fT1
    #log(V1) - log(V2) = D(ft1-ft2)
    fT1 = (1-T1/Tc)^(2/7)
    fT2 = (1-T2/Tc)^(2/7)
    D = log(vl1/vl2)/(fT1-fT2)
    C = log(vl1) - D*fT1
    vl0_interp = exp(C+D*(1-T/Tc)^(2/7))
    vl0_lb = x0_volume(model.model,P0,T,SA[1.0],phase=:l)
    #vl00 = 0.9*vl0_interp + 0.1*vl0_lb
    #vl0 = volume_compress(model.model,P0,T,V0=vl00)
    vl0 = vl0_interp
    return [log10(vl0),log10(vv0)]
end

function saturation_pressure_approx(model::CACHED_SAT_PURE_APROX{invTlogPInterpolation},T)
    if !haskey(model.cache,:saturation_pressure)
        return saturation_pressure(model,T)
    end
    T_vec,P_vec,Vl_vec,Vv_vec = model.cache[:saturation_pressure]
    isone(length(T_vec)) && return saturation_pressure(model,T) #use underlying x0_saturation_pressure
    Tmin,Tc = first(T_vec),last(T_vec)
    nan = zero(T)/zero(T)
    T > Tc && return (nan,nan,nan) #error, because T>Tc
    T < Tmin && return saturation_pressure(model,T) #use underlying x0_saturation_pressure
    pos = searchsorted(T_vec,T)
    
    if length(pos) == 1 #the value is already there
        vl = only(Vl_vec[pos])
        vv = only(Vv_vec[pos])
        p0 = only(P_vec[pos])
        return p0,vl,vv
    end
    adyacent_pos = last(pos):first(pos)
    T1,T2 = T_vec[adyacent_pos]
    P1,P2 = P_vec[adyacent_pos]

    #invTlogP: log(Pi)  ≈ A-B/Ti
    #log(P1) - log(P2) = -B(1/T1 - 1/T2)
    #log(P1/P2)/(1/T1 - 1/T2)
    B = -log(P1/P2)/(1/T1 - 1/T2)
    A = log(P1) + B/T1
    
    P0 = exp(A - B/T)
    vv0 = volume_virial(model.model,P0,T)

    vl1,vl2 = Vl_vec[adyacent_pos]
    #rackett approximation
    #log(V1) = C + D(1-Tr1)^(2/7) = C+D*fT1
    #log(V1) - log(V2) = D(ft1-ft2)
    fT1 = (1-T1/Tc)^(2/7)
    fT2 = (1-T2/Tc)^(2/7)
    D = log(vl1/vl2)/(fT1-fT2)
    C = log(vl1) - D*fT1
    vl0_interp = exp(C+D*(1-T/Tc)^(2/7))
    vl0_lb = x0_volume(model.model,P0,T,SA[1.0],phase=:l)
    vl0 = vl0_interp 
    #vl0 = volume_compress(model.model,P0,T,V0=vl00)
    return (P0,vl0,vv0)
end
function saturation_pressure_cached(::TPInterpolation,model::CachedEoS,T)
    res = saturation_pressure!(model,T)
end

export saturation_pressure_approx, TPInterpolation

