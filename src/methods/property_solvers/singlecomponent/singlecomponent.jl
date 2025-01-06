"""
    check_valid_sat_pure(model,P_sat,Vl,Vv,T,ε0 = 5e7)

Checks that a saturation method converged correctly. it checks:
- That both volumes are mechanically stable
- That both volumes are different, with a difference of at least `ε0` epsilons
"""
function check_valid_sat_pure(model,P_sat,V_l,V_v,T,ε0 = 5e7)
   return check_valid_eq2(model,model,P_sat,V_l,V_v,T,ε0)
end

function check_valid_eq2(model1,model2,p,V1,V2,T,ε0 = 5e7)
    ε = abs(V1-V2)/(eps(typeof(V1-V2)))
    ε <= ε0 && return false
    p1,dpdv1 = p∂p∂V(model1,V1,T,SA[1.0])
    p2,dpdv2 = p∂p∂V(model2,V2,T,SA[1.0])
    return  (dpdv1 <= 0)        && #mechanical stability of phase 1
            (dpdv2 <= 0)        && #mechanical stability of phase 2
            T > zero(T)         && #positive temperature
            p > zero(p)         && #positive pressure
            p1 > zero(p1)       && #positive pressure at phase 1
            p2 > zero(p2)          #positive pressure at phase 2
end

function check_valid_2ph_input(v1,v2,p,T)
    isfinite(v1) && isfinite(v2) | isfinite(T) | isfinite(p) | (T < zero(T)) | (p < zero(p))
end

function μp_equality1_p(model1,model2,v1,v2,T,ps,μs)
    z = SA[1.0]
    RT = Rgas(model1)*T
    f1(V) = a_res(model1,V,T,z)
    f2(V) = a_res(model2,V,T,z)
    A1,Av1 = Solvers.f∂f(f1,v1)
    A2,Av2 =Solvers.f∂f(f2,v2)
    p1,p2 = RT*(-Av1 + 1/v1),RT*(-Av2 + 1/v2)
    Δμᵣ = A1 - v1*Av1 - A2 + v2*Av2 + log(v2/v1)
    Fμ = Δμᵣ
    Fp = (p1 - p2)*ps
    return SVector(Fμ,Fp)
end

function μp_equality1_T(model1,model2,v1,v2,p,T,ps,μs)
    z = SA[1.0]
    RT = Rgas(model1)*T
    f1(V) = a_res(model1,V,T,z)
    f2(V) = a_res(model2,V,T,z)
    A1,Av1 = Solvers.f∂f(f1,v1)
    A2,Av2 =Solvers.f∂f(f2,v2)
    p1,p2 = RT*(-Av1 + 1/v1),RT*(-Av2 + 1/v2)
    Δμᵣ = A1 - v1*Av1 - A2 + v2*Av2 + log(v2/v1)
    Fμ = Δμᵣ
    Fp1 = (p1 - p)*ps
    Fp2 = (p2 - p)*ps
    return SVector(Fμ,Fp1,Fp2)
end

function try_2ph_pure_pressure(model,T,v10,v20,ps,mus,method)
    return try_2ph_pure_pressure(model,model,T,v10,v20,ps,mus,method)
end

function try_2ph_pure_pressure(model1,model2,T,v10,v20,ps,mus,method)
    f(x) = μp_equality1_p(model1,model2,exp(x[1]),exp(x[2]),T,ps,mus)
    TT = T*oneunit(eltype(model1))*oneunit(eltype(model2))
    V0 = svec2(log(v10),log(v20),TT)

    if !check_valid_2ph_input(v10,v20,true,T)
        _0 = zero(V0[1])
        nan = _0/_0
        fail = (nan,nan,nan)
        return fail,false
    end

    sol = Solvers.nlsolve2(f,V0,Solvers.Newton2Var(),NEqOptions(method))
    v1 = exp(sol[1])
    v2 = exp(sol[2])
    p_eq = pressure(model2,v2,T)
    res = (p_eq,v1,v2)
    valid = check_valid_eq2(model1,model2,p_eq,v1,v2,T)
    return res,valid
end


function try_2ph_pure_temperature(model1,model2,p,T0,v10,v20,ps,mus,method)
    f(x) = μp_equality1_T(model1,model2,exp(x[2]),exp(x[3]),p,x[1],ps,mus)
    pp = p*oneunit(eltype(model1))*oneunit(eltype(model2))
    V0 = svec3(T0,log(v10),log(v20),pp)

    if !check_valid_2ph_input(v10,v20,p,T0)
        _0 = zero(V0[1])
        nan = _0/_0
        fail = (nan,nan,nan)
        return fail,false
    end

    if !isfinite(V0[2]) | !isfinite(V0[2]) | !isfinite(p) | (p < zero(p)) | (T0 < zero(T0))
        _0 = zero(V0[1])
        nan = _0/_0
        fail = (nan,nan,nan)
        return fail,false
    end

    sol = Solvers.nlsolve2(f,V0,Solvers.Newton2Var(),NEqOptions(method))
    T_eq = sol[1]
    v1 = exp(sol[2])
    v2 = exp(sol[3])
    res = (T_eq,v1,v2)
    valid = check_valid_eq2(model1,model2,p,v1,v2,T_eq)
    return res,valid
end

function try_2ph_pure_temperature(model,p,T0,v10,v20,ps,mus,method)
    return try_2ph_pure_temperature(model,model,p,T0,v10,v20,ps,mus,method)
end

include("saturation/saturation.jl")
include("crit_pure.jl")
include("triple_point.jl")
include("sublimation.jl")
include("melting.jl")
include("widom.jl")

export saturation_pressure, saturation_liquid_density, saturation_temperature
export crit_pure, enthalpy_vap, acentric_factor
export triple_point, sublimation_pressure, melting_pressure, sublimation_temperature, melting_temperature
export widom_pressure, widom_temperature
export ciic_pressure, ciic_temperature