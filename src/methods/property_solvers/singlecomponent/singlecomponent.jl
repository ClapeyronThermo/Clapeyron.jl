"""
    check_valid_sat_pure(model,P_sat,Vl,Vv,T,z = SA[1.0])

Checks that a saturation method converged correctly. It checks:
- That both volumes are mechanically stable.
- That both volumes are different, with a difference of at least `ε0` epsilons.
"""
function check_valid_sat_pure(model,P_sat,V_l,V_v,T,z = SA[1.0])
   return check_valid_eq2(model,model,P_sat,V_l,V_v,T,z)
end

_p∂p∂V(model,V,T,z,p) = p∂p∂V(model,V,T,z)

function _p∂p∂V(model::GibbsBasedModel,V,T,z,p)
    _,dvdp = V∂V∂p(model,p,T,z)
    return p,1/dvdp
end

_is_positive(x::Number) = isfinite(x) && x > zero(x)
_is_positive(x::Tuple) = all(_is_positive,x)

function check_valid_eq2(model1,model2,p,V1,V2,T,z = SA[1.0],ε0 = 5e7)
    ε = abs(V1-V2)/(eps(typeof(V1-V2)))
    ε <= ε0 && return false
    p1,dpdv1 = _p∂p∂V(model1,V1,T,z,p)
    p2,dpdv2 = _p∂p∂V(model2,V2,T,z,p)
    return  (dpdv1 <= 0)                    && #mechanical stability of phase 1
            (dpdv2 <= 0)                    && #mechanical stability of phase 2
            _is_positive((p1,p2,V2,V2,T,p)) #positive and finite pressures and volumes
end

function a∂a∂V(model,V,T,z::AbstractVector)
    f = @deferred_V(a_res,a∂a∂V)
    a,∂a∂V = Solvers.f∂f(f,V)
    return SVector(a,∂a∂V)
end

function μp_equality1_p(model1,model2,v1,v2,T,ps,μs,z = SA[1.0])
    RT = Rgas(model1)*T
    A1,Av1 = a∂a∂V(model1,v1,T,z)
    A2,Av2 = a∂a∂V(model2,v2,T,z)
    p1,p2 = RT*(-Av1 + 1/v1),RT*(-Av2 + 1/v2)
    Δμᵣ = A1 - v1*Av1 - A2 + v2*Av2 + log(v2/v1)
    Fμ = Δμᵣ
    Fp = (p1 - p2)*ps
    return SVector(Fμ,Fp)
end

function μp_equality1_p(model,v1,v2,T,z = SA[1.0]) 
    ps,μs = equilibria_scale(model,z)
    μp_equality1_p(model,model,v1,v2,T,ps,μs,z)
end

function μp_equality1_T(model1,model2,v1,v2,p,T,ps,μs,z = SA[1.0])
    RT = Rgas(model1)*T
    A1,Av1 = a∂a∂V(model1,v1,T,z)
    A2,Av2 = a∂a∂V(model2,v2,T,z)
    p1,p2 = RT*(-Av1 + 1/v1),RT*(-Av2 + 1/v2)
    Δμᵣ = A1 - v1*Av1 - A2 + v2*Av2 + log(v2/v1)
    Fμ = Δμᵣ
    Fp1 = (p1 - p)*ps
    Fp2 = (p2 - p)*ps
    return SVector(Fμ,Fp1,Fp2)
end

#variant used for edge_temperature with multicomponent.
#somehow is more stable
function μp_equality1_T2(model,p,z,x,Ts)
    lnv1,lnv2,T1,T2 = x
    n = sum(z)
    v1,v2 = exp(lnv1),exp(lnv2)
    RT1,RT2 = n*Rgas(model)*T1,n*Rgas(model)*T2
    A1,Av1 = a∂a∂V(model,V1,T1,z)
    A2,Av2 = a∂a∂V(model,V2,T2,z)
    p1,p2 = RT1*(-Av1 + 1/v1),RT2*(-Av2 + 1/v2)
    Δμᵣ = A1 - v1*Av1 - A2 + v2*Av2 + log(v2/v1)
    Fμ = Δμᵣ
    Fp1 = (p1 - p)/p
    Fp2 = (p2 - p)/p
    FT = (T1 - T2)/Ts
    return SVector(Fμ,Fp1,Fp2,FT)
end
struct μequality1_obj{TYPE,M1,M2,TT,PS,MUS,Z}
    type::Val{TYPE}
    model1::M1
    model2::M2
    X::TT
    ps::PS
    μs::MUS
    z::Z
end


function edge_pressure_objective(model1,model2,T,ps,μs,z = SA[1.0])
    return μequality1_obj(Val(:Psat),model1,model2,T,ps,μs,z)
end

function edge_temperature_objective(model1,model2,T,ps,μs,z = SA[1.0])
    return μequality1_obj(Val(:Tsat),model1,model2,T,ps,μs,z)
end

function edge_temperature_objective2(model1,model2,T,ps,μs,z = SA[1.0])
    return μequality1_obj(Val(:Tsat2),model1,model2,T,ps,μs,z)
end

function (obj::μequality1_obj{:Psat})(x)
    lnv1,lnv2 = x
    v1,v2 = exp(lnv1),exp(lnv2)
    return μp_equality1_p(obj.model1,obj.model2,v1,v2,obj.X,obj.ps,obj.μs,obj.z)
end

function (obj::μequality1_obj{:Tsat})(x)
    T,lnv1,lnv2 = x
    v1,v2 = exp(lnv1),exp(lnv2)
    return μp_equality1_T(obj.model1,obj.model2,v1,v2,obj.X,T,obj.ps,obj.μs,obj.z)
end

function (obj::μequality1_obj{:Tsat2})(x)
    return μp_equality1_T2(obj.model,obj.X,obj.z,x,obj.μs)
end

StaticForwardDiffTags.inner_function(f::μequality1_obj{:Psat}) = μp_equality1_p
StaticForwardDiffTags.inner_function(f::μequality1_obj{:Tsat}) = μp_equality1_T
StaticForwardDiffTags.inner_function(f::μequality1_obj{:Tsat2}) = μp_equality1_T2
StaticForwardDiffTags.deferred_valtype(f::μequality1_obj{TYPE,M1,M2,TT,PS,MUS,Z}) where {TYPE,M1,M2,TT,PS,MUS,Z} = Base.promote_eltype(f.model1,f.model2,f.X,f.ps,f.μs,f.z)

function try_2ph_edge_pressure(model1,model2,T,v10,v20,ps,μs,z,method)
    f = WithContext(edge_pressure_objective(model1,model2,T,ps,μs,z),∂Tag{∂₁f}())
    TT = T*oneunit(eltype(model1))*oneunit(eltype(model2))
    V0 = svec2(log(v10),log(v20),TT)
    if !_is_positive((v10,v20,T))
        _0 = zero(V0[1])
        nan = _0/_0
        fail = (nan,nan,nan)
        return fail,false
    end
    neq_options = method === nothing ? NEqOptions() : NEqOptions(method)
    sol = Solvers.nlsolve2(f,V0,Solvers.Newton2Var(),neq_options)
    v1 = exp(sol[1])
    v2 = exp(sol[2])
    p_eq = pressure(model2,v2,T)
    res = (p_eq,v1,v2)
    valid = check_valid_eq2(model1,model2,p_eq,v1,v2,T,z)
    return res,valid
end

function try_2ph_pure_temperature(model1,model2,p,T0,v10,v20,ps,μs,method)
    f = WithContext(edge_temperature_objective(model1,model2,p,ps,μs,z),∂Tag{∂₁f}())
    pp = p*oneunit(eltype(model1))*oneunit(eltype(model2))
    V0 = svec3(T0,log(v10),log(v20),pp)

    if !_is_positive((v10,v20,p,T0))
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