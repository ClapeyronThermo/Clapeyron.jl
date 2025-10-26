"""
    eutectic_point(model::CompositeModel, p = 1e5; x0 = x0_eutectic_point(model,p))
    eutectic_point(model::SolidHFus, p = 1e5; x0 = x0_eutectic_point(model,p))

Calculates the eutectic point of a binary mixture at a given pressure.

# Arguments
- `model`: A `CompositeModel` (with solid and liquid models specified) or a `SolidHFusModel`
- `p`: Pressure in `[Pa]` (default: `1e5`)
- `x0`: Initial guess as a vector or tuple `[T, x₁]` containing the temperature (`[K]`) and mole fraction of the first component (default: `x0_eutectic_point(model, p)`)

# Returns
A tuple `(T_eutectic, x_eutectic)` containing:
- Eutectic temperature in `[K]`
- Liquid composition at the eutectic point as a vector `[x₁, x₂]`

"""
function eutectic_point(model::CompositeModel,p=1e5; x0=nothing)
    p = p*one(eltype(model))
    binary_component_check(eutectic_point,model)
    solid = solid_model(model)
    fluid = fluid_model(model)
    f(x) = obj_eutectic_point(solid,fluid,p,1/x[1],FractionVector(exp(x[2])))
    if x0 === nothing 
        T0,x10 = x0_eutectic_point(model,p)
    else
        T0,x10 = x0
    end
    _1 = Base.promote_eltype(fluid_model(model),p) |> one
    v0 = svec2(1/T0,log(x10),_1)
    sol  = Solvers.nlsolve2(f,v0,Solvers.Newton2Var())
    #sol = Solvers.x_sol(r)
    #!all(<(r.options.f_abstol),r.info.best_residual) && (sol .= NaN)
    T = 1/sol[1]
    x = FractionVector(exp(sol[2]))
    return T,x
end

function obj_eutectic_point(solid,liquid,p,T,x)
    μsol = chemical_potential(solid,p,T,x)
    γliq = activity_coefficient(liquid,p,T,x; phase=:l)
    RT = Rgas()*T
    dμ1 = log(γliq[1]*x[1]) - μsol[1]/RT
    dμ2 = log(γliq[2]*x[2]) - μsol[2]/RT
    return SVector(dμ1,dμ2)
    #F[2] = dμ2
    #F[1:2] = μliq .- μsol
    #return F
end

function x0_eutectic_point(model::CompositeModel,p)
    solid = solid_model(model)
    fluid = fluid_model(model)
    R = Rgas(fluid)

    if solid isa SolidHfusModel
        Tm,Hfus = solid.params.Tm.values,solid.params.Hfus.values
        T1,T2 = Tm[1],Tm[2]
        K1,K2 = -Hfus[1]/R,-Hfus[2]/R
    else
        pure = split_pure_model(model)
        m1,m2 = pure[1],pure[2]
        T1,vs1,vl1 = melting_temperature(m1,p)
        T2,vs2,vl2 = melting_temperature(m2,p)
        K1 = -dpdT_saturation(m1,solid,vl1,vs1,T1)*T1*T1/p
        K2 = -dpdT_saturation(m2,solid,vl2,vs2,T2)*T2*T2/p
    end

    #account for non-idealities due to activity coefficients
    #taylor expansion of log(γ1) around xi = 1, T = Tm
    #=
    ϵx = sqrt(eps())
    ϵT1,ϵT2 = sqrt(eps(T1)),sqrt(eps(T2))
    γ1ϵ = activity_coefficient(fluid,p,T1 - ϵT1,FractionVector(1.0 - ϵx,2); phase=:l)[1]
    γ2ϵ = activity_coefficient(fluid,p,T2 - ϵT2,FractionVector(1.0 - ϵx,1); phase=:l)[2]
    ∂γ1γT1,∂γ2γT2 = (1 - γ1ϵ)/ϵT1,(1 - γ2ϵ)/ϵT2
    ∂γ1γx1,∂γ2γx2 = (1 - γ1ϵ)/ϵx,(1 - γ2ϵ)/ϵx
    ∂logγ1γw1,∂logγ2γw2 = ∂γ1γx1,∂γ2γx2
    ∂logγ1γτ1,∂logγ2γτ2 = -∂γ1γT1*T1*T1,-∂γ2γT2*T2*T2
    @show γ1ϵ,γ2ϵ
    @show ∂logγ1γw1,∂logγ2γw2
    K̄1 = (K1 - ∂logγ1γτ1)/(1 + ∂logγ1γw1)
    K̄2 = (K2 - ∂logγ2γτ2)/(1 + ∂logγ2γw2)
    @show K̄1,K̄2
    @show K1,K2 =#
    #return ideal_eutectic_solver(K̄1,K̄2,T1,T2)
    return ideal_eutectic_solver(K1,K2,T1,T2)
end

ideal_eutectic_solver(K1,K2,T1,T2) = ideal_eutectic_solver(promote(K1,K2,T1,T2)...)

function ideal_eutectic_solver(K1::T,K2::T,T1::T,T2::T) where T
    #solves the problem:
    #=
    log(x1) = K1*(τ - τ1)
    log(x2) = K2*(τ - τ2)
    x1 + x2 = 1
    τ(T) = 1/T

    for ideal eutectic models, K = -Hfus/R
    =#
    τ1,τ2 = 1/T1,1/T2
    f(τ) = 1 - exp(K1*(τ - τ1)) - exp(K2*(τ - τ2))
    τmin = max(τ1,τ2)
    τmax = 2*τmin
    for i in 1:100
        fx = f(τmax)
        fx > 0 && break
        τmax *= 2 
    end    
    prob = Roots.ZeroProblem(f,(τmin,τmax))
    τ = Roots.solve(prob)
    T0 = 1/τ
    x0 = exp(K1*(τ-τ1))
    return (T0,x0)
end