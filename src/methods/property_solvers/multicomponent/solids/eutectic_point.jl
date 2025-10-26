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


    Te,x1e = ideal_eutectic_solver(K1,K2,T1,T2)

    #=
    successive substitution method with bounds in T
    
    1D function:
    f(τ) = 1 - exp(K1*(τ - τ1))/γ1(τ,x1) - exp(K2*(τ - τ2))/γ2(τ,x1)
    
    0. start with Te,x1e from ideal eutectic solution -> calculate γ1,γ2
    
    at each step:
        1. τe_ss = (log((1 - x1e)*γ2)/K2 + τ2)
        2. τe update: check if τe_ss is in bounds, if τe_ss is outside bounds, use bisection step
        3. x1e = exp(K1*(τe - τ1))/γ1
        4. update γ1,γ2
        5. check bounds to see if we break the iteration.
    
    
    For SolidHfus, this is equivalent to the 2d eutectic solver
    =#



    γ1,γ2 = activity_coefficient(fluid,p,Te,FractionVector(x1e,2); phase=:l)
    τe = 1/Te
    τ1,τ2 = 1/T1,1/T2
    τmin = max(τ1,τ2)*one(τe+γ1)
    τmax = Inf*one(τe+γ1)
    fτ = 1 - exp(K1*(τe - τ1))/γ1 - exp(K2*(τe - τ2))/γ2
    if fτ > 0
        τmax = τe
    else
        τmin = τe
    end

    for i in 1:20
        τex = (log((1 - x1e)*γ2)/K2 + τ2)

        if τex < τmin
            τe = 0.5*τe + 0.5*τmin 
        elseif τex > τmax
            τe = 0.5*τe + 0.5*τmax
        else
            τe = τex
        end

        x1e = exp(K1*(τe - τ1))/γ1
        Te = 1/τe
        γ1,γ2 = activity_coefficient(fluid,p,Te,FractionVector(x1e,2); phase=:l)
        fτ = 1 - exp(K1*(τe - τ1))/γ1 - exp(K2*(τe - τ2))/γ2
        if fτ > 0
            τmax = τe
        else
            τmin = τe
        end
        abs(τmax-τmin)/τe < 1e-4 && break
    end
    return Te,x1e
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