function slle_solubility(model::CompositeModel,p,T;x0=nothing)
    solute = model.components[end]
    idx_sol = solute.==model.components

    solid_r,idx_sol_r = index_reduction(model.solid,idx_sol)
    μsol = chemical_potential(solid_r,p,T,[1.])
    
    if isnothing(x0)
        x0 = x0_slle_solubility(model,p,T,μsol)
    end
    f!(F,x) = obj_slle_solubility(F,model.liquid,p,T,[exp10(x[1]),1-exp10(x[1])-exp10(x[2]),exp10(x[2])],[exp10(x[3]),1-exp10(x[3])-exp10(x[4]),exp10(x[4])],μsol)
    results = Solvers.nlsolve(f!,x0)
    sol = exp10.(Solvers.x_sol(results))
    x1 = [sol[1],1-sol[1]-sol[2],sol[2]]
    x2 = [sol[3],1-sol[3]-sol[4],sol[4]]
    return (x1,x2)
end

function obj_slle_solubility(F,liquid,p,T,x1,x2,μsol)
    nc = length(liquid)
    γliq1 = activity_coefficient(liquid,p,T,x1)
    μliq1 = @. Rgas()*T*log(γliq1*x1)
    γliq2 = activity_coefficient(liquid,p,T,x2)
    μliq2 = @. Rgas()*T*log(γliq2*x2)

    F[1] = μliq1[end] - μsol[1]
    F[2:nc+1] = @. μliq2 - μliq1
    return F
end

function x0_slle_solubility(model,p,T,μsol)
    z0 = [0.5,0.5,1e-10]
    x0 = zeros(length(model),2)
    for i in 1:length(model)-1
        z∞ = ones(length(model))*1e-3
        z∞[i] = 1.0
        z∞ ./= sum(z∞)
        γ∞ = activity_coefficient(model.liquid,p,T,z∞)
        z∞[end] = exp(μsol[1]/(Rgas()*T)-log(γ∞[end]))
        z∞ ./= sum(z∞)
        γ∞ = activity_coefficient(model.liquid,p,T,z∞)
        x = zeros(length(model))
        x[i] = 1.0
        x[end] = exp(μsol[1]/(Rgas()*T)-log(γ∞[end]))
        x[findfirst(x.==0)] = 1/γ∞[findfirst(x.==0)]
        x ./= sum(x)

        # x = x0_lle_init(model.liquid,p,T,z0,x)

        x0[:,i] = x
    end
    println(x0)

    K0 = x0[:,1]./x0[:,2]
    z0 = (x0[:,1]+x0[:,2])/2
    println(z0)
    println(K0)

    (x,n,G) = tp_flash(model.liquid,p,T,z0,RRTPFlash(K0=K0,equilibrium=:lle))
    x0 = zeros(2*(length(model)-1))
    x0[1:2] = [x[1,1],x[1,3]]
    x0[3:end] = [x[2,1],x[2,3]]

    return log10.(x0)
end