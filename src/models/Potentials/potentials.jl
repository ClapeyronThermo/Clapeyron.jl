
#a_potential(smodel,ρ,T̄,data)
abstract type PotentialModel end
abstract type LJPotentialModel <: PotentialModel end
abstract type MiePotentialModel <: PotentialModel end


function a_res_ecs(refmodel::PotentialModel,scalingmodel,V,T,z)
    f,h = shape_factors(scalingmodel,refmodel,T,z)
    T0 = T/f
    rho0 = sum(z)*h/V
    data = ecs_data(refmodel,scalingmodel,V,T,z)
    return a_potential(refmodel,rho0,T0,z,data)
end

function a_res(model::PotentialModel,V,T,z)
    scaling = NoECSScaling()
    return a_res_ecs(model,scaling,V,T,z)
end

Rgas(::PotentialModel) = 1.0
T_scale(::PotentialModel,z) = 1.0

#default
ecs_data(refmodel,scalingmodel,V,T,z) = nothing 
a_potential(model,rho,T) = a_potential(model,rho,T,nothing)
is_splittable(::PotentialModel) = false
idealmodel(::PotentialModel) = BasicIdeal()

include("scalings/scalings.jl")