
#a_potential(smodel,V,T̄)
abstract type PotentialModel <: EoSModel end
abstract type LJPotentialModel <: PotentialModel end
abstract type LJTSPotentialModel <: LJPotentialModel end
abstract type MiePotentialModel <: LJPotentialModel end
abstract type Mie6PotentialModel <: MiePotentialModel end

function a_res_ecs(refmodel::PotentialModel,scalingmodel,V,T,z)
    f,h = shape_factors(scalingmodel,refmodel,T,z)
    T0 = T/f
    V0 = V/(sum(z)*h)
    return a_potential(refmodel,V0,T0)
end

function a_res(model::PotentialModel,V,T,z)
    scaling = NoECSScaling()
    return a_res_ecs(model,scaling,V,T,z)
end

#default
Rgas(::PotentialModel) = 1.0
T_scale(::PotentialModel,z) = 1.0
is_splittable(::PotentialModel) = false
idealmodel(::PotentialModel) = BasicIdeal(1.0) #basic ideal with R = 1
Base.length(::PotentialModel) = 1
#defaults for LJ model

function x0_crit_pure(model::LJPotentialModel,z)
    return (1.32, -log10(0.31))
end

include("scalings/scalings.jl")
