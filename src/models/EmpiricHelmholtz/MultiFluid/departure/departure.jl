abstract type MultiFluidDepartureModel <: EoSModel end

#returns ∑xᵢ*arᵢ(δ,τ)
function multiparameter_a_res0(model,V,T,z,δ,τ,lnδ = log(δ),lnτ = log(τ),∑z = sum(z))
    pures = model.pures
    aᵣ = zero(Base.promote_eltype(model,V,T,z))
    for i in @comps
        mᵢ = pures[i]
        aᵣ += z[i]*reduced_a_res(mᵢ,δ,τ,lnδ,lnτ)
    end
    return aᵣ/∑z
end

recombine_departure!(model::MultiFluid,mixing::MultiFluidDepartureModel) = nothing

include("GEDeparture.jl")
include("EmpiricDeparture.jl")
include("QuadraticDeparture.jl")
include("TL.jl")