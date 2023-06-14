
#returns ∑xᵢ*arᵢ(δ,τ)
function multiparameter_a_res0(model,V,T,z,δ,τ,∑z)
    pures = model.pures
    aᵣ = zero(δ+τ)
    lnδ = log(δ)
    lnτ = log(τ)
    for i in @comps
        mᵢ = pures[i]
        aᵣ += z[i]*reduced_a_res(mᵢ,δ,τ,lnδ,lnτ)
    end
    return aᵣ/∑z
end

"""
    calculate_missing_mixing(params::EmpiricMultiFluidParam,mixing::MixingRuleModel)

Calculates missing values, using the parameters stored in `params`. modifies `mixing` implace. this function is called at `EmpiricMultiFluid` model creation time.

"""
calculate_missing_mixing!(params,mixing) = nothing

include("Asymmetric.jl")
include("LorentzBerthelotMixing.jl")