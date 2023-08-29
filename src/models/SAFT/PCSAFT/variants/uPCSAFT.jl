struct uPCSAFTParam <: EoSParam
    Mw::SingleParam
    segment::SingleParam
    sigma::PairParam
    epsilon::PairParam
    epsilon_assoc::AssocParam
    bondvol::AssocParam
end

abstract type uPCSAFTModel <: PCSAFTModel end
@newmodel uPCSAFT uPCSAFTModel uPCSAFTParam

export uPCSAFT

function d(model::uPCSAFTModel, V, T, z)
    ϵᵢᵢ = diagvalues(model.params.epsilon)
    σᵢᵢ = diagvalues(model.params.sigma)
    return @. σᵢᵢ*(1 - 0.12*exp(-3ϵᵢᵢ/ T))
end