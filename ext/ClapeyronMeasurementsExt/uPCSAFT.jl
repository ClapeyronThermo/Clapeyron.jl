using Clapeyron: PCSAFTModel, @newmodel

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