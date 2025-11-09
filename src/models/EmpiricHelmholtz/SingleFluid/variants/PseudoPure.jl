#at the implementation level, pseudo-pure fluid models are different from pure fluid models
#just by storing an additional saturation pressure ancillary and by setting pseudo_pure = true.

struct SingleFluidFlashMethod <: FlashMethod end

is_pseudo_pure(model::SingleFluid) = model.properties.pseudo_pure

