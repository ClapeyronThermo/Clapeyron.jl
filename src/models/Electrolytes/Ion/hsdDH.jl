abstract type hsdDHModel <: DHModel end

struct hsdDH{ϵ} <: hsdDHModel
    components::Array{String,1}
    RSPmodel::ϵ
    references::Array{String,1}
end

"""
    hsdDH(solvents::Array{String,1},
        ions::Array{String,1};
        RSPmodel = ConstRSP,
        userlocations = String[],
        RSPmodel_userlocations = String[],
        verbose = false)

## Input models
- `RSPmodel`: Relative Static Permittivity Model

## Description
This function is used to create a Debye-Hückel model (using the hard-sphere diameter). The Debye-Hückel term gives the excess Helmholtz free energy to account for the electrostatic interactions between ions in solution.
`hsdDH` is the default ion model in `ePCSAFT`

## References
1. Debye, P., Huckel, E. (1923). Phys. Z. 24, 185.
"""
hsdDH

export hsdDH

function hsdDH(solvents,ions; RSPmodel=ConstRSP, userlocations=String[], RSPmodel_userlocations=String[], verbose=false)
    components = deepcopy(ions)
    prepend!(components,solvents)
    references = String[]
    init_RSPmodel = @initmodel RSPmodel(solvents,ions,userlocations = RSPmodel_userlocations, verbose = verbose)
    model = hsdDH(components, init_RSPmodel,references)
    return model
end

function hard_sphere_diameter(model::CKSAFTModel, V, T, z, _data = @f(data))
    return first(_data)
end

function hard_sphere_diameter(model::PCSAFTModel, V, T, z, _data = @f(data))
    return first(_data)
end

function hard_sphere_diameter(model::SAFTVRMieModel, V, T, z, _data = @f(data))
    return first(_data)
end

function hard_sphere_diameter(model::SAFTgammaMieModel, V, T, z, _data = @f(data))
    _,_,vrdata = _data
    return first(vrdata)
end

function hard_sphere_diameter(model::CPAModel, V, T, z, _data = @f(data))
    σ = zeros(Base.promote_eltype(model,T,z),length(model))
    for i in 1:length(model)
        zi = FillArrays.OneElement(i, length(model))
        bi = lb_volume(model.cubicmodel,T,zi)
        σ[i] = cbrt((3/2/N_A/π)*bi)
    end
    return σ
end


function get_sigma(ionmodel::hsdDHModel, V, T, z, model, neutral_data = @f(data))
    return hard_sphere_diameter(model,V,T,z,neutral_data)
end