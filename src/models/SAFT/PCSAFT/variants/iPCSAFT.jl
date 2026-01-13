abstract type iPCSAFTModel <: PCSAFTModel end

struct iPCSAFTParam{T} <: EoSParam
    Mw::SingleParam{T}
    segment::SingleParam{T}
    v_shift::SingleParam{T}
    sigma::PairParam{T}
    epsilon::PairParam{T}
    epsilon_assoc::AssocParam{T}
    bondvol::AssocParam{T}
end

function iPCSAFTParam(Mw,segment,c,sigma,epsilon,epsilon_assoc,bondvol)
    return build_parametric_param(iPCSAFTParam,Mw,segment,c,sigma,epsilon,epsilon_assoc,bondvol)
end

Base.eltype(p::iPCSAFTParam{T}) where T = T

@newmodel iPCSAFT iPCSAFTModel iPCSAFTParam{T}

"""
    PCSAFTModel <: SAFTModel

    PCSAFT(components;
    idealmodel = BasicIdeal,
    userlocations = String[],
    ideal_userlocations = String[],
    reference_state = nothing,
    verbose = false,
    assoc_options = AssocOptions())

## Input parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `v_shift`: Single Parameter (`Float64`) - Volume shift `[m³·mol⁻¹]`
- `sigma`: Single Parameter (`Float64`) - Segment Diameter `[Å]`
- `epsilon`: Single Parameter (`Float64`) - Reduced dispersion energy `[K]`
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Parameter (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m³]`

## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `v_shift`: Single Parameter (`Float64`) - Volume shift `[m³·mol⁻¹]`
- `sigma`: Pair Parameter (`Float64`) - Mixed segment Diameter `[m]`
- `epsilon`: Pair Parameter (`Float64`) - Mixed reduced dispersion energy`[K]`
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m³]`

## Input models
- `idealmodel`: Ideal Model

## Description
Perturbed-Chain SAFT (PC-SAFT), with parameters fitted to reproduce the critical point and the saturated liquid density at T = 0.8Tc

## References
1.  Moine, E., Piña-Martinez, A., Jaubert, J.-N., Sirjean, B., & Privat, R. (2019). I-PC-SAFT: An Industrialized Version of the Volume-Translated PC-SAFT Equation of State for Pure Components, Resulting from Experience Acquired All through the Years on the Parameterization of SAFT-Type and Cubic Models. Industrial & Engineering Chemistry Research, 58(45), 20815–20827. [10.1021/acs.iecr.9b04660](https://doi.org/10.1021/acs.iecr.9b04660)
"""
PCSAFT
export iPCSAFT

default_references(::Type{iPCSAFT}) = ["10.1021/acs.iecr.9b04660"]
default_locations(::Type{iPCSAFT}) = ["SAFT/PCSAFT/iPCSAFT/iPCSAFT_like.csv"]

function lb_volume(model::iPCSAFTModel, T, z)
    c = model.params.v_shift.values
    c̄ = dot(z, c)
    lb_v0 = lb_volume_saft(model, T, z)
    return lb_v0 + c̄
end

function x0_volume_liquid(model::iPCSAFTModel,p,T,z)
    lb_v = lb_volume(model,T,z)
    Ts = T_scale(model,z)
    if T > 0.9Ts
        return 1.25*lb_v
    else
        return x0_volume_liquid_lowT(model,p,T,z)
    end
end

function transform_params(::Type{iPCSAFT},params)
    sigma = params["sigma"]
    sigma.values .*= 1E-10
    return saft_lorentz_berthelot(params)
end

function a_res(model::iPCSAFTModel, V0, T, z)
    c̄ = dot(z,model.params.v_shift.values)
    V = V0 + c̄
    _data = @f(data)
    #=
    the a_v_shift term accounts for translation in the ideal model
    #Δideal = log(V0/V) = -log(V/V0) = -log(1 + c̄/V0)

    all translated models should have this term
    i'm almost sure cubics have it if we rewrite the a_res expression
    =#
    a_v_shift = -log1p(c̄/V0)
    return @f(a_hc,_data) + @f(a_disp,_data) + @f(a_assoc,_data) + a_v_shift
end

#optimization
function  Δ(model::iPCSAFT, V, T, z,_data=@f(data))
    ϵ_assoc = model.params.epsilon_assoc.values
    κ = model.params.bondvol.values
    σ = model.params.sigma.values
    Δout = assoc_similar(κ,typeof(V+T+first(z)+one(eltype(model))))
    Δout.values .= false  #fill with zeros, maybe it is not necessary?
    for (idx,(i,j),(a,b)) in indices(Δout)
        κijab = κ[idx]
        if κijab != 0
            gij = @f(g_hs,i,j,_data)
            Δout[idx] = gij*σ[i,j]^3*(expm1(ϵ_assoc[i,j][a,b]/T))*κ[idx]
        end
    end
    return Δout
end
