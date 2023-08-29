abstract type RackettLiquidModel <: LiquidVolumeModel end


struct RackettLiquidParam <: EoSParam
    Tc::SingleParam{Float64}
    Pc::SingleParam{Float64}
    Zc::SingleParam{Float64}
end

@newmodelsimple RackettLiquid RackettLiquidModel RackettLiquidParam

"""
    RackettLiquid(components;
    userlocations::Vector{String}=String[],
    verbose::Bool=false)

## Input parameters

- `Tc`: Single Parameter (Float64) - Critical Temperature [K]
- `Pc`: Single Parameter (Float64) - Critical Pressure [Pa]
- `Vc`: Single Parameter (`Float64`) - Critical Volume `[m³/mol]`

## Model Parameters

- `Tc`: Single Parameter (Float64) - Critical Temperature [K]
- `Pc`: Single Parameter (Float64) - Critical Pressure [Pa]
- `Zc`: Single Parameter (Float64) - Critical Compressibility Factor

## Description

Rackett Equation of State for saturated liquids. it is independent of the pressure.
```
Tr = T/Tc
V = (R̄Tc/Pc)Zc^(1+(1-Tr)^(2/7))
```
## References
- Rackett, H. G. (1970). Equation of state for saturated liquids. Journal of Chemical and Engineering Data, 15(4), 514–517. [doi:10.1021/je60047a012](https://doi.org/10.1021/je60047a012)
"""
RackettLiquid
default_locations(::Type{RackettLiquid}) = critical_data()
default_references(::Type{RackettLiquid}) = ["10.1021/je60047a012"]
function transform_params(::Type{RackettLiquid},params,components)
    Tc = params["Tc"]
    Pc = params["Pc"]
    vc = params["Vc"]
    _zc = Pc.values .* vc.values ./ (R̄ .* Tc.values)
    Zc = SingleParam("Critical Compressibility factor",components,_zc)
    params["Zc"] = Zc
    return params
end

function volume_impl(model::RackettLiquidModel,p,T,z=SA[1.0],phase=:unknown,threaded=false,vol0 = 0.0)
    tci = model.params.Tc.values
    pci = model.params.Pc.values
    zci = model.params.Zc.values
    Zcm = zero(eltype(z))
    a = zero(eltype(z))
    b = zero(eltype(z))
    ∑z = sum(z)
    checkbounds(tci,length(z))
    for i ∈ @comps
        zi = z[i]
        Tcᵢ = tci[i]
        Pcᵢ = pci[i]
        bi = (R̄*Tcᵢ)/Pcᵢ
        ai = (R̄*Tcᵢ)*bi
        zii = zi*zi
        a += zii*ai
        b += zii*bi
        Zcm += z[i]*zci[i]
        for j in 1:(i-1)
            zj = z[j]
            Tcⱼ = tci[i]
            Pcⱼ = pci[i]
            bj = (R̄*Tcⱼ)/Pcⱼ
            aj = (R̄*Tcⱼ)*bi
            aij = sqrt(ai*aj)
            bij = 0.5*(bi+bj)
            zij = zi*zj
            b += 2*zij*bij
            a += 2*zij*aij
        end
    end
    Tcm = a/b/R̄
    Pcm_inv = (b/(∑z*∑z))/(R̄*Tcm)
    Zcm = Zcm/∑z
    Tr = T/Tcm
    return ∑z*R̄*Tcm*Pcm_inv*Zcm^(1+(1-Tr)^(2/7))
end

function volume_impl(model::RackettLiquidModel,p,T,z::SingleComp,phase=:unknown,threaded=false,vol0 = 0.0)
    Tc = only(model.params.Tc.values)
    Pc = only(model.params.Pc.values)
    Pc_inv = 1/Pc
    Zc = only(model.params.Zc.values)
    ∑z = only(z)
    Tr = T/Tc
    return ∑z*R̄*Tc*Pc_inv*Zc^(1+(1-Tr)^(2/7))
end


"""
    YamadaGunnLiquid(components;
                userlocations::Vector{String}=String[],
                verbose::Bool=false)

## Input parameters

- `Tc`: Single Parameter (Float64) - Critical Temperature `[K]`
- `Pc`: Single Parameter (Float64) - Critical Pressure `[Pa]`
- `Vc`: Single Parameter (`Float64`) - Critical Volume `[m³/mol]`
- `acentricfactor`: Single Parameter (`Float64`) - Acentric Factor

## Model Parameters

- `Tc`: Single Parameter (Float64) - Critical Temperature `[K]`
- `Pc`: Single Parameter (Float64) - Critical Pressure `[Pa]`
- `Zc`: Single Parameter (Float64) - Critical Compressibility Factor

## Description

The Yamada-Gunn equation of state is a modification of the Rackett equation of state that uses a different approach to calculate the compressibility factor `Zc`:
```
Tr = T/Tc
Zc = 0.29056 - 0.08775ω
V = (R̄Tc/Pc)Zc^(1+(1-Tr)^(2/7))
```
## References
- Rackett, H. G. (1970). Equation of state for saturated liquids. Journal of Chemical and Engineering Data, 15(4), 514–517. [doi:10.1021/je60047a012](https://doi.org/10.1021/je60047a012)
- Gunn, R. D., & Yamada, T. (1971). A corresponding states correlation of saturated liquid volumes. AIChE Journal. American Institute of Chemical Engineers, 17(6), 1341–1345. [doi:10.1002/aic.690170613](https://doi.org/10.1002/aic.690170613)
"""
function YamadaGunnLiquid(components; userlocations=String[], verbose::Bool=false)
    _components = format_components(components)
    params = getparams(_components, ["properties/critical.csv"]; userlocations=userlocations, verbose=verbose)
    acentricfactor = params["acentricfactor"]
    Tc = params["Tc"]
    Pc = params["Pc"]
    _zc = 0.29056 .- 0.08775 .* acentricfactor.values
    Zc = SingleParam("Critical Compressibility factor",_components,_zc)
    packagedparams = RackettLiquidParam(Tc,Pc,Zc)
    references =["10.1021/je60047a012","10.1002/aic.690170613"]
    model = RackettLiquid(_components,packagedparams,references)
    return model
end

export YamadaGunnLiquid,RackettLiquid