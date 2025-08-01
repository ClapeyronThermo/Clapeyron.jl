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

- `Tc`: Single Parameter (Float64) - Critical Temperature `[K]`
- `Pc`: Single Parameter (Float64) - Critical Pressure `[Pa]`
- `Vc`: Single Parameter (`Float64`) - Critical Volume `[m³·mol⁻¹]`

## Model Parameters

- `Tc`: Single Parameter (Float64) - Critical Temperature `[K]`
- `Pc`: Single Parameter (Float64) - Critical Pressure [Pa]
- `Zc`: Single Parameter (Float64) - Critical Compressibility Factor

## Description

Rackett Equation of State for saturated liquids. It is independent of the pressure.
```
Tr = T/Tc
V = (R̄Tc/Pc)Zc^(1+(1-Tr)^(2/7))
```

## Model Construction Examples
```julia
# Using the default database
model = RackettLiquid("water") #single input
model = RackettLiquid(["water","ethanol"]) #multiple components

# User-provided parameters, passing files or folders
model = RackettLiquid(["neon","hydrogen"]; userlocations = ["path/to/my/db","critical.csv"])

# User-provided parameters, passing parameters directly

model = RackettLiquid(["neon","hydrogen"];
        userlocations = (;Tc = [44.492,33.19],
                        Vc = [4.25e-5, 6.43e-5],
                        Pc = [2679000, 1296400])
                    )
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

function volume_impl(model::RackettLiquidModel,p,T,z,phase,threaded,vol0)
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

function volume_impl(model::RackettLiquidModel,p,T,z::SingleComp,phase,threaded,vol0)
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
                verbose::Bool=false)::RackettLiquid

## Input parameters

- `Tc`: Single Parameter (Float64) - Critical Temperature `[K]`
- `Pc`: Single Parameter (Float64) - Critical Pressure `[Pa]`
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
It can be used as a substitute of `RackettLiquid` when `Vc` is not known.

## Model Construction Examples
```julia
# Using the default database
model = YamadaGunnLiquid("water") #single input
model = YamadaGunnLiquid(["water","ethanol"]) #multiple components

# User-provided parameters, passing files or folders
model = YamadaGunnLiquid(["neon","hydrogen"]; userlocations = ["path/to/my/db","critical.csv"])

# User-provided parameters, passing parameters directly

model = YamadaGunnLiquid(["neon","hydrogen"];
        userlocations = (;Tc = [44.492,33.19],
                        Pc = [2679000, 1296400],
                        acentricfactor = [-0.03,-0.21])
                    )
```

## References
- Rackett, H. G. (1970). Equation of state for saturated liquids. Journal of Chemical and Engineering Data, 15(4), 514–517. [doi:10.1021/je60047a012](https://doi.org/10.1021/je60047a012)
- Gunn, R. D., & Yamada, T. (1971). A corresponding states correlation of saturated liquid volumes. AIChE Journal. American Institute of Chemical Engineers, 17(6), 1341–1345. [doi:10.1002/aic.690170613](https://doi.org/10.1002/aic.690170613)
"""
function YamadaGunnLiquid(components; userlocations = String[], verbose::Bool=false)
    _components = format_components(components)
    params = getparams(_components, ["properties/critical.csv"]; userlocations = userlocations, verbose = verbose)
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