abstract type CPCAMixingModel <: MixingRule end

struct CPCAMixingParam <:EoSParam
    kb1::SingleParam{Float64}
    kb2::SingleParam{Float64}
    segment::SingleParam{Float64}
end

struct CPCAMixing <: CPCAMixingModel
    components::Array{String,1}
    params::CPCAMixingParam
    references::Array{String,1}
end

function CPCAMixing(components; activity = nothing, userlocations=String[],activity_userlocations=String[], verbose::Bool=false)
    params = getparams(components,String[], userlocations=userlocations, verbose=verbose)
    k1 = params["kb1"]
    k2 = params["kb2"]
    segment = params["segment"]
    param = CPCAMixingParam(kb1,kb2,segment)
    return CPCAMixing(format_components(components), param, ["10.1021/acs.iecr.2c00902"])
end

function mixing_rule(model::ABCubicModel,V,T,z,mixing_model::CPCAMixing,α,a,b,c)
    kb1 = mixing_model.params.kb1.values
    kb2 = mixing_model.params.kb2.values
    m = mixing_model.params.segment.values
    
    Tc = model.params.Tc.values
    m̄ = zero(first(z))
    n = sum(z)
    invn = (one(n)/n)
    invn2 = invn^2
    #b̄ = dot(z,Symmetric(b),z) * invn2
    m̄²ā = zero(T+first(z))
    m̄b̄ = zero(first(z))
    for i in 1:length(z)
        zi,mi,αi = z[i],m[i],α[i]
        αi = α[i]
        m̄ += zi*mi
        βi = kb1[i]*exp(-kb2[i]*T/Tc[i])
        m̄b̄ += b[i,i]*zi*mi*βi
        m̄²ā += mi*mi*a[i,i]*αi*zi^2
        for j in 1:(i-1)
            zj,mj,αj = z[j],m[j],α[j]
            m̄²ā += 2*a[i,j]*sqrt(αi*α[j])*zi*zj*mi*mj
        end
    end
    m̄inv = 1/m̄
    ā = m̄²ā*m̄inv*m̄inv
    b̄ = m̄b̄*m̄inv
    c̄ = dot(z,c)*invn

    return ā,b̄,c̄
end

const CPCACubic = RK{I,CPAAlpha,T,CPCAMixing} where {I,T}

function cubic_ab(model::CPCACubic,V,T,z=SA[1.0],n=sum(z))
    a = model.params.a.values
    b = model.params.b.values
    T = T*float(one(T))
    α = @f(α_function,model.alpha)
    c = @f(translation,model.translation)
    length(z) > 1 && 
    return @f(mixing_rule,model.mixing,α,a,b,c)
    ā = a[1,1]*α[1]
    kb1 = model.mixing.params.kb1.values[1]
    kb2 = model.mixing.params.kb2.values[1]
    m = model.mixing.params.segment.values[1]
    Tc = model.params.Tc.values[1]
    β = kb1*exp(-kb2*T/Tc)
    b̄ = b[1,1]*β
    c̄ = c[1]
    return ā ,b̄, c̄
end

function lb_volume(model::CPCACubic,z=SA[1.0])
    mix = model.mixing.params
    m = mix.segment.values
    kb1 = mix.kb1.values
    kb2 = mix.kb2.values
    b = model.params.b.values
    res = zero(first(z) + 1.0)
    for i in @comps
        kb1i,kb2i = kb1[i],kb2[i]
        #=
        this is an exponential, so it diverges at high T
        βi = kb1*exp(-kb2*T/Tc)
        kb1 is > 1, and increasing with Mw
        kb2 < 1, increases with Mw until crossing to > 1
        as a result, βi
        if we do:
        ```  
        f1(x) = 0.7112*x^0.1502
        f2(x) = 0.1549*log(x) - 0.3224
        ff(Mw,Tr) = f1(Mw)*exp(-f2(Mw)*Tr)
        ```
        f1(x,1) ≈ 0.96, slighly decreasing but mostly stable
        f1(x,2) ≈ decreases a lot with temperature, tends to 0
        this function does not work with infinite T!!!
        we fix βi = 0.5, anything lower should fail the volume check.
        =#
        if kb1i > 1.0 && 0 < kb2i
            res += 0.5*b[i,i]*z[i]
        else
            #TODO: handle this case. useful in case of hydrogen (not covered in the correlations.)
            res += b[i,i]*z[i]
        end
    end
    return res
end

abstract type CPCAModel <: CPAModel end

struct CPCAParam <: EoSParam
    Mw::SingleParam{Float64}
    Tc::SingleParam{Float64}
    segment::SingleParam{Float64}
    a::PairParam{Float64}
    b::PairParam{Float64}
    c1::SingleParam{Float64}
    kb1::SingleParam{Float64}
    kb2::SingleParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64} 
end

struct CPCA{T <: IdealModel,c <: CubicModel} <: CPCAModel
    components::Array{String,1}
    radial_dist::Symbol
    cubicmodel::c
    params::CPCAParam
    sites::SiteParam
    idealmodel::T
    assoc_options::AssocOptions
    references::Array{String,1}
end

"""
    CPCAModel <: CPAModel

    function CPCA(components;
        idealmodel=BasicIdeal,
        translation = NoTranslation
        userlocations=String[],
        ideal_userlocations=String[],
        translation_userlocations=String[],
        verbose=false,
        assoc_options = AssocOptions())

## Input parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `a`: Single Parameter (`Float64`) - Atraction parameter `[m^6*Pa/mol]`
- `b`: Single Parameter (`Float64`) - Covolume `[m^3/mol]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Paramater (no units)
- `l`: Pair Parameter (`Float64`) (optional) - Binary Interaction Paramater (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m^3]`

## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `a`: Pair Parameter (`Float64`) - Mixed Atraction Parameter `[m^6*Pa/mol]`
- `b`: Pair Parameter (`Float64`) - Mixed Covolume `[m^3/mol]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[J]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m^3]`

## Input models
- `idealmodel`: Ideal Model
- `translation`: Translation model

## Description

Cubic + Chain + Association (CPCA) EoS. Consists in the addition of a cubic part, a chain part and an association part:
```
a_res(model::CPA) = a_res(model::Cubic) + a_chain(model) + a_assoc(model)
```

## References
1. Kontogeorgis, G. M., Michelsen, M. L., Folas, G. K., Derawi, S., von Solms, N., & Stenby, E. H. (2006). Ten years with the CPA (cubic-plus-association) equation of state. Part 1. Pure compounds and self-associating systems. Industrial & Engineering Chemistry Research, 45(14), 4855–4868. [doi:10.1021/ie051305v](https://doi.org/10.1021/ie051305v)
"""
function CPCA(components;
    idealmodel=BasicIdeal,
    translation = NoTranslation,
    userlocations=String[],
    ideal_userlocations=String[],
    translation_userlocations=String[],
    verbose=false,
    assoc_options = AssocOptions())
    locs = ["SAFT/CPA/CPCA", "properties/molarmass.csv"]
    params = getparams(components, locs; userlocations=userlocations, verbose=verbose)
    _components = format_components(components)
    
    segment = params["segment"]
    Mw  = params["Mw"]

    sites = get!(params,"sites") do
        SiteParam(components)
    end

    Pc = get!(params,"Pc") do
        SingleParam("Pc",_components,zeros(length(_components)))
    end
    k = get(params,"k",nothing)
    l = get(params,"l",nothing)
    #Tc is created from a,b
    ai = params["a"]
    bi = params["b"]
    Ωa,Ωb =  ab_consts(RK)
    ar = ai ./ Ωa
    br = bi ./ Ωb
    Tc = SingleParam("Tc",_components,@. ar/br/R̄)
    a  = epsilon_LorentzBerthelot(ai, k)
    b  = sigma_LorentzBerthelot(bi, l)

    epsilon_assoc = get!(params,"epsilon_assoc") do
        AssocParam("epsilon_assoc",components)
    end

    bondvol = get!(params,"bondvol") do
        AssocParam("bondvol",components)
    end

    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)

    c1 = get!(params,"c1") do
        SingleParam("c1",_components,@.(0.2425*Mw^0.2829))
    end
    alpha_param = CPAAlphaParam(c1)
    init_alpha = CPAAlpha(_components,alpha_param,String[]) 

    kb1 = get!(params,"kb1") do
        SingleParam("kb1",_components,@.(0.7112*(Mw)^0.1502))
    end

    kb2 = get!(params,"kb2") do
        SingleParam("kb2",_components,@.(0.1549*log(Mw) - 0.3224))
    end
    mix_param = CPCAMixingParam(kb1,kb2,segment)
    init_mixing = CPCAMixing(_components,mix_param,String[])

    init_translation = init_model(translation,components,translation_userlocations,verbose)
    cubicparams = ABCubicParam(a, b, Tc, Pc, Mw) #PR, RK, vdW
    init_cubicmodel = RK(_components,init_alpha,init_mixing,init_translation,cubicparams,init_idealmodel,String[])

    bondvol,epsilon_assoc = assoc_mix(bondvol,epsilon_assoc,cbrt.(b),assoc_options)
    packagedparams = CPCAParam(Mw, Tc, segment, a, b, c1, kb1, kb2, epsilon_assoc, bondvol)
    #init cubic model
    

    references = ["10.1021/acs.iecr.3c02774"]
    model = CPCA(_components, :KG, init_cubicmodel, packagedparams, sites, init_idealmodel, assoc_options, references)
    return model
end

function recombine_impl!(model::CPCAModel)
    assoc_options = model.assoc_options
    a = model.params.a
    b = model.params.b
    model.cubicmodel.alpha.params.c1.values .= model.params.c1.values
    model.cubicmodel.mixing.params.segment.values .= model.params.segment.values
    model.cubicmodel.mixing.params.kb1.values .= model.params.kb1.values
    model.cubicmodel.mixing.params.kb2.values .= model.params.kb2.values
    a  = epsilon_LorentzBerthelot!(a)
    b  = sigma_LorentzBerthelot!(b)
    epsilon_assoc = model.params.epsilon_assoc
    bondvol = model.params.bondvol
    bondvol,epsilon_assoc = assoc_mix(bondvol,epsilon_assoc,cbrt.(b),assoc_options) #combining rules for association

    model.params.epsilon_assoc.values.values[:] = epsilon_assoc.values.values
    model.params.bondvol.values.values[:] = bondvol.values.values
    return model
end

function data(model::CPCAModel, V, T, z)
    nabc = data(model.cubicmodel,V,T,z)
    n,ā,b̄,c̄ = nabc
    m = model.params.segment.values
    m̄n = dot(z,m)
    β = m̄n*b̄/V
    m̄ = m̄n/n
    return nabc,β,m̄
end

function x0_volume_liquid(model::CPCAModel,T,z)
    nabc = data(model.cubicmodel,0.0,T,z)
    n,ā,b̄,c̄ = nabc
    return (1.25b̄ + c̄)*n
end

function a_res(model::CPCAModel, V, T, z, _data = @f(data))
    nabc,β,m̄ = _data
    n,ā,b̄,c̄ = nabc
    return m̄*a_res(model.cubicmodel,V,T,z,nabc) + a_chain(model,V+c̄*n,T,z,_data) + a_assoc(model,V+c̄*n,T,z,_data)
end

function a_chain(model::CPCAModel,V,T,z,_data = @f(data))
    nabc,β,m̄ = _data
    n,ā,b̄,c̄ = nabc
    res = zero(V + first(z))
    m = model.params.segment.values
    ln_ghs = -log1p(-0.475*β) #log(1/(1-0.475β))
    for i in @comps
        res = z[i]*(m[i] - 1)
    end
    return -res*ln_ghs/n
end

function Δ(model::CPCAModel, V, T, z, i, j, a, b, _data = @f(data))
    nabc,β,m̄ = _data
    ϵ_associjab = model.params.epsilon_assoc.values[i,j][a,b]
    βijab = model.params.bondvol.values[i,j][a,b]
    b = model.params.b.values
    rdf = model.radial_dist
    bi0,bj0 = b[i,i],b[j,j]
    _1_m_l = 2*b[i,j]/(bi0 + bj0) #bij = (bi + bj)*0.5*(1 - lj)
    m = model.params.segment.values
    kb1 = model.params.kb1.values
    kb2 = model.params.kb2.values
    Tc = model.cubicmodel.params.Tc.values
    βi = kb1[i]*exp(-kb2[i]*T/Tc[i])
    βj = kb1[j]*exp(-kb2[j]*T/Tc[j])
    bi = bi0 * βi
    bj = bj0 * βj
    #@show βi,βj
    bij = (bi + bj)*_1_m_l*0.5
    g = 1/(1-0.475β)
    return g*expm1(ϵ_associjab/T)*βijab*bij/N_A
end

export CPCA,CPCAMixing