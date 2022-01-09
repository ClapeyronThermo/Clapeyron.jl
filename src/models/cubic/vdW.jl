struct vdWParam <: EoSParam
    a::PairParam{Float64}
    b::PairParam{Float64}
    Tc::SingleParam{Float64}
    Pc::SingleParam{Float64}
    Mw::SingleParam{Float64}
end

abstract type vdWModel <: ABCubicModel end

struct vdW{T <: IdealModel,α,c,M} <: vdWModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    alpha::α
    mixing::M
    translation::c
    params::vdWParam
    idealmodel::T
    absolutetolerance::Float64
    references::Array{String,1}
end

@registermodel vdW
export vdW

function vdW(components::Vector{String}; idealmodel=BasicIdeal,
    alpha = NoAlpha,
    mixing = vdW1fRule,
    activity=nothing,
    translation=NoTranslation,
    userlocations=String[], 
    ideal_userlocations=String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
     verbose=false)
    params = getparams(components, ["properties/critical.csv", "properties/molarmass.csv","SAFT/PCSAFT/PCSAFT_unlike.csv"]; userlocations=userlocations, verbose=verbose)
    k  = params["k"]
    pc = params["pc"]
    Mw = params["Mw"]
    Tc = params["Tc"]
    init_mixing = init_model(mixing,components,activity,mixing_userlocations,activity_userlocations,verbose)
    a,b = ab_premixing(vdW,init_mixing,Tc,pc,k)
    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    init_translation = init_model(translation,components,translation_userlocations,verbose)
    init_alpha = init_model(alpha,components,alpha_userlocations,verbose)
    icomponents = 1:length(components)
    packagedparams = vdWParam(a,b,Tc,pc,Mw)
    references = String[]
    model = vdW(components,icomponents,init_alpha,init_mixing,init_translation,packagedparams,init_idealmodel,1e-12,references)
    return model
end

function ab_consts(::Type{<:vdWModel})
    Ωa =  27/64
    Ωb =  1/8
    return Ωa,Ωb
end


function cubic_abp(model::vdWModel, V, T, z)
    n = ∑(z)
    a,b,c = cubic_ab(model,V,T,z,n)
    v = V/n+c
    p = R̄*T/(v-b) - a/(v^2)
    return a,b,p
end

function cubic_poly(model::vdWModel,p,T,z)
    n = sum(z)
    a,b,c = cubic_ab(model,p,T,z,n)
    _1 = one(a+b)
    RT⁻¹ = 1/(R̄*T)
    A = a*p*RT⁻¹*RT⁻¹
    B = b*p*RT⁻¹
    return (-A*B, A, -B-_1, _1),c
end



function a_res(model::vdWModel, V, T, z)
    n=sum(z)
    ā,b̄,c̄ = cubic_ab(model,V,T,z,n)
    RT⁻¹ = 1/(R̄*T)
    ρt = (V/n+c̄)^(-1) # translated density
    ρ  = n/V
    return -log(1+(c̄-b̄)*ρ) - ā*ρt*RT⁻¹
    #
    #return -log(V-n*b̄) - ā*n/(R̄*T*V) + log(V)
end   

cubic_zc(::vdWModel) = 3/8
