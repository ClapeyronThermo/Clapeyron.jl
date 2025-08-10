abstract type LiquidVolumeModel <: EoSModel end


#differentials
function ∂𝕧∂T(model,p,T,z::AbstractVector)
    v(∂T) = simple_volume(model,p,∂T,z)
    return Solvers.derivative(v,T)
end

function ∂𝕧∂p(model,p,T,z::AbstractVector)
    v(∂p) = simple_volume(model,∂p,T,z)
    return Solvers.derivative(v,p)
end

function ∂𝕧(model,p,T,z)
    f(∂p,∂T) = simple_volume(model,∂p,∂T,z)
    _f,_df = Solvers.fgradf2(f,p,T)
    return _df,_f
end

function ∂𝕧_vec(model,p,T,z::AbstractVector)
    _df,_f = ∂𝕧(model,p,T,z)
    return SVector(_f,_df[1],_df[2])
end

function 𝕧∂𝕧dp(model,p,T,z::AbstractVector)
    f(x) = simple_volume(model,x,T,z)
    V,∂V∂p = Solvers.f∂f(f,p)
    return SVector(V,∂V∂p)
end

function 𝕧∂𝕧dT(model,p,T,z::AbstractVector)
    f(x) = simple_volume(model,p,x,z)
    V,∂V∂T = Solvers.f∂f(f,T)
    return SVector(V,∂V∂T)
end

function ∂2𝕧(model,p,T,z)
    f(_p,_T) = simple_volume(model,_p,_T,z)
    _f,_∂f,_∂2f = Solvers.∂2(f,p,T)
    return (_∂2f,_∂f,_f)
end
 
function 𝕧_hess(model,p,T,z)
    f(w) = simple_volume(model,first(w),last(w),z)
    p,T = promote(p,T)
    pT_vec = SVector(p,T)
    return Solvers.hessian(f,pT_vec)
end

function ∂²𝕧∂T²(model,p,T,z)
    V(x) = simple_volume(model,p,x,z)
    ∂V∂T(x) = Solvers.derivative(V,x)
    ∂²V∂T²(x) = Solvers.derivative(∂V∂T,x)
    return ∂²V∂T²(T)
end

function PT_property(model::LiquidVolumeModel,p,T,z,phase,threaded,vol0,f::F,USEP::Val{UseP}) where {F,UseP}
    return PT_property_volume_based(model,p,T,z,f)
end

function PT_property_volume_based(model,p,T,z,f::typeof(VT_isothermal_compressibility),USEP::Val{UseP}) where {UseP}
    v,dvdp = 𝕧∂𝕧dp(model,p,T,z)
    return -dvdp/v
end

function PT_property_volume_based(model,p,T,z,f::typeof(VT_isobaric_expansivity),USEP::Val{UseP}) where {UseP}
    v,dvdT = 𝕧∂𝕧dT(model,p,T,z)
    return -dvdT/v
end

#we suppose we have an expression for a_res here
function PT_property_volume_based(model,p,T,z,f,USEP::Val{UseP}) where {UseP}
    v = simple_volume(model,v,T,z)
    if UseP
        return f(model,v,T,z,p)
    else
        return f(model,v,T,z)
    end
end

include("NaNLiquid/NaNLiquid.jl")
include("Rackett/RackettLiquid.jl")
include("COSTALD/costald.jl")
include("DIPPR105Liquid/DIPPR105Liquid.jl")
include("PolExpLiquid/PolExpLiquid.jl")
include("GrenkeElliott/GrenkeElliottWater.jl")