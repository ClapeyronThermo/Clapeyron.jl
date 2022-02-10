
#original
#f(z) = eos(model,V,T,z)
#H(z) = ForwardDiff.hessian(f,z)/(R̄*T)
#L(z) = det(H(z))
#dL(z) = ForwardDiff.gradient(L,z)
#M(z) = [H(z)[1:end-1,:];transpose(dL(z))]
"""
    LdetM(model,V,T,z)

Used for UCEP and UCST_mix. 
with `a(x)` the reduced `(A/RT)` helmholtz energy dependent on composition `xᵢ` for `i` ∈ `1:n`, returns `L` and `det(M)`, where `L` and `M` are defined as:
```
L := det(ℍ(a)) (ℍ = hessian)
M := ℍ(a) for rows ∈ 1:n-1
  := ∇L for row n
```
"""
function LdetM(model,V,T,z)
    f(x) = eos(model,V,T,x)/(R̄*T)
    H(x) = ForwardDiff.hessian(f,x) #∂A/∂zᵢ∂zⱼ == ∂A/∂zⱼ∂zᵢ
    L(x) = det(Symmetric(H(x)))
    dL(x) = ForwardDiff.gradient(L,x)
    HH = H(z)   
    LL = det(HH)    
    Mᵢ = @view(HH[end,:])
    Mᵢ .=  dL(z)
    MM = HH
    #M(x) = [HH[1:end-1,:];transpose(dL(x))]
    return LL , det(MM)
end


include("rachford_rice.jl")
include("bubble_point.jl")
include("dew_point.jl")
include("azeotrope_point.jl")
include("LLE_point.jl")
include("VLLE.jl")
include("crit_mix.jl")
include("UCEP.jl")
include("UCST_mix.jl")
include("tp_flash.jl")
export bubble_pressure,    dew_pressure,    LLE_pressure,    azeotrope_pressure, VLLE_pressure
export bubble_temperature, dew_temperature, LLE_temperature, azeotrope_temperature, VLLE_temperature
export crit_mix, UCEP_mix, UCST_mix
