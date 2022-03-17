function isstable(model,V,T,z)
    stable = true
    β = VT_isothermal_compressibility(model,V,T,z)
    if β<0
        @warn "StabilityWarning: Phase is mechanically unstable"
        stable = false
    end
    A(x) = eos(model,V,T,x)
    Hf = ForwardDiff.hessian(A,z)
    (Λ,U)=eigen(Hf)
    λ = minimum(Λ)
    if λ <0
        @warn "StabilityWarning: Phase is diffusively unstable"
        stable = false
    end
    return stable
end