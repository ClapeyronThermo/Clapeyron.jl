# The unnormalized equations and binary pivot adapt Clapeyron's MIT-licensed
# fixed-composition criticality code, covered by the notice in ../../LICENSE.md.
module NormalizedCriticality
using Clapeyron,LinearAlgebra
const FD=Clapeyron.ForwardDiff
export binary_constraints
"""Binary-only scaling by the separated positive transverse eigenvalue.
Same roots on the simple-zero-mode manifold. Does not certify stability.
Returns raw and normalized constraints explicitly; their tolerances are not interchangeable.
"""
function binary_constraints(model,V,T,z)
 length(z)==2 || throw(ArgumentError("binary only"))
 length(model)==2 || throw(ArgumentError("binary model required"))
 isfinite(V) && V>0 && isfinite(T) && T>0 || throw(ArgumentError("finite positive volume and temperature required"))
 all(x->isfinite(x) && x>0,z) || throw(ArgumentError("finite strictly positive amounts required"))
 fr(x)=sum(x)*Clapeyron.a_res(model,V,T,x)
 function B(x)
  sx=sqrt.(x);BB=FD.hessian(fr,x).*(sx*transpose(sx))
  BB[1,1]+=one(eltype(BB));BB[2,2]+=one(eltype(BB));BB
 end
 L(x)=det(Symmetric(B(x)))
 H=B(z);ell=(H[1,1]+H[2,2]+sqrt((H[1,1]-H[2,2])^2+4H[1,2]^2))/2
 ell>1e-5 || throw(DomainError(ell,"positive transverse eigenvalue required"))
 ll=det(Symmetric(H));grad=FD.gradient(L,z).*sqrt.(z)
 pivot=abs(H[2,2])>=abs(H[1,1]) ? 1 : 2
 H[pivot,:].=grad;mm=det(H)
 (;normalized=[ll/ell,mm/ell^2],raw=[ll,mm],lambda_transverse=ell,pivot)
end
end
