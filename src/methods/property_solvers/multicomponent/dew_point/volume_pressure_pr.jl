# Private per-call opt-in; no new exported method and no change to the default.
struct _PRVolumeBubblePressure{M<:ChemPotBubblePressure} <: BubblePointMethod
    base::M
end

struct _PRVolumeDewPressure{M<:ChemPotDewPressure} <: DewPointMethod
    base::M
end

struct _PRVolumeJacobian end
const _PRVolumePressure = Union{_PRVolumeBubblePressure,_PRVolumeDewPressure}

_pr_volume_no_exclusions(m::ChemPotBubblePressure) = m.nonvolatiles === nothing
_pr_volume_no_exclusions(m::ChemPotDewPressure) = m.noncondensables === nothing

# Explicit allowlist: constant shape, fixed-T independent component alpha,
# vdW1f quadratic pairs and a component-additive translation. Do not admit
# arbitrary DeltaCubicModel, alpha, mixing or translation subclasses.
function _constshape_volume_model(model)
    # This file is loaded before the concrete cubic models.
    return model isa Union{
        PR{BasicIdeal,PRAlpha,<:Union{NoTranslation,ConstantTranslation,PenelouxTranslation},vdW1fRule},
        RK{BasicIdeal,SoaveAlpha,NoTranslation,vdW1fRule},
        vdW{BasicIdeal,NoAlpha,NoTranslation,vdW1fRule}
    }
end

function _pr_volume_supported(m, model, T, x)
    _constshape_volume_model(model) || return false
    has_dual(model) && return false
    m isa Union{ChemPotBubblePressure{Float64},ChemPotDewPressure{Float64}} || return false
    return T isa Float64 && isfinite(T) && T > 0 && eltype(x) === Float64 &&
        2 <= length(x) == length(model) <= 5 && all(t -> isfinite(t) && t > 0, x) &&
        !m.verbose && _pr_volume_no_exclusions(m) && issymmetric(model.params.b.values) &&
        m.max_iters == 1000 && m.atol == 1e-8 && m.rtol == 1e-12 && iszero(m.f_limit)
end

# Gate before public normalization, model/phase reduction and primal conversion.
_select_volume_pressure_method(m, model, T, x) = m
_select_volume_pressure_method(m::_PRVolumePressure, model, T, x) =
    _pr_volume_supported(m.base, model, T, x) ? m : m.base

get_verbosity(m::_PRVolumePressure) = get_verbosity(m.base)
index_reduction(m::_PRVolumeBubblePressure, idx) =
    all(idx) ? m : index_reduction(m.base, idx)
index_reduction(m::_PRVolumeDewPressure, idx) =
    all(idx) ? m : index_reduction(m.base, idx)
Solvers.primalval(m::_PRVolumeBubblePressure) = _PRVolumeBubblePressure(primalval(m.base))
Solvers.primalval(m::_PRVolumeDewPressure) = _PRVolumeDewPressure(primalval(m.base))

bubble_pressure_impl(model::EoSModel, T, x, m::_PRVolumeBubblePressure) =
    bubble_pressure_impl(model, T, x, m.base, _PRVolumeJacobian())
dew_pressure_impl(model::EoSModel, T, y, m::_PRVolumeDewPressure) =
    dew_pressure_impl(model, T, y, m.base, _PRVolumeJacobian())

_volume_pressure_solve(::Nothing, f!, z, lb, ub, m, model, T, x, mask, pivot, side) =
    _bubbledew_chempot_solve(f!, z, lb, ub, m, Val(min(length(z), 8)))

function _volume_pressure_solve(::_PRVolumeJacobian, f!, z, lb, ub, m, model, T, x,
                                mask, pivot, ::Val{Side}) where Side
    # An in-place attempt copies the original seed before the bounded Newton call.
    n = length(x)
    if 2 <= n <= 5
        candidate = _pr_volume_try(f!, model, T, x, z, lb, ub, m, mask, pivot,
                                   Side === :bubble, Val(n))
        candidate !== nothing && return candidate, true
    end
    return _volume_pressure_solve(nothing, f!, z, lb, ub, m, model, T, x, mask, pivot, Val(Side))
end

# For admitted fixed-shape cubics, use actual A and B pair coefficients
# (including kij/lij) in untranslated volume. For c=n/V′,
# Φ(c)=A_res^untranslated(1,T,c)/(RT), g=∇Φ, H=∇²Φ:
#     p/(RT) = Σc + c⋅g - Φ,      ∇c[p/(RT)] = 1 + H*c.
# Fixed composition needs only H*c; the free-phase columns need H explicitly.
struct _PRVolumeDomainError <: Exception end
# Shape and pair coefficients are independent of the phase densities and
# compositions on the allowlist. Rebuild per fixed-T call.
function _pr_volume_coeff(model,T,xx,::Val{N}) where N
    α=α_function(model,1.0,T,xx,model.alpha); a=model.params.a.values
    b=model.params.b.values; RT=Rgas(model)*T
    A=SMatrix{N,N}(ntuple(Val(N*N)) do k
       i=rem(k-1,N)+1;j=div(k-1,N)+1
       a[max(i,j),min(i,j)]*sqrt(α[i]*α[j])/RT
    end)
    B=SMatrix{N,N}(ntuple(Val(N*N)) do k
       i=rem(k-1,N)+1;j=div(k-1,N)+1
       b[max(i,j),min(i,j)]
    end)
    u,w=cubic_uwT(model,T,xx)
    (;A,B,u,w)
end
function _pr_volume_jet(coeff,c,::Val{wantH}) where wantH
    (; A,B,u,w)=coeff; s=sum(c); b=B*c; a=A*c
    χ=dot(c,b)/s
    isfinite(s) && s > 0 && 0 <= χ < 1 || throw(_PRVolumeDomainError())
    q=dot(c,a); k=inv(1-χ); L=-log1p(-χ)
    D=1+χ*(u+w*χ)
    isfinite(D) && D>0 || throw(_PRVolumeDomainError())
    if wantH === true
        K,K1,K2=_cubic_attraction_K(χ,u,w,Val(2))
    else
        K,K1=_cubic_attraction_K(χ,u,w,Val(1))
    end
    du=(2*b .- χ)/s; dq=2*a; g1=s*k-q*K1
    g=L .+ g1*du .- K*dq
    p=s*k-q/D
    wantH === :gp && return g,p
    h=(k-1) .+ (s*k*k+q*(u+2w*χ)/(D*D))*du .- dq/D
    r=s+dot(c,h)
    if wantH
        H=_pr_volume_H(B,A,du,dq,g1,s,k,q,K,K1,K2)
        return g,h,p,r,H
    else
        return g,h,p,r
    end
end
function _pr_volume_H(B::SMatrix{N,N},A,du,dq,g1,s,k,q,K,K1,K2) where N
    SMatrix{N,N}(ntuple(Val(N*N)) do n
        i=rem(n-1,N)+1;j=div(n-1,N)+1
        k*(du[i]+du[j])+(s*k*k-q*K2)*du[i]*du[j]+g1*(2B[i,j]-du[i]-du[j])/s -
          K1*(dq[i]*du[j]+du[i]*dq[j])-2K*A[i,j]
    end)
end
function _pr_volume_lower(coeff,w)
    Bw=coeff.B*w; n=sum(w); l=dot(Bw,w)/n^2
    return l, (2/n^2)*(Bw .- n*l)
end
function _pr_volume_objective(f!,model,T,x,pivot,nn::Val{N},bubble) where N
    xx=SVector{N}(x); logx=log.(xx); is=inv(p_scale(model,xx))
    # Stock η = log(v_phys - b̄ + c̄) = log(v′ - b̄), with v′=v_phys+c̄.
    # Construct v′ inside the kernel; the stock caller decodes the same η
    # into v_phys and uses physical volumes for the acceptance and return.
    lox=cubic_lb_volume(model,T,xx)/sum(xx)
    RT=Rgas(model)*T; rtps=RT*is; coeff=_pr_volume_coeff(model,T,xx,nn)
    ix,iw=bubble ? (1,2) : (2,1)
    function fj!(F,J,z)
        w=SVector{N}(FractionVector(@view(z[3:end]),pivot))
        all(>(0),w) || throw(_PRVolumeDomainError())
        lo,dl=_pr_volume_lower(coeff,w)
        ex,ew=exp(z[ix]),exp(z[iw]); vx,vw=ex+lox,ew+lo
        isfinite(vx) && isfinite(vw) && vx > max(0, lox) && vw > max(0, lo) ||
            throw(_PRVolumeDomainError())
        logv=log(vw)-log(vx)
        if J === nothing
            g0,p0=_pr_volume_jet(coeff,xx/vx,Val(:gp))
            g1,p1=_pr_volume_jet(coeff,w/vw,Val(:gp))
            @inbounds for i in 1:N
                F[i]=g0[i]-g1[i]+logv+logx[i]-log(w[i])
            end
            F[N+1]=(p0-p1)*rtps
            all(isfinite,F) || throw(_PRVolumeDomainError())
            return F
        end
        g0,h0,p0,r0=_pr_volume_jet(coeff,xx/vx,Val(false))
        g1,h1,p1,r1,H=_pr_volume_jet(coeff,w/vw,Val(true))
        rx=ex/vx; rw=ew/vw
        @inbounds for i in 1:N
            F[i]=g0[i]-g1[i]+logv+logx[i]-log(w[i])
            J[i,ix]=-(1+h0[i])*rx
            J[i,iw]=(1+h1[i])*rw
        end
        F[N+1]=(p0-p1)*rtps
        J[N+1,ix]=-r0*rx*rtps
        J[N+1,iw]=r1*rw*rtps
        # v′ = exp(η) + b̄(w). Moving w_j also moves the pivot with
        # derivative -1, so b̄ changes by ∂b̄/∂w_j - ∂b̄/∂w_pivot.
        k=2
        @inbounds for j in 1:N
            j==pivot && continue
            k+=1; dv=dl[j]-dl[pivot]; rel=dv/vw
            for i in 1:N
                J[i,k]=(H[i,pivot]-H[i,j])/vw - ((i==j)-(i==pivot))/w[i] + (1+h1[i])*rel
            end
            J[N+1,k]=rtps*((h1[pivot]-h1[j])/vw+r1*rel)
        end
        all(isfinite, F) && all(isfinite, J) || throw(_PRVolumeDomainError())
        return F,J
    end
    ff=similar(x,N+1)
    j!(J,z)=(fj!(ff,J,z); J)
    # Solve chemo-potential equality in untranslated coordinates. This omits
    # -c_i*(p_fixed-p_free)/(RT) from each physical chemical residual; it is
    # zero at pressure equality. Use the original f! for acceptance/fallback.
    gp!(F,z)=fj!(F,nothing,z)
    return NLSolvers.VectorObjective(gp!,j!,fj!,nothing)
end
function _pr_volume_solve(f!,model,T,x,z,lb,ub,m,pivot,bubble,n)
    obj=_pr_volume_objective(f!,model,T,x,pivot,n,bubble); seed=copy(z)
    # Regular admitted cases converge in a few steps. A difficult trial
    # must not spend the full stock iteration budget before a stock retry.
    opts = NLSolvers.NEqOptions(f_limit=m.f_limit, f_abstol=m.atol,
                                 f_reltol=m.rtol, maxiter=min(m.max_iters, 12))
    return Solvers.nlsolve(NLSolvers.NEqProblem(obj; inplace=true), seed,
                           LineSearch(Newton2(seed), BoundedLineSearch(lb,ub)), opts)
end
function _pr_volume_try(f!, model, T, x, z, lb, ub, m, mask, pivot, bubble, n::Val)
    _pr_volume_supported(m, model, T, x) && all(mask) && z isa Vector{Float64} &&
        all(isfinite, z) || return nothing
    try
        r=_pr_volume_solve(f!,model,T,x,z,lb,ub,m,pivot,bubble,n)
        __check_convergence(r) || return nothing
        s=Solvers.x_sol(r)
        w=FractionVector(@view(s[3:end]),pivot)
        all(t->isfinite(t)&&t>0,w) || return nothing
        ix,iw=bubble ? (1,2) : (2,1)
        vx=v_from_η(model,s[ix],T,x); vw=v_from_η(model,s[iw],T,w)
        isfinite(vx) && isfinite(vw) && min(vx,vw)>0 || return nothing
        all(t->abs(t)<m.atol,f!(similar(z),s)) || return nothing
        p=pressure(model,vx,T,x); pw=pressure(model,vw,T,w)
        isfinite(pw) && pw>0 && abs(p-pw)<=1e-6*p || return nothing
        (bubble ? vx<vw : vx>vw) && isfinite(p) && p>0 && bubbledew_check(model,p,T,vw,vx,w,x) || return nothing
        return s
    catch err
        (err isa _PRVolumeDomainError || err isa DomainError || err isa LinearAlgebra.SingularException) || rethrow()
        return nothing
    end
end
