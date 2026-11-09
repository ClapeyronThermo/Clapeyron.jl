function material_balance_rr_converged(w,z,β::Number,n = sum(z),ztol = sqrt(eps(eltype(β))))
    βx = SVector(1 - β,β)
    return material_balance_rr_converged(w,z,βx,n,ztol)
end

function material_balance_rr_converged(w,z,β::AbstractVector,n = sum(z),ztol = sqrt(eps(eltype(β))))
    rk = zero(Base.promote_eltype(w[1],z,β))
    np = length(β)
    #material balance test from https://github.com/WhitsonAS/Rachford-Rice-Contest
    for i in 1:length(z)
        rk1 = zero(rk)
        rk2 = zero(rk)
        zi = z[i]/n
        for j in 1:np
            Vi = β[j]*w[j][i]
            rk1 +=  Vi
            rk2 += abs(Vi)
        end
        if !iszero(zi)
            rk = max(rk,abs(rk1 - zi)/abs(rk2 + zi))
        end
    end
    return rk <= ztol
end
