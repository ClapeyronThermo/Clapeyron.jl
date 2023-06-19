struct QuadraticDepartureParam <: EoSParam
    k0::PairParam{Float64}
    k1::PairParam{Float64}
end

@newmodelsimple QuadraticDeparture MultiFluidDepartureModel QuadraticDepartureParam

function QuadraticDeparture(components::AbstractVector, userlocations=String[], verbose::Bool=false)
    params = getparams(components,["Empiric/departure/Quadratic_departure_unlike.csv"];userlocations = userlocations,verbose = verbose)
    k0 = get(params,"k0",nothing)
    k1 = get(params,"k1",nothing)
    k0 === nothing && (k0 = PairParam("k0",components))
    k1 === nothing && (k1 = PairParam("k0",components))
    pkgparams = QuadraticDepartureParam(k0,k1)
    references = ["10.1021/acs.iecr.1c01186","10.1016/j.fluid.2018.04.015"]
    return QuadraticDeparture(pkgparams,references)
end

function multiparameter_a_res(model,V,T,z,departure::QuadraticDeparture,δ,τ,∑z = sum(z))
    lnδ = log(δ)
    lnτ = log(τ)
    _0 = zero(lnδ+lnτ)
    n = length(model)
    aᵣₖ = fill(_0,length(model))
    m = model.pures
    Rinv = 1/Rgas(model)
    for i in 1:n
        mᵢ = m[i]
        aᵣᵢ[i] = reduced_a_res(mᵢ,δ,τ,lnδ,lnτ)*Rinv*Rgas(mᵢ)
    end
    k₀ = departure.params.k0
    k₁ = departure.params.k1
    aᵣ = _0
    for i in 1:n
        aᵢ = aᵣₖ[i]
        zᵢ = z[i]
        aᵣ += aᵢ*zᵢ*zᵢ
        for j in 1:(i-1)
            aᵢⱼ += 2*zᵢ*z[j]*0.5*(aᵢ + aᵣₖ[j])*(1 - k₀[i,j] - k₁[i,j]*T)
        end
    end
    return aᵣ/(∑z*∑z)
end