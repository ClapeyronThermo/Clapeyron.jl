mw(model::EoSModel) = paramvals(model.params.Mw)

function group_molecular_weight(groups::GroupParam,mw,z = @SVector [1.])
    n = groups.n_flattenedgroups
    res = zero(first(z))
    Σz = sum(z)
    @inbounds for i in 1:length(groups.components)
        ni = n[i]
        gi = groups.i_groups[i]
        mwi = zero(res)
        for idx in 1:length(gi)
            j = gi[idx]
            mwi += mw[j]*ni[j]
        end
        res +=z[i]*mwi
    end
    return 0.001*res/Σz
end


comp_molecular_weight(mw,z = @SVector [1.]) = 0.001*mapreduce(*,+,mw,z)

include("initial_guess.jl")
include("differentials.jl")
include("VT.jl")
include("property_solvers/property_solvers.jl")
include("pT.jl")
include("unitful_base.jl")
include("unitful_methods.jl")

