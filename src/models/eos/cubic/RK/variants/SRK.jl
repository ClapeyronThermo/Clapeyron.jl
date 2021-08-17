function cubic_ab(model::RK{<:Any,SoaveAlpha},T,z=SA[1.0],n=sum(z))
    a = model.params.a.values
    b = model.params.b.values
    αmodel = model.alpha
    ω = αmodel.params.acentricfactor.values
    Tc = model.params.Tc.values
    _1 = one(T+n)
    invn = one(n)/n
    αx = @. (_1+(0.480+1.547*ω-0.176*ω^2)*(1-√(T/Tc))) * z * invn 
    āᾱ =dot(αx, Symmetric(a), αx)
    b̄ = dot(z, Symmetric(b), z) * invn * invn
    return āᾱ ,b̄
end
#just a function, no struct
function SRK(components::Vector{String}; idealmodel=BasicIdeal,
    activity = nothing,
    userlocations=String[], 
    ideal_userlocations=String[],
    alpha_userlocations = String[],
    activity_userlocations = String[],
     verbose=false)

     return RK(components;
     idealmodel = idealmodel,
     alpha = SoaveAlpha,
     activity=activity,
     ideal_userlocations = ideal_userlocations,
     alpha_userlocations = alpha_userlocations,
     activity_userlocations = activity_userlocations,
     verbose = verbose)
end
export SRK