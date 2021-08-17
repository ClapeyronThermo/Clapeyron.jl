
#TODO
# - standarize chemical group notation to make it the same as SAFTgammaMie
# - add a database of group Mw

struct JobackIdealParam <: EoSParam
    Mw::SingleParam{Float64}
    N_a::SingleParam{Int}
    T_c::SingleParam{Float64}
    P_c::SingleParam{Float64}
    V_c::SingleParam{Float64}
    T_b::SingleParam{Float64}
    T_m::SingleParam{Float64}
    H_form::SingleParam{Float64}
    G_form::SingleParam{Float64}
    a::SingleParam{Float64}
    b::SingleParam{Float64}
    c::SingleParam{Float64}
    d::SingleParam{Float64}
    H_fusion::SingleParam{Float64}
    H_vap::SingleParam{Float64}
    eta_a::SingleParam{Float64}
    eta_b::SingleParam{Float64}
end


abstract type JobackIdealModel <: IdealModel end
@newmodelgc JobackIdeal JobackIdealModel JobackIdealParam

export JobackIdeal
function JobackIdeal(components;userlocations=String[], verbose=false)
    groups = GroupParam(components,["ideal/JobackIdeal_Groups.csv"], verbose=verbose)
    params = getparams(groups, ["ideal/JobackIdeal.csv","properties/molarmass_groups.csv"]; userlocations=userlocations, verbose=verbose)
    Mw = params["Mw"]::SingleParam{Float64}
    N_a = params["N_a"]::SingleParam{Int}
    T_c = params["T_c"]::SingleParam{Float64}
    P_c = params["P_c"]::SingleParam{Float64}
    V_c = params["V_c"]::SingleParam{Float64}
    T_b = params["T_b"]::SingleParam{Float64}
    T_m = params["T_m"]::SingleParam{Float64}
    H_form = params["H_form"]::SingleParam{Float64}
    G_form = params["G_form"]::SingleParam{Float64}
    a = params["a"]::SingleParam{Float64}
    b = params["b"]::SingleParam{Float64}
    c = params["c"]::SingleParam{Float64}
    d = params["d"]::SingleParam{Float64}
    H_fusion = params["H_fusion"]::SingleParam{Float64}
    H_vap = params["H_vap"]::SingleParam{Float64}
    eta_a = params["eta_a"]::SingleParam{Float64}
    eta_b = params["eta_b"]::SingleParam{Float64}
    packagedparams = JobackIdealParam(
    Mw,
    N_a,
    T_c,
    P_c,
    V_c,
    T_b,
    T_m,
    H_form,
    G_form,
    a,
    b,
    c,
    d,
    H_fusion,
    H_vap,
    eta_a,
    eta_b)
    references = ["10.1080/00986448708960487"]
    sites = SiteParam(groups.components)
    model = JobackIdeal(packagedparams, groups, sites, BasicIdeal, references=references, verbose=verbose)
    return model
end


function C_p(model::JobackIdeal,T,z=SA[1.0])
    a = model.params.a.values
    b = model.params.b.values
    c = model.params.c.values
    d = model.params.d.values
    n = model.groups.n_flattenedgroups
    res = zero(T+first(z))
    Σz = sum(z)
    @inbounds for i in @comps
        ni = n[i]
        c0 = ∑(a[j]*ni[j] for j in @groups(i)) - 37.93
        c1 = ∑(b[j]*ni[j] for j in @groups(i)) + 0.210
        c2 = ∑(c[j]*ni[j] for j in @groups(i)) - 3.91e-4
        c3 = ∑(d[j]*ni[j] for j in @groups(i)) + 2.06e-7
        pol = (c0,c1,c2,c3)
        res +=z[i]*evalpoly(T,pol)
    end
    return res/Σz
end



function T_b(model::JobackIdeal)
    n = model.groups.n_flattenedgroups
    T_b = model.params.T_b.values
    res = 0.0
    @inbounds begin
        ni = n[1]
        res += ∑(T_b[j]*ni[j] for j in @groups(1))
    end
    res = res + 198.2
end

function T_c(model::JobackIdeal,Tb=T_b(model))
    n = model.groups.n_flattenedgroups
    T_c = model.params.T_c.values
    res = 0.0
    @inbounds begin
        ni = n[1]
        ΣT_ci = ∑(T_c[j]*ni[j] for j in @groups(1))
        res += Tb * (0.584+0.965*ΣT_ci - ΣT_ci^2)^-1
    end
    return res
end

#this style of writing is uglier but faster,
#ideally they should be of the same speed
function V_c(model::JobackIdeal)
    n = model.groups.n_flattenedgroups
    V_c = model.params.V_c.values
    res = 0.0
    groups = @groups(1)
    ni = n[1]
    @inbounds begin
        ΣV_ci = zero(res)
        for idx in 1:length(groups)
            j = groups[idx]
            ΣV_ci +=V_c[j]*ni[j]
        end
        res += 17.5 + ΣV_ci
    end
    #result in mL/mol, converting to m3/mol
    return res*1e-6
end

function P_c(model::JobackIdeal)
    n = model.groups.n_flattenedgroups
    P_c = model.params.P_c.values
    N_a = model.params.N_a.values
    groups = @groups(1)
    ni = n[1]
    @inbounds begin
        ΣP_ci  = 0.0
        ΣN_a = 0
        for idx in 1:length(groups)
            j = groups[idx]
            ΣP_ci += P_c[j]*ni[j]
            ΣN_a += N_a[j]*ni[j]
        end
        res = (0.113 + 0.0032*ΣN_a - ΣP_ci)^-2
    end
    #result in bar, converting to Pa
    return res*100000.0
end

function crit_pure(model::JobackIdeal)
    return (T_c(model),P_c(model),V_c(model))
end

function a_ideal(model::JobackIdealModel, V, T, z)
    a = model.params.a.values
    b = model.params.b.values
    c = model.params.c.values
    d = model.params.d.values
    n = model.groups.n_flattenedgroups
    res = zero(V+T+first(z))
    Σz = sum(z)
    x = z/Σz
    @inbounds for i in @comps
        ni = n[i]
        c0 = ∑(a[j]*ni[j] for j in @groups(i)) - 37.93
        c1 = ∑(b[j]*ni[j] for j in @groups(i)) + 0.210
        c2 = ∑(c[j]*ni[j] for j in @groups(i)) - 3.91e-4
        c3 = ∑(d[j]*ni[j] for j in @groups(i)) + 2.06e-7
        polycoeff = [c0,c1,c2,c3]
        res +=x[i]*(log(z[i]/V) + 1/(R̄*T)*(sum(polycoeff[k]/k*(T^k-298^k) for k in 1:4)) -
        1/R̄*((polycoeff[1]-R̄)*log(T/298)+sum(polycoeff[k]/(k-1)*(T^(k-1)-298^(k-1)) for k in 2:4)))
    end
    return res
end

export JobackIdeal
#acetone = [("acetone",["−CH3"=>1,">C=O (non-ring)"=>1])]