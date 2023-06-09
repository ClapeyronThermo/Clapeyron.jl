abstract type XiangDeitersModel <: EoSModel end

struct XiangDeitersParam <: EoSParam
    Mw::SingleParam{Float64}
    Tc::SingleParam{Float64}
    Pc::SingleParam{Float64}
    Vc::SingleParam{Float64}
    acentricfactor::SingleParam{Float64}
end

struct XiangDeiters{T <: IdealModel} <: XiangDeitersModel
    components::Array{String,1}
    params::XiangDeitersParam
    idealmodel::T
    consts::EmpiricSingleFluidResidualParam
    references::Array{String,1}
end

@registermodel XiangDeiters

function XiangDeitersConsts(ω,θ)
    d = [1, 1, 1, 2, 3, 7, 1, 1, 2, 5, 1, 1, 4, 2]
    t = [0.25, 1.25, 1.5, 1.375, 0.25, 0.875, 0, 2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5]
    l = [1, 1, 1, 1, 2, 2, 2, 3]
    a0 = [8.5740489E-01,  -3.2863233E+00, 1.6480939E+00,  -5.4524817E-02, 6.1623592E-02, 2.7389266E-04,  -6.0655087E-02,
    -3.1811852E-02, -1.1550422E-01, -1.8610466E-02, -1.8348671E-01, 5.5071325E-03, -1.2268039E-02, -5.0433436E-03]
    a1 = [5.6200117E-01, 3.2439544E+00, -4.9628768E+00, -2.2132851E-01, 9.3388356E-02, 2.4940171E-05,  -1.7519204E-01,
    8.9325660E-01, 2.9886613E+00, 1.0881387E-01,  -6.7166746E-01, 1.4477326E-01, -2.8716809E-01, -1.1478402E-01]
    a2 = [-8.1680511E+01, 4.6384732E+02,  -2.7970850E+02, 2.9317364E+01, -2.2324825E+01, -5.0932691E-02, -7.2836590E+00,
    -2.2063100E+02, -3.0435126E+02, 5.8514719E+00,  1.7995451E+02, -1.0178400E+02, 4.0848053E+01,  1.2411984E+01]
    n = a0 .+ a1 .* ω + a2 .* θ
    return EmpiricSingleFluidResidualParam(n,t,d,l)
end

function XiangDeiters(components::Vector{String}; idealmodel=BasicIdeal,
    userlocations=String[],
    ideal_userlocations=String[],
    verbose=false)
    
    @assert length(components) == 1

    params = getparams(components, ["properties/critical.csv", "properties/molarmass.csv"]; userlocations=userlocations, verbose=verbose)
    Pc = params["Pc"]
    Vc = params["Vc"]
    Mw = params["Mw"]
    Tc = params["Tc"]
    ω = params["acentricfactor"]

    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    packagedparams = XiangDeitersParam(Mw,Tc,Pc,Vc,ω)
    Zc = Pc[1]*Vc[1]/(R̄*Tc[1])
    θ = (Zc - 0.29)^2
    consts = XiangDeitersConsts(ω[1],θ)
    references = String["10.1016/j.ces.2007.11.029"]
    model = XiangDeiters(components,packagedparams,init_idealmodel,consts,references)
    return model
end

function a_res(model::XiangDeitersModel,V,T,z = SA[1.0])
    Tc = model.params.Tc[1]
    Vc = model.params.Vc[1]
    N = only(z)
    rho = (N/V)
    δ = rho*Vc
    τ = Tc/T
    return reduced_a_res(model.consts,δ,τ)
end

lb_volume(model::XiangDeitersModel,z = SA[1.0]) = only(model.params.Vc.values)*sum(z)/4

function x0_sat_pure(model::XiangDeitersModel,T,z=SA[1.0])
    #nice trick, if the model has Tc,Pc,acentricfactor, we can do this
    sat_p = LeeKeslerSat(model)
    psat,_,_ = saturation_pressure(sat_p,T)
    
    vl = volume(model,psat,T,phase = :l)
    vv0 = R̄*T/psat
    vv = volume(model,psat,T,vol0 = vv0)
    return vl,vv
end

function x0_volume_liquid(model::XiangDeitersModel,T,z = SA[1.0])
    lb_v = lb_volume(model,z)
    Tc = only(model.params.Tc.values)
    Pc = only(model.params.Pc.values)
    Pc_inv = 1/Pc
    Vc = only(model.params.Vc.values)
    ∑z = only(z)
    Zc = Pc*Vc/(∑z*R̄*Tc)
    Tr = T/Tc
    
    rackket_v = ∑z*R̄*Tc*Pc_inv*Zc^(1+(1-Tr)^(2/7))
    0.4*lb_v + 0.6*rackket_v
end

p_scale(model::XiangDeitersModel,z = SA[1.0]) = dot(model.params.Pc,z)
T_scale(model::XiangDeitersModel,z = SA[1.0]) = dot(model.params.Tc,z)

function crit_pure(model::XiangDeitersModel) 
    single_component_check(crit_pure,model)
    return (model.params.Tc[1],model.params.Pc[1],model.params.Vc[1])
end
function methane_a_ideal(δi,τi)
#gerg v2 part
    n0 = (19.597508817,
    -83.959667892,
    3.00088,
    0.76315,  #C
    0.00460,  #E
    8.74432,  #G
    -4.46921  #I
) #methane
    
    v0 = (4.306474465, #*Tc = D
    0.936220902,
    5.577233895,
    5.722644361)                #methane
    
    a1 = n0[1]
    a2 = n0[2]
    c0 = n0[3]

    #δi = ρ / ρc
    #τi = Tc * Tinv
    n1,n2,n3,n4,n5,n6,n7 = n0
    ϑ1,ϑ2,ϑ3,ϑ4 = v0  
    ai = 0.0
    iszero(n4) || (ai += n4*log(@show abs(sinh(ϑ1*τi))))
    iszero(n5) || (ai -= n5*log(@show cosh(ϑ2*τi)))
    iszero(n6) || (ai += n6*log(@show abs(sinh(ϑ3*τi)))) 
    iszero(n7) || (ai -= n7*log(@show cosh(ϑ4*τi)))

    @show ai
    println()
    Tc = 190.564
    xx = gerg_to_gpe(n0,v0,Tc)

    aj = term_a0_gpe(τi,log(τi),0.,xx.n,xx.t,xx.c,xx.d)
    println()
    @show aj
    return (ai,aj)
end

function gerg_to_gpe(n0,v0,Tc)
    n1,n2,n3,n4,n5,n6,n7 = n0
    v1,v2,v3,v4 = v0.*Tc
    a1,a2,c0 = n1,n2,n3
    #alylee1 = [0.,n4,v1,n5,v2]
    #alylee2 = [0.,n6,v3,n7,v4]
    
    #@show alylee1
    #@show alylee2
    params = [4.00088,0.76315,820.659,0.00460,178.410,8.74432,1062.82,-4.46921,1090.53]
    #params2 = [4.00088,n4,v1,n5,v2,n6,v3,n7,v4]
    #@show norm(params - params2)
    alylee1 = params[1:5]
    alylee2 = vcat(0,params[6:9])
    #@show (alylee2[2:5] -params[6:9])
    n = Float64[]
    t = Float64[]
    c = Float64[]
    d = Float64[]

    for vv in (alylee1,alylee2)
        _,B,C,D,E = vv
        push!(n,B) 
        push!(t,-2C/Tc)
        push!(c,1)
        push!(d,-1)
        #sinh((C/Tc)*tau) = 1 - exp((-2C/Tc)*tau)
        push!(n,-D)
        push!(t,-2E/Tc)
        push!(c,1)
        push!(d,1)
    end

    return (;n,t,c,d,a1,a2,c0)
end

#=
sinh(vi*tau) = (exp(2vi*tau) - 1)/2exp(vi*tau)

log(sinh(vi*tau)) = log((exp(2vi*tau) - 1)) - log(2exp(vi*tau))
= 

=#

