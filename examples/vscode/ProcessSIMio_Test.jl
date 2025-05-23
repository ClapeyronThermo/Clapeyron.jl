using ModelingToolkit, DifferentialEquations, Clapeyron, NLsolve
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkit: scalarize
using Suppressor

fluid = ["carbon dioxide","nitrogen","water"]
nfluid = length(fluid)
model = GERG2008(fluid)

Molecular_Weight(model::EoSModel,z) = Clapeyron.molecular_weight(model::EoSModel,z)
@register_symbolic Molecular_Weight(model::EoSModel,z::AbstractVector)

Enthalpy(model::EoSModel,p,T,z) = Clapeyron.enthalpy(model::EoSModel,p,T,z)
@register_symbolic Enthalpy(model::EoSModel,p,T,z::AbstractVector)

Entropy(model::EoSModel,p,T,z) = Clapeyron.entropy(model::EoSModel,p,T,z)
@register_symbolic Entropy(model::EoSModel,p,T,z::AbstractVector)

Density(model::EoSModel,p,T,z) = Clapeyron.mass_density(model::EoSModel,p,T,z)
@register_symbolic Density(model::EoSModel,p,T,z::AbstractVector)

function FlashTp(model::EoSModel,p,T,z)
    verbose = false
    beta,xp,vp,Gsol = Clapeyron.tp_flash_impl(model,p,T,z, HELDTPFlash(verbose = verbose))
    hp = Vector{Float64}(undef,0)
    for ip = eachindex(beta)
        push!(hp, Clapeyron.VT_enthalpy(model,vp[ip],T,xp[ip]))
    end
    sp = Vector{Float64}(undef,0)
    for ip = eachindex(beta)
        push!(sp, Clapeyron.VT_entropy(model,vp[ip],T,xp[ip]))
    end
    res = fill(0.0,(1 + 3*3 + 3*nfluid))
    ires = 1
    res[ires] = length(beta)
    for ip = eachindex(beta)
        ires += 1
        res[ires] =  beta[ip]
    end
    for ip = eachindex(beta)
        for ic = eachindex(xp[ip])
            ires += 1
            res[ires] = xp[ip][ic]
        end
    end
    for ip = eachindex(beta)
        ires += 1
        res[ires] = hp[ip]
    end
    for ip = eachindex(beta)
        ires += 1
        res[ires] = sp[ip]
    end
    return res
end
@register_array_symbolic FlashTp(model::EoSModel,p,T,z::AbstractVector) begin
    size = (1 + 3*3 + 3*nfluid,)
    eltype= eltype(z)
end
#@register_symbolic FlashTp(model::EoSModel,p,T,z::AbstractVector)

function Get_Beta_FlashTp(res,phase)
    np = round(Int64, res[1])
    ires = 1
    beta = Vector{Float64}(undef,0)
    for ip = 1:np
        ires += 1
        push!(beta, res[ires])    
    end
    return beta[phase]
end
@register_symbolic Get_Beta_FlashTp(res::AbstractVector, phase::Int64)

function Get_Enthalpy_FlashTp(res,phase)
    np = round(Int64, res[1])
    ires = 1+np+np*nfluid
    hp = Vector{Float64}(undef,0)
    for ip = 1:np
        ires += 1
        push!(hp, res[ires])
    end
    return hp[phase]
end
@register_symbolic Get_Enthalpy_FlashTp(res::AbstractVector, phase::Int64)

function Get_Entropy_FlashTp(res,phase)
    np = round(Int64, res[1])
    ires = 1+np+np*nfluid+np
    sp = Vector{Float64}(undef,0)
    for ip = 1:np
        ires += 1
        push!(sp, res[ires])
    end
    return sp[phase]
end
@register_symbolic Get_Entropy_FlashTp(res::AbstractVector, phase::Int64)

function Get_Composition_FlashTp(res,phase)
    np = round(Int64, res[1])
    ires = 1+np
    xp = Vector{Vector{Float64}}(undef,0)
    for ip = 1:np
        x = fill(0.0,nfluid)
        for ic = 1:nfluid
            ires += 1
            x[ic] = res[ires]
        end
        push!(xp, x)
    end
    return xp[phase]
end
@register_array_symbolic Get_Composition_FlashTp(res::AbstractVector, phase::Int64) begin
    size = (nfluid,)
    eltype= eltype(res)
end
#@register_symbolic Get_Composition_FlashTp(res::AbstractVector, phase::Int64)

@connector function fluidPort(;name, pg, Tg, zg)
    vars = @variables begin
        p(t), [input = true, description = "Pressure (Pa)", guess = pg]
        T(t), [input = true, description = "Temperature (T)", guess = Tg]
        mdot(t), [input = true, description = "mass flow rate (kg/s)"]
        (z(t))[1:nfluid], [input = true,description  ="Mole fraction vector", guess = zg]
    end
    ODESystem(Equation[], t, vars, [];name=name)
end

@component function MassSource(;name, ps, Ts, mdots, zs)
    @named port = fluidPort(pg = ps, Tg = Ts, zg = zs)
    vars = @variables begin
        s(t)
        T(t)
        p(t)
        h(t)
        (z(t))[1:nfluid], [description  = "Mole fraction vector",  guess = zs]
        (flash(t))[1:(1 + 3*3 + 3*nfluid)], [guess = FlashTp(model,ps,Ts,zs)]
    end

    para = @parameters begin
    #   T_source, [description = "Temperature at source (K)"]
    #   p_source, [description = "pressure at source (Pa)"]
    #   mdot_source, [description = "mass flow rate at source (kg/s)"]
     end

    eqs = [
        scalarize(flash .~ FlashTp(model,ps,Ts,zs))
        scalarize(z .~ zs)
        T ~ Ts
        p ~ ps
        h ~ Get_Enthalpy_FlashTp(flash,1)*Get_Beta_FlashTp(flash,1) + Get_Enthalpy_FlashTp(flash,2)*Get_Beta_FlashTp(flash,2)
        s ~ Get_Entropy_FlashTp(flash,1)*Get_Beta_FlashTp(flash,1) + Get_Entropy_FlashTp(flash,2)*Get_Beta_FlashTp(flash,2)

        port.T ~ T
        port.p ~ p
        port.mdot ~ mdots
        scalarize(port.z .~ z)
    ]

#    eqs = Symbolics.scalarize.(reduce(vcat, Symbolics.scalarize.(eqs)))
    compose(ODESystem(eqs, t, collect(Iterators.flatten(vars)), para;name=name),port)

end

@component function MassSink(;name, ps, Ts, zs)
    @named port = fluidPort(pg = ps, Tg = Ts, zg = zs)
    vars = @variables begin
       s(t)
       T(t)
       p(t)
       h(t)
       (z(t))[1:nfluid], [description  = "Mole fraction vector",  guess = zs]
       (flash(t))[1:(1 + 3*3 + 3*nfluid)], [guess = FlashTp(model,ps,Ts,zs)]
    end
    para = @parameters begin

     end
     eqs = [
        scalarize(flash .~ FlashTp(model,port.p,port.T,port.z))
        scalarize(z .~ port.z)
        T ~ port.T
        p ~ port.p
        h ~ Get_Enthalpy_FlashTp(flash,1)*Get_Beta_FlashTp(flash,1) + Get_Enthalpy_FlashTp(flash,2)*Get_Beta_FlashTp(flash,2)
        s ~ Get_Entropy_FlashTp(flash,1)*Get_Beta_FlashTp(flash,1) + Get_Entropy_FlashTp(flash,2)*Get_Beta_FlashTp(flash,2)

     ]

#     eqs = Symbolics.scalarize.(reduce(vcat, Symbolics.scalarize.(eqs)))
     compose(ODESystem(eqs, t, collect(Iterators.flatten(vars)), para;name=name),port)
end

function PolytropicCompression(model::EoSModel,T_in,p_in,z,πc,η)
    @assert πc >= 1
    @assert η  <= 1
    @assert η  >  0
    mw = Molecular_Weight(model,z)
    d_in = Density(model,p_in,T_in,z)
    p_out = p_in*πc
    h_in = Enthalpy(model,p_in,T_in,z)/mw
    function Compressor(F,x)
        d_out = Density(model,p_out,x[1],z)
        h_out = Enthalpy(model,p_out,x[1],z)/mw
        npoly = log(p_out/p_in)/log(d_out/d_in)
        F[1] = h_out - h_in - npoly/(npoly-1)*(p_out/d_out - p_in/d_in)/η
    end
    Tpoly = nlsolve(Compressor , [T_in])
    return Tpoly.zero[1]
end
@register_symbolic PolytropicCompression(model::EoSModel,T_in::Float64,p_in::Float64,z::Array{Float64, nfluid},πc::Float64,η::Float64)

@component function Compressor(;name, ps, Ts, zs, pd, Td)
    @named inport = fluidPort(pg = ps, Tg = Ts, zg = zs)
    @named outport = fluidPort(pg = pd, Tg = Td, zg = zs)
    vars = @variables begin
       s_in(t)
       T_in(t)
       p_in(t)
       h_in(t)
       (z_in(t))[1:nfluid]

       s_out(t)
       T_out(t)
       p_out(t)
       h_out(t)
       (z_out(t))[1:nfluid]

    end
    para = @parameters begin
        πc = 4.528, [description = "Pressure ratio (-)"]
        η = 0.84, [description = "Polytropic Efficiency (-)"]
    end
    eqs = [

        scalarize(z_in .~ inport.z)
        T_in ~ inport.T
        p_in ~ inport.p
        h_in ~ Enthalpy(model,p_in,T_in,z_in)
        s_in ~ Entropy(model,p_in,T_in,z_in)

        scalarize(z_out .~ z_in)
        p_out ~ πc*p_in
        T_out ~ PolytropicCompression(model,T_in,p_in,z_in,πc,η)
        s_out ~ Entropy(model,p_out,T_out,z_out)
        h_out ~ Enthalpy(model,p_out,T_out,z_out)

        outport.T ~ T_out
        outport.p ~ p_out
        outport.mdot ~ inport.mdot
        scalarize(outport.z .~ z_out)
    ]

#    eqs = Symbolics.scalarize.(reduce(vcat, Symbolics.scalarize.(eqs)))    
    compose(ODESystem(eqs, t, collect(Iterators.flatten(vars)), para;name=name),inport,outport)

end

@component function Separator(;name, pg, Tg, zg)
    @named inport = fluidPort(pg = pg, Tg = Tg, zg = zg)
    @named outport1 = fluidPort(pg = pg, Tg = Tg, zg = zg) # on top
    @named outport2 = fluidPort(pg = pg, Tg = Tg, zg = zg) # on bottom
    vars = @variables begin
       mdot_in(t)
       s_in(t)
       T_in(t)
       p_in(t)
       h_in(t)
       (z_in(t))[1:nfluid]

       mdot_out1(t)
       s_out1(t)
       T_out1(t)
       p_out1(t)
       h_out1(t)
       (z_out1(t))[1:nfluid]

       mdot_out2(t)
       s_out2(t)
       T_out2(t)
       p_out2(t)
       h_out2(t)
       (z_out2(t))[1:nfluid]
       (flash(t))[1:(1 + 3*3 + 3*nfluid)]

    end
    para = @parameters begin end
    eqs = [

        mdot_in ~ inport.mdot

        scalarize(flash .~ FlashTp(model,p_in,T_in,z_in))
        scalarize(z_in .~ inport.z)
        T_in ~ inport.T
        p_in ~ inport.p
        h_in ~ Get_Enthalpy_FlashTp(flash,1)*Get_Beta_FlashTp(flash,1) + Get_Enthalpy_FlashTp(flash,2)*Get_Beta_FlashTp(flash,2)
        s_in ~ Get_Entropy_FlashTp(flash,1)*Get_Beta_FlashTp(flash,1) + Get_Entropy_FlashTp(flash,2)*Get_Beta_FlashTp(flash,2)

        scalarize(z_out1 .~ Get_Composition_FlashTp(flash,1))
        p_out1 ~ p_in
        T_out1 ~ T_in
        h_out1 ~ Get_Enthalpy_FlashTp(flash,1)
        s_out1 ~ Get_Entropy_FlashTp(flash,1)

        scalarize(z_out2 .~ Get_Composition_FlashTp(flash,2))
        p_out2 ~ p_in
        T_out2 ~ T_in
        h_out2 ~ Get_Enthalpy_FlashTp(flash,2)
        s_out2 ~ Get_Entropy_FlashTp(flash,2)
        
        outport1.T ~ T_out1
        outport1.p ~ p_out1
        mdot_out1 ~ mdot_in*Get_Beta_FlashTp(flash,1) 
        outport1.mdot ~ mdot_out1
        scalarize(outport1.z .~ z_out1)

        outport2.T ~ T_out2
        outport2.p ~ p_out2
        mdot_out2 ~ mdot_in*Get_Beta_FlashTp(flash,2)
        outport2.mdot ~ mdot_out2
        scalarize(outport2.z .~ z_out2)
   
    ]

#    eqs = Symbolics.scalarize.(reduce(vcat, Symbolics.scalarize.(eqs)))
    compose(ODESystem(eqs, t, collect(Iterators.flatten(vars)), para;name=name),inport,outport1,outport2)
end

p = 1.4e5
T = 25.0+273.15
zdry=[0.9981, 0.0019]
zwater=0.025
z=append!(zdry*(1.0-zwater),zwater)

FlashResult = FlashTp(model,p,T,z)
println("FlashTp $(FlashResult)")

xpp = Get_Composition_FlashTp(FlashResult,1)
println("FlashTp $(xpp)")

hpp = Get_Enthalpy_FlashTp(FlashResult,1)
println("FlashTp $(hpp)")

spp = Get_Entropy_FlashTp(FlashResult,1)
println("FlashTp $(spp)")

bpp = Get_Beta_FlashTp(FlashResult,1)
println("FlashTp $(bpp)")

@suppress_err begin

    local src, separator, compressor, sink, eqs, systems, flowsheet, sys, u0, para, prob

    @named src = MassSource(ps = p, Ts = T, mdots = 50, zs = z)
    @named separator = Separator(pg = p, Tg = T, zg = z)
    @named compressor = Compressor(ps = p, Ts = T, zs = z, pd = 4.528*p, Td = 161.0+2731.5)
    @named sink =  MassSink(ps = p, Ts = T, zs = z)

    eqs =   [
            connect(src.port,separator.inport)
            connect(separator.outport1,compressor.inport)
            connect(compressor.outport,sink.port)
            ] # Define connections

    systems=[src,separator, compressor,sink] # Define system

    @named flowsheet = System(eqs, t, systems=systems)
    sys = structural_simplify(flowsheet)
    u0 = []
    para =  [   
                sys.compressor.πc =>4.528,
                sys.compressor.η =>0.84 
            ] # boundary conditions and parameters


    prob = SteadyStateProblem(sys,u0,para)
    sol = solve(prob)

    println("Compressor discharge temperature = $(round(sol[compressor.T_out][1]-273.15,digits=2)) deg C")
    comp_h_in = sol[compressor.h_in][1]
    println("Compressor h_in = $(round(comp_h_in,digits=2)) J/mol")
    comp_h_out = sol[compressor.h_out][1]
    println("Compressor h_out = $(round(comp_h_out,digits=2)) J/mol")
    comp_mdot = sol[compressor.inport.mdot][1]
    println("Compressor mdot = $(round(comp_mdot,digits=2)) kg/s")
    comp_mw =  Molecular_Weight(model,sol[compressor.z_out])
    println("Compressor mw = $(round(comp_mw*1000,digits=2)) kg/kmol")
    power = comp_mdot*(comp_h_out - comp_h_in)/comp_mw/1000.0/1000.0
    println("Compressor power = $(round(power,digits=2)) MW")

end

