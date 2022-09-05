abstract type CPCRuleModel <: MixingRule end

function omega_fit(Tc,Pc,m,Δ1,Δ2,g)
    #http://dx.doi.org/10.1021/acs.iecr.9b00436
    #supplementary information
    #should work with any cubic EoS with constant Δ1,Δ2 at the moment
    δ1,δ2 = -Δ1,-Δ2
    function ∂lng∂β(β)
        _g,_∂g = Solvers.f∂f(g,β)
        return _∂g/_g
    end

    function f0(β)
        ∂1lng, ∂2lng, ∂3lng = Solvers.f∂f∂2f(∂lng∂β,β)
        mdm = ((1-m)/m)
        k1 = 2*(2*β - 1)/(1 - β)^3
        k2 = mdm*(β*β*∂2lng + 2*β*∂1lng + 1) + 1/(1 - β)^2
        f0 = (1 + δ1*β)*(1 + δ2*β)
        f1 = β/f0
        f2 = (1 - δ1*δ2*β*β)/(f0*f0)
        f3 = (β*(δ1*δ2*β)^2 - 3*δ1*δ2*β - δ1 - δ2)/(f0*f0*f0)
        k3 = 2/(β*f2 + f1)
        k4 = (β*β*f3 - f1)
        k5 = mdm*(β*β*β*∂3lng + 2*β*β*∂2lng - 2*β*∂1lng - 2)
        return k1 - k2*k3*k4 + k5
    end

    prob = Roots.ZeroProblem(f0,0.5) #TODO, look for better initial point
    βc = Roots.solve(prob)
    f0 = (1 + δ1*βc)*(1 + δ2*βc)
    f1 = βc/f0
    f2 = (1 - δ1*δ2*βc*βc)/(f0*f0)
    f3 = (βc*(δ1*δ2*βc)^2 - 3*δ1*δ2*βc - δ1 - δ2)/(f0*f0*f0)

    ∂1lng, ∂2lng, ∂3lng = Solvers.f∂f∂2f(∂lng∂β,βc)
    λmon = - f1 + (1 - βc)*(βc*f2 + f1)
    λchain = 1 + βc*∂1lng + (1 - βc)*(-βc*βc*∂2lng - 2*βc*∂lng - 1)
    Zc_mon = 1/((1 - βc)*(1 + f1/λmon))
    Zc_chain = (f1*λchain/λmon + 1 + βc*∂1lng)/(1 + f1/λmon)
    Zc = m*Zc_mon + (1 - m)*Zc_chain
    Ωa = (1/m^2)*(βc*Zc*Zc/λmon + (m - 1)*βc*Zc*λchain/λmon)
    Ωb = βc*Zc/m
    return Ωa,Ωb
end

function ab_premixing(::Type{AB},mixing::CPCRuleModel,Tc,pc,kij) where AB <: <:ABCubicModel
    Ωa, Ωb = ab_consts(AB)
    components = pc.components
    n = length(components)
    m = mixing.params.segment
    omega_a = mixing.params.omega_a
    omega_b = mixing.params.omega_b
    g(β) = g_rdf(mixing.rdf,β)
    for i in 1:n
        Ωa,Ωb = omega_fit(Tc[i],Pc[i],m,Δ1,Δ2,g)
        omega_a[i] = Ωa
        omega_b[i] = Ωb
    end
    ai = omega_a .* R̄^2*Tc^2/pc
    bi = omega_b*R̄*Tc/pc

    a = epsilon_LorentzBerthelot(SingleParam("a",components,ai),kij)
    b = sigma_LorentzBerthelot(SingleParam("b",components,bi))
    return a,b
end


include("mixing/CPCRule.jl")
