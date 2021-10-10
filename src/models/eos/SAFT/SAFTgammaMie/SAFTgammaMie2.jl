#just a holder for the z partitions.
#to allow split_model to work correctly
struct γMieZ <: EoSModel
    components::Vector{String}
    z::SingleParam{Vector{Float64}}
end

struct GCSAFTgammaMieParam <: EoSParam
    segment::SingleParam{Int}
    mixedsegment::SingleParam{Vector{Float64}}
    shapefactor::SingleParam{Float64}
    lambda_a::PairParam{Float64}
    lambda_r::PairParam{Float64}
    sigma::PairParam{Float64}
    epsilon::PairParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
end

@registermodel γMieZ

struct SAFTGMie{I,VR} <: SAFTgammaMieModel
    components::Vector{String}
    groups::GroupParam
    sites::SiteParam
    params::GCSAFTgammaMieParam
    idealmodel::I
    vrmodel::VR
    mie_zfractions::γMieZ
    absolutetolerance::Float64
    references::Array{String,1}
end

function SAFTgammaMie2(components; 
    idealmodel=BasicIdeal,
    userlocations=String[],
    ideal_userlocations=String[],
    verbose=false,
    absolutetolerance = 1e-12)

    groups = GroupParam(components, ["SAFT/SAFTgammaMie/SAFTgammaMie_groups.csv"]; verbose=verbose)
    params,sites = getparams(groups, ["SAFT/SAFTgammaMie","properties/molarmass_groups.csv"]; userlocations=userlocations, verbose=verbose)
    components = groups.components
    gc = groups.i_flattenedgroups
    comps = 1:length(components)
    gc_segment = params["vst"]
    shapefactor = params["S"]
    S = shapefactor.values
    vst = gc_segment.values
    v  = groups.n_flattenedgroups
    comp_segment = zeros(Float64,length(comps))
    for i ∈ comps
        res_i = 0.0
        vi = v[i]
        groups_i = groups.i_groups[i]
        for idx in 1:length(groups_i)
            k = groups_i[idx]
            res_i += vi[k]*S[k]*vst[k]
        end
        comp_segment[i] = res_i
    end
    segment = SingleParam("segment",components,comp_segment)

    #used in x_S:
    #x_S(group i) = dot(z,mixsegment[i])/dot(z,m_vr)
    mixsegment =  [[v[i][k]*vst[k]*S[k] for i ∈ comps] for k ∈ gc]
    gc_mixedsegment = SingleParam("mixed segment",components,mixsegment,[false for i in gc],String[],String[])

    gc_sigma = sigma_LorentzBerthelot(params["sigma"])
    gc_sigma.values .*= 1E-10
    gc_epsilon = epsilon_HudsenMcCoubrey(params["epsilon"], gc_sigma)
    ϵ = gc_epsilon.values
    σ = gc_sigma.values
     #helper functions
     function ẑ(i, k)
        return v[i][k]*vst[k]*S[k] / ∑(v[i][l]*vst[l]*S[l] for l ∈ gc)
    end

    zz = [[0.0 for k in 1:length(v[i])] for i in 1:length(v)]
    zzparam = SingleParam("z fraction",components,zz,[false for i in 1:length(components)],String[],String[])
    for i in 1:length(zz)
        zzi = zz[i]
        for k in 1:length(zzi)
            zzi[k] = ẑ(i, k)
        end
    end
    
    function σ̄(i)
        return cbrt(∑(∑(ẑ(i,k)*ẑ(i,l)*σ[k,l]^3 for l ∈ gc) for k ∈ gc))
    end
    
    function σ̄(i, j)
        return (σ̄(i) + σ̄(j))/2
    end

    function ϵ̄(i)
        return ∑(∑(ẑ(i,k)*ẑ(i,l)*ϵ[k,l] for l ∈ gc) for k ∈ gc)
    end
    
    function ϵ̄(i, j)
        if i == j
            return ϵ̄(i)
        else
            return sqrt(σ̄(i)*σ̄(j))/σ̄(i,j) * sqrt(ϵ̄(i)*ϵ̄(j))
        end
    end

    comp_ϵ = [ϵ̄(i, j) for (i,j) in Iterators.product(comps,comps)]
    epsilon = PairParam("epsilon",components,comp_ϵ)
    
    comp_σ = [σ̄(i, j) for (i,j) in Iterators.product(comps,comps)]
    sigma = PairParam("sigma",components,comp_σ)

    gc_lambda_a = lambda_LorentzBerthelot(params["lambda_a"])
    gc_lambda_r = lambda_LorentzBerthelot(params["lambda_r"])
    #@show gc_lambda_a
    gc_λa = gc_lambda_a.values
    gc_λr = gc_lambda_r.values
    
    function λi(ll,i)
        return ∑(∑(ẑ(i,k)*ẑ(i,l)*ll[k,l] for l ∈ gc) for k ∈ gc)
    end
    comp_lambda_a = [λi(gc_λa,i) for i in comps]
    comp_lambda_r = [λi(gc_λr,i) for i in comps]
    lambda_a = lambda_LorentzBerthelot(SingleParam("lambda_a",components,comp_lambda_a))
    lambda_r = lambda_LorentzBerthelot(SingleParam("lambda_r",components,comp_lambda_r))

    comp_epsilon_assoc = AssocParam{Float64}("epsilon assoc",components)
    comp_bondvol = AssocParam{Float64}("bondvol",components)
    
    comp_mw = zeros(Float64,length(components))
    gc_mw = params["Mw"].values

    for i in 1:length(components)
        vi = v[i]
        gi = groups.i_groups[i]
        mwi = 0.0
        for idx in 1:length(gi)
            j = gi[idx]
            mwi += gc_mw[j]*vi[j]
        end
        comp_mw[i] =mwi
    end  
    _mw = SingleParam("molecular_weight",components,comp_mw)
    gcparams = GCSAFTgammaMieParam(gc_segment,gc_mixedsegment, shapefactor,gc_lambda_a,gc_lambda_r,gc_sigma,gc_epsilon,params["epsilon_assoc"],params["bondvol"])
    vrparams = SAFTVRMieParam(segment,sigma,lambda_a,lambda_r,epsilon,comp_epsilon_assoc,comp_bondvol,_mw)
    vr = SAFTVRMie(vrparams, SiteParam(components), BasicIdeal(); ideal_userlocations=ideal_userlocations, verbose=verbose)
    idmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    mie_z = γMieZ(components,zzparam)
    γmierefs = ["10.1063/1.4851455", "10.1021/je500248h"]
    gmie = SAFTGMie(components,groups,sites,gcparams,idmodel,vr,mie_z,absolutetolerance,γmierefs)
    return gmie
end
@registermodel SAFTGMie

function a_res(model::SAFTGMie, V, T, z)
    _data = @f(data)
    vrmodel = model.vrmodel
    _a_hs = a_hs(vrmodel,V,T,z,_data)
    _a_dispchain = a_dispchain(vrmodel,V,T,z,_data)
    return _a_hs+_a_dispchain 
end

function data(model::SAFTGMie, V, T, z)
    _d_gc = @f(d_gc)
    _d = @f(d,_d_gc)
    ζi = @f(ζ0123,_d_gc)
    _ζ_X,σ3x = @f(ζ_X_σ3,_d_gc)
    _ρ_S = @f(ρ_S)
    _ζst = _ζst = σ3x*_ρ_S*π/6
    return (_d_gc,(_d,_ρ_S,ζi,_ζ_X,_ζst,σ3x))
end

function d_gc(model::SAFTGMie, V, T, z)
    _ϵ = model.params.epsilon.diagvalues
    _σ = model.params.sigma.diagvalues
    _λa = model.params.lambda_a.diagvalues
    _λr = model.params.lambda_r.diagvalues
    u = SAFTγMieconsts.u
    w = SAFTγMieconsts.w
    _0 = zero(T+first(z))
    n = length(_ϵ)
    _d_gc = zeros(typeof(_0),n)
    for k ∈ 1:n
        ϵ = _ϵ[k]
        σ = _σ[k]
        λa = _λa[k]
        λr = _λr[k]
        θ = Cλ(model.vrmodel,V,T,z,λa,λr)*ϵ/T
        di = _0
        λrinv = 1/λr
        λaλr = λa/λr
        for j ∈ 1:5
            di += w[j]*(θ/(θ+u[j]))^λrinv*(exp(θ*(1/(θ/(θ+u[j]))^λaλr-1))/(u[j]+θ)/λr)
        end
        _d_gc[k] = σ*(1-di)
    end
    return _d_gc
end

function d(model::SAFTGMie,V,T,z,_d_gc = @f(d_gc))
    _d = zeros(eltype(_d_gc),length(z))
    _0 = zero(eltype(_d))
    _zz = model.mie_zfractions.z.values
    for i ∈ @comps
        _z = _zz[i]
        m = length(_z)
        di = _0
        for k ∈ 1:m
            zk = _z[k]
            iszero(zk) && continue
            dk = _d_gc[k]
            di += zk*zk*dk^3
            for l in 1:k - 1
                di += 0.25*zk*_z[l]*(dk + _d_gc[l])^3
            end
        end
        _d[i] = cbrt(di)
    end
    #d̂[i] = cbrt(∑(∑(@f(ẑ,i,k)*@f(ẑ,i,l)*@f(d,k,l)^3 for l ∈ @groups) for k ∈ @groups))
    return _d
end


function ζ0123(model::SAFTGMie, V, T, z,_d_gc=@f(d_gc))
    vrmodel = model.vrmodel
    m = vrmodel.params.segment.values
    m̄ = dot(z, m)
    m̄inv = 1/m̄
    _0 = zero(V+T+first(z))
    ζ0,ζ1,ζ2,ζ3 = _0,_0,_0,_0
    _ϵ = model.params.epsilon.diagvalues
    σ = model.params.sigma.diagvalues
    _λr = model.params.lambda_r.diagvalues
    _λa = model.params.lambda_a.diagvalues
    u = SAFTVRMieconsts.u
    _mi = model.params.mixedsegment.values
    w = SAFTVRMieconsts.w
    for i ∈ @groups
        λa = _λa[i]
        λr = _λr[i]
        ϵ = _ϵ[i]  
        di =_d_gc[i]
        mi = _mi[i]
        xS = dot(mi,z)*m̄inv
        ζ0 += xS
        ζ1 += xS*di
        ζ2 += xS*di*di
        ζ3 += xS*di*di*di
    end
    c = π/6*N_A*m̄/V
    ζ0,ζ1,ζ2,ζ3 = c*ζ0,c*ζ1,c*ζ2,c*ζ3
    return ζ0,ζ1,ζ2,ζ3 
end

#return π/6*@f(ρ_S)*∑(@f(x_S,k)*@f(d,k)^m for k ∈ @groups)

function ρ_S(model::SAFTGMie, V, T, z)
    return ρ_S(model.vrmodel,V,T,z)
end

function ζ_X_σ3(model::SAFTGMie, V, T, z,_d = @f(d_gc))
    vrmodel = model.vrmodel
    m = vrmodel.params.segment.values
    m̄ = dot(z, m)
    m̄inv = 1/m̄
    σ = model.params.sigma.values
    _mi = model.params.mixedsegment.values
    ρS = N_A/V*m̄
    _ζ_X = zero(first(_d))
    kρS = ρS* π/6
    σ3_x = _ζ_X
    r1 =_ζ_X
    for i ∈ @groups
        mi = _mi[i]
        x_Si = dot(mi,z)*m̄inv
        σ3_x += x_Si*x_Si*(σ[i,i]^3)
        _ζ_X += x_Si*x_Si*(_d[i])^3   
        for j ∈ 1:i-1
            mj = _mi[j]
            x_Sj = dot(mj,z)*m̄inv
            σ3_x += 2*x_Si*x_Sj*(σ[i,j]^3)
            dij = 0.5*(_d[i] + _d[j])
            r1 = 2*x_Si*x_Sj*dij^3          
            _ζ_X += r1
        end
    end

    #return π/6*@f(ρ_S)*∑(@f(x_S,i)*@f(x_S,j)*(@f(d,i)+@f(d,j))^3/8 for i ∈ comps for j ∈ comps)
    return kρS*_ζ_X,σ3_x
end

function  a_hs(model::SAFTGMie, V, T, z,_data = @f(data))
    _,vrdata = _data
    return a_hs(model.vrmodel,V,T,z,vrdata)
end

function a_chain(model::SAFTGMie, V, T, z,_data = @f(data))
    _,vrdata = _data
    return a_chain(model.vrmodel,V,T,z,vrdata)
end

function a_disp(model::SAFTGMie, V, T, z,_data = @f(data))
    groups = @groups
    gcmodel = model
    model = gcmodel.vrmodel
    _d,vrdata = _data
    _,ρS,ζi,_ζ_X,_ζst,_ = vrdata
    comps = @comps
    l = length(comps)
    ∑z = ∑(z)
    m = model.params.segment.values
    _mi = gcmodel.params.mixedsegment.values
    _ϵ = gcmodel.params.epsilon.values
    _λr = gcmodel.params.lambda_r.values
    _λa = gcmodel.params.lambda_a.values
    _σ = gcmodel.params.sigma.values
    m̄ = dot(z, m)
    m̄inv = 1/m̄
    a₁ = zero(V+T+first(z))
    a₂ = a₁
    a₃ = a₁
    _ζst5 = _ζst^5
    _ζst8 = _ζst^8
    _KHS = @f(KHS,_ζ_X,ρS)
    for i ∈ groups
        j = i
        mi = _mi[i]
        x_Si = dot(mi,z)*m̄inv
        ϵ = _ϵ[i,j]
        λa = _λa[i,j]
        λr = _λr[i,j] 
        σ = _σ[i,j]
        _C = @f(Cλ,λa,λr)
        dij = _d[i]
        x_0ij = σ/dij
        dij3 = dij^3
        x_0ij = σ/dij
        #calculations for a1 - diagonal
        aS_1_a = @f(aS_1,λa,_ζ_X)
        aS_1_r = @f(aS_1,λr,_ζ_X)
        B_a = @f(B,λa,x_0ij,_ζ_X)
        B_r = @f(B,λr,x_0ij,_ζ_X)
        a1_ij = (2*π*ϵ*dij3)*_C*ρS*
        (x_0ij^λa*(aS_1_a+B_a) - x_0ij^λr*(aS_1_r+B_r))

        #calculations for a2 - diagonal
        aS_1_2a = @f(aS_1,2*λa,_ζ_X)
        aS_1_2r = @f(aS_1,2*λr,_ζ_X)
        aS_1_ar = @f(aS_1,λa+λr,_ζ_X)
        B_2a = @f(B,2*λa,x_0ij,_ζ_X)
        B_2r = @f(B,2*λr,x_0ij,_ζ_X)
        B_ar = @f(B,λr+λa,x_0ij,_ζ_X)
        α = _C*(1/(λa-3)-1/(λr-3))
        f1,f2,f3,f4,f5,f6 = @f(f123456,α)
        _χ = f1*_ζst+f2*_ζst5+f3*_ζst8
        a2_ij = π*_KHS*(1+_χ)*ρS*ϵ^2*dij3*_C^2 *
        (x_0ij^(2*λa)*(aS_1_2a+B_2a)
        - 2*x_0ij^(λa+λr)*(aS_1_ar+B_ar)
        + x_0ij^(2*λr)*(aS_1_2r+B_2r))
        
        #calculations for a3 - diagonal
        a3_ij = -ϵ^3*f4*_ζst * exp(f5*_ζst+f6*_ζst^2)
        #adding - diagonal
        a₁ += a1_ij*x_Si*x_Si
        a₂ += a2_ij*x_Si*x_Si
        a₃ += a3_ij*x_Si*x_Si
        for j ∈ 1:i-1
            mj = _mi[j]
            x_Sj = dot(mj,z)*m̄inv
            ϵ = _ϵ[i,j]
            λa = _λa[i,j]
            λr = _λr[i,j] 
            σ = _σ[i,j]
            _C = @f(Cλ,λa,λr)
            dij = 0.5*(_d[i]+_d[j])
            x_0ij = σ/dij
            dij3 = dij^3
            x_0ij = σ/dij
            #calculations for a1
            a1_ij = (2*π*ϵ*dij3)*_C*ρS*
            (x_0ij^λa*(@f(aS_1,λa,_ζ_X)+@f(B,λa,x_0ij,_ζ_X)) - x_0ij^λr*(@f(aS_1,λr,_ζ_X)+@f(B,λr,x_0ij,_ζ_X)))

            #calculations for a2
            α = _C*(1/(λa-3)-1/(λr-3))
            f1,f2,f3,f4,f5,f6 = @f(f123456,α)
            _χ = f1*_ζst+f2*_ζst5+f3*_ζst8
            a2_ij = π*_KHS*(1+_χ)*ρS*ϵ^2*dij3*_C^2 *
            (x_0ij^(2*λa)*(@f(aS_1,2*λa,_ζ_X)+@f(B,2*λa,x_0ij,_ζ_X))
            - 2*x_0ij^(λa+λr)*(@f(aS_1,λa+λr,_ζ_X)+@f(B,λa+λr,x_0ij,_ζ_X))
            + x_0ij^(2*λr)*(@f(aS_1,2λr,_ζ_X)+@f(B,2*λr,x_0ij,_ζ_X)))
            
            #calculations for a3
            a3_ij = -ϵ^3*f4*_ζst * exp(f5*_ζst+f6*_ζst^2)
            #adding
            a₁ += 2*a1_ij*x_Si*x_Sj
            a₂ += 2*a2_ij*x_Si*x_Sj
            a₃ += 2*a3_ij*x_Si*x_Sj            
        end
    end
    a₁ = a₁*m̄/T/∑z
    a₂ = a₂*m̄/(T*T)/∑z
    a₃ = a₃*m̄/(T*T*T)/∑z
    #@show (a₁,a₂,a₃)
    adisp =  a₁ + a₂ + a₃ 
    return adisp
end
