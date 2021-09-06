#=struct SAFTgammaMieParam <: EoSParam
    segment::SingleParam{Int} #done
    shapefactor::SingleParam{Float64} #aux
    lambda_a::PairParam{Float64} #done
    lambda_r::PairParam{Float64} #done
    sigma::PairParam{Float64} #done
    epsilon::PairParam{Float64} #done
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
end
=#
function SAFTgammaMie2(components; idealmodel=BasicIdeal, userlocations=String[], ideal_userlocations=String[], verbose=false)
    groups = GroupParam(components, ["SAFT/SAFTgammaMie/SAFTgammaMie_groups.csv"]; verbose=verbose)
    params,sites = getparams(groups, ["SAFT/SAFTgammaMie"]; userlocations=userlocations, verbose=verbose)
    components = groups.components
    gc = model.groups.i_flattenedgroups
    comps = 1:length(components)
    _segment = params["vst"]
    shapefactor = params["S"]
    S = shapefactor.values
    vst = _segment.values
    v  = model.groups.n_flattenedgroups
    segment = SingleParam("segment",components,[∑(v[i][k]*vst[k]*S[k] for k ∈ gc[i]) for i in comps])
    
    
    
    function ẑ(i,k) 
        return v[i][k]*vst[k]*S[k] / ∑(v[i][l]*vst[l]*S[l] for l ∈ gc[i])
    end

    _epsilon = epsilon_HudsenMcCoubrey(params["epsilon"], sigma)
    params["sigma"].values .*= 1E-10
    _sigma = sigma_LorentzBerthelot(params["sigma"])
    
    ϵ = _epsilon.values
    
    function σ̄(i)
        σ = model.params._sigma.values
        return cbrt(∑(∑(ẑ(i,k)*ẑ,(i,l)*σ[k,l]^3 for l ∈ gc) for k ∈ gc))
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
            return sqrt(@f(σ̄,i)*@f(σ̄,j))/@f(σ̄,i,j) * sqrt(ϵ̄(i)*ϵ̄(j))
        end
    end

    comp_ϵ = [ϵ̄(i, j) for (i,j) in Iterators.product(comps,comps)]
    epsilon = PairParam("epsilon",components,comp_ϵ)
    
    comp_σ = [σ̄(i, j) for (i,j) in Iterators.product(comps,comps)]
    epsilon = PairParam("sigma",components,comp_σ)


    gc_λa = lambda_LorentzBerthelot(params["lambda_a"]).values
    gc_λr = lambda_LorentzBerthelot(params["lambda_r"]).values
    
    comp_lambda_a = [∑(∑(ẑ(i,k)*ẑ(j,l)*gc_λa[k,l] for l ∈ gc) for k ∈ gc) for (i,j) in Iterators.product(comps,comps)]
    comp_lambda_r = [∑(∑(ẑ(i,k)*ẑ(j,l)*gc_λr[k,l] for l ∈ gc) for k ∈ gc) for (i,j) in Iterators.product(comps,comps)]

    lambda_a = PairParam("lambda_a",components,comp_lambda_a)
    lambda_r = PairParam("lambda_r",components,comp_lambda_r)
    
    
    
    epsilon_assoc = params["epsilon_assoc"]
    bondvol = params["bondvol"]
   #=

   SEGMENT
    #by comp
    x = z/∑(z)
    m = model.params.segment.values
    return -∑(x[i]*(log(@f(g_Mie,i))*(m[i]-1)) for i ∈ @comps)

    #by GC
    x = z/∑(z)
    v  = model.groups.n_flattenedgroups
    vst = model.params.segment.values
    S = model.params.shapefactor.values
    return -∑(x[i] * (∑(v[i][k]*vst[k]*S[k] for k ∈ @groups(i))-1) * log(@f(g_Mie,i)) for i ∈ @comps)
    
    LAMBDA
    function λ̄a(model::SAFTgammaMieModel, V, T, z, i)
        λa = model.params.lambda_a.values
        return ∑(∑(@f(ẑ,i,k)*@f(ẑ,i,l)*λa[k,l] for l ∈ @groups) for k ∈ @groups)
    end

    function λ̄r(model::SAFTgammaMieModel, V, T, z, i)
        λr = model.params.lambda_r.values
        return ∑(∑(@f(ẑ,i,k)*@f(ẑ,i,l)*λr[k,l] for l ∈ @groups) for k ∈ @groups)
    end

    
    =#
    #sites = SiteParam(Dict("e1" => params["n_e1"], "e2" => params["n_e2"], "H" => params["n_H"]))
    packagedparams = SAFTgammaMieParam(segment, shapefactor, lambda_a, lambda_r, sigma, epsilon, epsilon_assoc, bondvol)
    references = ["10.1063/1.4851455", "10.1021/je500248h"]

    #model = SAFTgammaMie(packagedparams, groups, sites, idealmodel; ideal_userlocations=ideal_userlocations, references=references, verbose=verbose)
    #return model
end