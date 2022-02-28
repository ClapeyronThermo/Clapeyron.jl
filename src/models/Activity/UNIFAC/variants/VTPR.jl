abstract type VTPRUNIFACModel <: UNIFACModel end


struct VTPRUNIFACCache <: EoSModel
    components::Vector{String}
    m::Vector{Float64}
end

struct VTPRUNIFAC{c<:EoSModel} <: VTPRUNIFACModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    groups::GroupParam
    params::UNIFACParam
    puremodel::Vector{c}
    references::Array{String,1}
    unifac_cache::VTPRUNIFACCache
end

@registermodel VTPRUNIFAC
export VTPRUNIFAC

function VTPRUNIFAC(components; puremodel=PR,
    userlocations=String[], 
     verbose=false)
    groups = GroupParam(components, ["Activity/UNIFAC/VTPR/VTPR_groups.csv"]; verbose=verbose)

    params = getparams(groups, ["Activity/UNIFAC/VTPR/VTPR_like.csv", "Activity/UNIFAC/VTPR/VTPR_unlike.csv"]; userlocations=userlocations,  asymmetricparams=["A","B","C"], ignore_missing_singleparams=["A","B","C"], verbose=verbose)
    A  = params["A"]
    B  = params["B"]
    C  = params["C"]
    Q  = params["Q"]
    R = deepcopy(Q)
    R.values .= 0
    icomponents = 1:length(components)
    cache = VTPRUNIFACCache(groups)
    init_puremodel = [puremodel([groups.components[i]]) for i in icomponents]
    gc_mixedsegment = mix_segment(groups)
    packagedparams = UNIFACParam(A,B,C,R,Q,gc_mixedsegment)
    references = String[]
    model = VTPRUNIFAC(components,icomponents,groups,packagedparams,init_puremodel,references,cache)
    return model
end

function VTPRUNIFACCache(groups::GroupParam)
    comps = 1:length(groups.components)
    comp_segment = zeros(length(comps))
    v = groups.n_flattenedgroups
    for i ∈ comps
        res_i = 0.0
        vi = v[i]
        groups_i = groups.i_groups[i]
        for idx ∈ 1:length(groups_i)
            k = groups_i[idx]
            res_i += vi[k]
        end
        comp_segment[i] = res_i
    end
    m = comp_segment
    return VTPRUNIFACCache(groups.components,m)
end

function excess_gibbs_free_energy(model::VTPRUNIFACModel,V,T,z)
    return excess_g_res(model,p,T,z)*R̄*T
end