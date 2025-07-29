##
## Support for Clapeyron.tp_flash interface
## MCFlashJL is defined as a stub in Clapeyron, here we define the initializer.
##
function C.MCFlashJL(;method = M.SSIFlash(),
    storage = nothing,
    V = NaN,
    K = nothing,
    tolerance = 1e-8,
    maxiter = 20000,
    verbose = false,
    check = true,
    z_min = nothing,
    update_forces = true)
    kwargs = (;tolerance,maxiter,verbose,check,z_min,update_forces)
    return C.MCFlashJL(method,storage,K,V,kwargs)
end

function M.flash_storage(model::C.EoSModel,p,T,z,method::C.MCFlashJL)
    cond = (p = p,T = T,z = z)
    mcf_method = method.method
    M.flash_storage(model,cond;method = mcf_method,method.kwarg...)
end

function C.tp_flash_impl(model,p,T,z,method::C.MCFlashJL)
    S = method.storage == nothing ? M.flash_storage(model,p,T,z,method) : method.storage
    K = method.K === nothing ? C.wilson_k_values!(zeros(typeof(p+T+one(eltype(model))),length(model)),model,p,T,S.crit) : method.K
    conditions = (p = p, T = T, z = z)
    res = M.flashed_mixture_2ph!(S, model, conditions, K; method.kwarg...)
    #β = M.flash_2ph!(S,K,model,conditions;method = method.method,method.kwarg...)
    M.flashed_mixture_2ph
    KK = res.K
    if res.state == M.single_phase_v
        β = one(eltype(KK))
    elseif res.state == M.single_phase_l
        β = zero(eltype(KK))
    elseif res.state == M.two_phase_lv
        β = res.V
    else
    end
    n = sum(z)
    x,y = res.liquid.mole_fractions,res.vapor.mole_fractions
    Zl,Zv = res.liquid.Z,res.vapor.Z
    nRT = n*C.Rgas(model)*T
    vl,vv = Zl*nRT/p,Zv*nRT/p
    volumes = [vl,vv]
    data = C.FlashData(p,T)
    comps = [x,y]
    βi = [n*(1-β),n*β]
    return C.FlashResult(comps,βi,volumes,data)
end
