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
    conditions = (p = p, T = T,z = z)
    β = M.flash_2ph!(S,K,model,conditions;method = method.method,method.kwarg...)
    x = M.liquid_mole_fraction!(S.x, z, K, β)
    y = M.vapor_mole_fraction!(S.y, x, K)
    g = C.__tpflash_gibbs_reduced(model,p,T,x,y,β,:vle)
    comps = [x,y]
    βi = [1-β,β]
    volumes = [volume(model,p,T,x,phase = :l),volume(model,p,T,y,phase = :v)]
    return comps,βi,volumes,g
end
