function eos(model::ABCubicModel, V, T, z=@SVector [1.0])
    n = sum(z)
    n⁻¹ = 1/n   
    x = z.*n⁻¹
    v = V/n
    return R̄*n*T * (a_ideal(idealmodel(model),V,T,z)+ a_resx(model,v,T,x))
end

function a_res(model::ABCubicModel, V, T, z=@SVector [1.0])
    n = sum(z)
    n⁻¹ = 1/n   
    x = z.*n⁻¹
    v = V*n⁻¹
    return a_resx(model,v,T,x)
end

#fast shortcut to evaluate cubics, pressure is known
function ∂f∂v(model::ABCubicModel,v,t,z)
    #@info "fast shortcut"
    a,b,p = cubic_abp(model,v,t,z)
    return -p
end

function x0_volume_sc(model::ABCubicModel,p,T,z)
    Zc = cubic_zc(model)
    return Zc*R̄*T/p
end

function volume(model::ABCubicModel,p,T,z=SA[1.0];phase="unknown")
    lb_v   =lb_volume(model,z)
    xx = z/sum(z)
    RTp = R̄*T/p
    _poly = cubic_poly(model,p,T,z)

    #sols = Solvers.poly3(_poly)
    sols = Solvers.roots3(_poly)
    #xa = real.(sols) .* RTp
    #@show xa
    #with this, real sol is first
    function imagfilter(x)
        absx = abs(imag(x)) 
        return absx <8*eps(typeof(absx))
    end
    
    isreal = map(imagfilter,sols)
    _sols = sort(real.(sols)) .* RTp
    
    #@show _sols[isreal]
    #@show model.params.b.values
    #@show sum(isreal)
    if sum(isreal) == 3 #3 roots
        sols = sort(real.(sols))
        vg = last(sols)*RTp
        _vl = first(sols)*RTp
        vl = ifelse(_vl>lb_v,_vl,vg)
    elseif  sum(isreal) == 1
        i = findfirst(imagfilter,sols)
        vl = real(sols[i])*RTp
        vg = real(sols[i])*RTp
    elseif  sum(isreal) == 0
        @show sols
        error("weird")
    end
 
    #fp(_v) = pressure(model,_v,T,z) - p
    #@show vl
    #@show vg
    function gibbs(v)
        _df,f =  ∂f(model,v,T,z)
        dv,dt = _df
        if abs((p+dv)/p) > 0.03
            return Inf
        else
            return f  +p*v
        end
    end
    #this catches the supercritical phase as well
    if vl ≈ vg
        return vl
    end
    if phase == "unknown"
        gg = gibbs(vg)
        gl = gibbs(vl)
        #@show vg,vl
        #@show gg,gl
        return ifelse(gg<gl,vg,vl)
    elseif is_liquid(phase)
        return vl
    elseif is_vapour(phase)
        return vg
    else
        gg = gibbs(vg)
        gl = gibbs(vl)
        #@show vg,vl
        #@show gg,gl
        return ifelse(gg<gl,vg,vl)
    end
end

function ab_consts(model::CubicModel)
    return ab_consts(typeof(model))
end


    


