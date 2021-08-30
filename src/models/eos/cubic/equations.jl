function cubic_ab(model::ABCubicModel,V,T,z=SA[1.0],n=sum(z))
    invn2 = (one(n)/n)^2
    a = model.params.a.values
    b = model.params.b.values
    α = @f(α_function,model.alpha)
    c = @f(translation,model.translation)
    if length(z)>1
        ā,b̄,c̄ = @f(mixing_rule,model.mixing,α,a,b,c)
    else
        ā = a[1,1]*α[1]
        b̄ = b[1,1]
        c̄ = c[1,1]
    end
    return ā ,b̄, c̄
end

function second_virial_coefficient(model::ABCubicModel,T,z = SA[1.0])
    #@info "fast shortcut"
    a,b,c = cubic_ab(model,1e-4,T,z)
    return b-a/(R̄*T)
end

# function x0_volume_sc(model::ABCubicModel,p,T,z)
#     Zc = cubic_zc(model)
#     return Zc*R̄*T/p
# end

function volume(model::ABCubicModel,p,T,z=SA[1.0];phase=:unknown,threaded=false)
    lb_v   =lb_volume(model,z)
    xx = z/sum(z)
    RTp = R̄*T/p
    _poly,c̄ = cubic_poly(model,p,T,z)

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
        return vl-c̄
    end
    if phase == :unknown
        gg = gibbs(vg-c̄)
        gl = gibbs(vl-c̄)
        #@show vg,vl
        #@show gg,gl
        return ifelse(gg<gl,vg-c̄,vl-c̄)
    elseif is_liquid(phase)
        return vl-c̄
    elseif is_vapour(phase)
        return vg-c̄
    # else
    #     gg = gibbs(vg)
    #     gl = gibbs(vl)
    #     #@show vg,vl
    #     #@show gg,gl
    #     return ifelse(gg<gl,vg,vl)
    end
end

function ab_consts(model::CubicModel)
    return ab_consts(typeof(model))
end


    


