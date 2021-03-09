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

molecular_weight(model::CubicModel,z = @SVector [1.]) = 0.001*mapreduce(+,*,paramvals(model.params.Mw),z)

function x0_volume_sc(model::ABCubicModel,p,T,z)
    Zc = cubic_zc(model)
    return Zc*R̄*T/p
end

