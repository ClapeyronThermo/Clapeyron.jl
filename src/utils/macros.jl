macro groups()
    keys(model.groups)
end

macro components()
    model.components
end

macro f(func, args...)
    args = [esc(arg) for arg in args]
    return :($(f)(model,z,V,T,$(args...)))
end



