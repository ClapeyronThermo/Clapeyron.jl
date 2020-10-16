macro groups(component)
    return :($(esc(:(keys(model.group_multiplicities[$(component)])))))
end

macro groups()
    return :($(esc(:(model.groups))))
end

macro comps()
    return :($(esc(:(model.components))))
end

macro f(func, args...)
    args = [esc(arg) for arg in args]
    return :($(func)($(esc(:model)),$(esc(:N)),$(esc(:V)),$(esc(:T)),$(args...)))
end



