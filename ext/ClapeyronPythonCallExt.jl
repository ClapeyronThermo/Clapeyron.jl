module ClapeyronPythonCallExt

using Clapeyron
using PythonCall

Clapeyron.format_gccomponents(x::PyList) = begin
    return [
        (xi isa Tuple) ? tuple([
            (xj isa PyDict) ? [k => v  for (k,v) in xj] : 
            (xj isa String) ? xj : 
            throw(error("Conversion failed from Python to Julia! Please open an issue."))
        for xj in xi]...) : xi
    for xi in x]
end

Clapeyron.colnames(x::PyDict) = begin 
    return keys(x)
end
end