from juliacall import Main as jl

jl.seval("using Clapeyron")

for n in jl.names(jl.Clapeyron):
    s = str(n)
    valid = not s[0] == '@' and not s[-1] == '!'
    if valid:
        exec(f'{n} = jl.{n}')
