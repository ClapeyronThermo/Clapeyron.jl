"""
    kij_mix(f,p::ClapeyronParam,k::PairParam)::PairParam
    kij_mix(f,p::ClapeyronParam)::PairParam

General combining rule for pair parameter with a `kᵢⱼ` interaction parameter. returns a pair parameter with non diagonal entries equal to:
```
pᵢⱼ = f(pᵢ,pⱼ,kᵢⱼ)
```
Where `f` is a 'combining' function that follows the rules:
```
f(pᵢ,pⱼ,0) = f(pⱼ,pᵢ,0)
f(pᵢ,pᵢ,0) = pᵢ
```
and `k` must follow:
```
kᵢᵢ = 0 
```
Ignores non-diagonal entries already set.

If a Single Parameter is passed as input, it will be converted to a Pair Parameter with `pᵢᵢ = pᵢ`.
"""
function kij_mix(f::F,P::SingleOrPair,K = nothing) where F
    param = PairParam(P.name,P.components,P.values,P.ismissingvalues,P.sourcecsvs,P.sources)
    return kij_mix!(f,param,K)
end


"""
    pair_mix(g,P::ClapeyronParam,Q::ClapeyronParam)::PairParam
    pair_mix(g,P::ClapeyronParam,Q::ClapeyronParam)::PairParam

General combining rule for a pair and a single parameter. returns a pair parameter `P` with non diagonal entries equal to:
```
Pᵢⱼ = g(Pᵢ,Pⱼ,Qᵢ,Qⱼ,Qᵢⱼ)
```
Where `f` is a 'combining' function that follows the rules:
```
Pᵢⱼ = Pⱼᵢ = g(Pᵢ,Pⱼ,Qᵢ,Qⱼ,Qᵢⱼ) = g(Pⱼ,Pᵢ,Qⱼ,Qᵢ,Qᵢⱼ)
g(Pᵢ,Pᵢ,Qᵢ,Qᵢ,Qᵢ) = Pᵢ
```
it is a more general form of `kij_mix`, where `kij_mix(f,P,Q) == pair_mix(g,P,Q)` is correct if:
```
f(Pᵢ,Pⱼ,Qᵢⱼ) = g(Pᵢ,Pⱼ,_,_,Qᵢⱼ)

If you pass a `SingleParam` or a vector as input for `Q`, then `Qᵢⱼ` will be considered 0.
```
"""
function pair_mix(f::F,P::SingleOrPair,Q::SingleOrPair) where F
    param = PairParam(P.name,P.components,P.values,P.ismissingvalues,P.sourcecsvs,P.sources)
    return pair_mix!(f,param,Q)
end
