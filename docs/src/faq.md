# Frequently Asked Questions

### Why there isn't a default model?

Every model has tradeoffs:
- SAFT equations are good with liquid volume predictions, but struggle with conditions near the critical point.
- Cubics are the fastest, but the liquid density deviates significantly from the experimental value if a translation rule is not used
- Activity coefficient models are not useful for high pressure conditions.
- Reference Helmholtz equations of state (like the ones used by REFPROP) cannot be used to calculate spinodal points, specially far from the critical point, and they have poor extrapolation behaviour.
- A saturation Correlation will be faster for calculating saturation conditions than a helmholtz equation of state.

In summary, the user should consider the conditions and type of properties when selecting an equation of state.

### Why are the default caloric properties so bad?

By Default, almost all models use `BasicIdeal` for their ideal model. 
This model assumes an isobaric heat capacity of `Cp = 2.5R` (the theoretical value for a monoatomic gas). 
If an user is only looking to calculate pressure‑dependent properties (equilibria, saturation, critical points, etc), the `BasicIdeal` model will not ask for additional parameters.
Note that `SingleFluid` models (and by extension, `MultiFluid` models) have a reference ideal model included.

### Why my flash/equilibria calculation failed?

There are a lot of motives, but a lot of those can be attributed by a bad initial point.
While we try to provide good initial points for our existing equilibria methods in a best‑effort basis, conditions far from the ideal (like azeotropes) will cause equilibria methods to fail if the initial point is too far from the solution. 
We reccommend to always try to provide an initial point for multicomponent properties, if that initial point is available.

### Why is my equation of state so slow to evaluate?

There are two main causes for a slow function in julia:
- Performance Problems: For example, mutating a global variable, or changing the type of a variable mid-computation. 
  Check [Julia's Performance Tips](https://docs.julialang.org/en/v1/manual/performance-tips/). In the `Clapeyron.jl` context, type instability is the main cause of performance problems, due to how derivatives are calculated (via `ForwardDiff.jl` dual numbers). 
  If the defined function is not stable for multiple number types, then it will be slower to evaluate:

  ```julia
  struct MyvdW <: EoSModel
      a::Float64
      b::Float64
  end

  #bad
  function Clapeyron.a_res(model::MyvdW,V,T,z)
      result = 0 # an integer type
      n = sum(z)
      a = model.a
      b = model.b
      rho  = n/V
      result =  result - log1p(b*rho) #will change to a value depending of the value of rho
      result = result - a*rho/8.314/T #will change to a value depending of the value of T
      return result
  end

  #good 
  function Clapeyron.a_res(model::MyvdW,V,T,z)
      result = zero(Base.promote_eltype(model,V,T,z)) #the result value already considers the types of the model, V,T and the input amounts
      n = sum(z)
      a = model.a
      b = model.b
      rho  = n/V
      result =  result - log1p(b*rho)
      result = result - a*rho/8.314/T
      return result
  end
  ```
  While julia is really good at inferring types in arithmetic operations (in the example above, `return log1p(b*rho) - a*rho/8.314/T` would work just as well), writting an equation of state in a type-stable way allows faster compilation and evaluation.

- Algorithm Problems: The most common of those problems is calculating a variable over an over again in a loop. 
  We recommend calculating everything that could be needed beforehand and passing those calculated variables along:

  ```julia
  function data(model::MyModel,V,T,z)
      d = Clapeyron.d(model,V,T,z)
      n = sum(z)
      m = LinearAlgebra.dot(model.params.segment.values,z)/n #∑mᵢzᵢ
      return (m,d,n)
  end

  function a_res(model::MyModel,V,T,z,mydata = data(model,V,T,z))
      return a1(model,V,T,z,mydata) + a2(model,V,T,z,mydata)
  end

  function a1(model,V,T,Z,mydata = data(model,V,T,z))
      m,d,n = mydata
      #definition for a1
  end

  function a2(model,V,T,Z,mydata = data(model,V,T,z))
      _,d,_ = mydata #in this function, only d is used
      #definition for a2
  end
  ```
