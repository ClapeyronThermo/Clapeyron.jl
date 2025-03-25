# v0.6.10

## New Features

- Experimental: support for using ForwardDiff through some equilibria procedures in some restricted cases:
  - `saturation_temperature`
  - `saturation_pressure`
  - `bubble_pressure` (only helmholtz EoS, without non-volatiles)
  - `bubble_temperature` (only helmholtz EoS, without non-volatiles)
  - `dew_pressure` (only helmholtz EoS, without non-condensables)
  - `dew_temperature` (only helmholtz EoS, without non-condensables)
  - `tp_flash` (only helmholtz EoS, without non-condensables nor non-volatiles)
  - `Pproperty`
  - `TProperty`
  - `X-Y` flashes (single component helmholtz EoS)

  For example, this function now works without the need to propagate dual information through any iterative solvers:

  ```julia
  function ad_function(model,Q)

    z = Clapeyron.SA[1.0]
    p = 1e5
    T = 350.0
    k = 10
    h1 = 6200.0
    T1 = PH.temperature(model,p,h1,z)
    h2 = h1 + Q
    T2 = PH.temperature(model,p,h2,z)
    return -Q + (T - (T1+T2)/2) * k
  end
  model = PR("air")
  ForwardDiff.derivative(Base.Fix1(ad_function,model),500.0)
  ```

## Bug Fixes

- Fix typo in composition return `ChemPotDewTemperature`
- incorrect scaling for second virial coefficient in the case of cubics
- fixing support for second order properties in activity + puremodel EoS calculation
