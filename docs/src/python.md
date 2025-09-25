# Python package `pyclapeyron`

`pyclapeyron` exposes the functionality of `Clapeyron.jl` to Python via [`juliacall`](https://juliapy.github.io/PythonCall.jl/stable/juliacall/), which provides the bridge to Julia under the hood. Unless otherwise noted, the public Python API mirrors the Julia API described in this documentation.

## Installation

Install from PyPI:

```
pip install pyclapeyron
```

## Example

```python
import pyclapeyron as cl
import numpy as np

pure = cl.SAFTVRMie("butanol")

Tc, pc, vc = cl.crit_pure(pure)
pv, vl, vv = cl.saturation_pressure(pure, 400)

rhol = cl.molar_density(pure, pv, 400, phase="liquid")
rhov = cl.molar_density(pure, pv, 400, phase="vapour")

mix = cl.PCSAFT(["propane", "ethanol"])
pv, vl, vv, _y = cl.bubble_pressure(mix, 300, np.array([0.5,0.5]))
y = np.array(_y)
```

## Notes & limitations

- Inputs: arrays passed to `Clapeyron.jl` functions must be NumPy arrays.
- Outputs: returned arrays are `juliacall` array objects; they can be converted to NumPy arrays with `np.array(...)` (see example above).
- Julia macros and mutating (in-place, `!`) functions are not directly callable from Python.
- Please report any issues on the `Clapeyron.jl` GitHub repository.
