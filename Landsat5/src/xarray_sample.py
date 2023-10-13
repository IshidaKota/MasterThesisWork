import numpy as np
import pandas as pd
import xarray as xr

data = xr.DataArray(np.random.randn(2, 3), dims=("x", "y"), coords={"x": [10, 20],"y": [10, 20,30]})

print(data)

print(f"data.values={data.values}")
print(data.dims)
print(data.coords)
print(data.attrs)
print(data.loc[10])