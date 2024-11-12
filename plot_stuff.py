import xarray as xr
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt

ds  = xr.open_dataset('final_state.nc')
irr = ds.irradiance.values
mem = ds.irr_memory.values
avg = ds.irr_mld.values
eps = ds.epsilon.values

#ratio = np.divide(mem - avg, irr - avg)
#eps = np.tile(eps,(10,1))

irr = irr - avg
mem = mem - avg

print(np.shape(eps))
print(np.shape(eps))

ratio = np.zeros((27))
for i in range(0,27):
    pf    = np.polyfit(irr[:,i],mem[:,i],1)
    ratio[i] = pf[0]

plt.scatter(ratio,eps)
plt.xlim(0,1)
plt.ylim(0,1)
plt.show()

print(np.polyfit(ratio,eps,1))
