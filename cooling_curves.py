import numpy as np
import matplotlib.pyplot as plt
import pathlib
home = str(pathlib.Path.home())
import os
plt.style.use('science')

data = np.loadtxt(os.path.join(os.getcwd(),'Eatson_cooling_curve_solar_logT4-9.dat'))
data2 = np.loadtxt(os.path.join(os.getcwd(),'Eatson_cooling_curve_WC_logT4-9.dat'))


fig, ax = plt.subplots(figsize=(5,3.5))
ax.plot(10**data[:,0],data[:,1],label='Solar')
ax.plot(10**data2[:,0],data2[:,1],label='WC')
ax.set_xlabel('Temperature (K)')
ax.set_ylabel(r'$\Lambda(T)$' + ' (erg cm$^{-3}$ s$^{-1}$)')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(10**4,10**9)
ax.legend()
fig.savefig(os.path.join(os.getcwd(),'cooling_curve.png'), dpi=300, bbox_inches='tight')