import numpy as np
import matplotlib.pyplot as plt
import sympy as smp
plt.style.use('science')
import os
import pathlib
home = str(pathlib.Path.home())

image_dir = os.path.join(home, "code/project/scripts/images")

def chi(v, d, Mdot):
    """
    Defines the chi parameter for a given wind velocity, distance, and mass loss rate.
    Reference: Stevens et al. 1992 
    """
    v_8 = v/1000
    d = d/1e12
    Mdot = Mdot/1e-7
    return (v_8**4 * d / Mdot)

def velocity(d, v_inf, R):
    """
    Calculates the wind velocity at a given distance from the star.
    """
    return v_inf * (1 - R/d)

# Defining values for WR140 system
D = 2.042e13 # Distance between stars at periastron in cm 
vO = 3100 # Velocity of O star wind in km/s
vWR = 2860 # Velocity of WR star wind in km/s
MdotO = 3.7e-7 # Mass loss rate of O star in Msun/yr
MdotWR = 1.7e-5 # Mass loss rate of WR star in Msun/yr
R_O = 1.74e12 # Radius of O star in cm
R_WR = 1.39e11  # Radius of WR star in cm

phi = np.sqrt(MdotO*vO/(MdotWR*vWR)) # Ratio of wind momentum rates

D_array = np.linspace(D, 5*D, 1000) 

r_WR_array = 1/(1+phi) * D_array 
r_O_array = D_array - r_WR_array

vels_O = velocity(r_O_array, vO, R_O)
chis_O = chi(vels_O, r_O_array, MdotO)

vels_WR = velocity(r_WR_array, vWR, R_WR)
chis_WR = chi(vels_WR, r_WR_array, MdotWR)

x = np.linspace(0, 10*D, 1000)


fig, ax = plt.subplots(3,1, figsize=(5,8))
ax[0].plot(D_array/D, r_O_array/D, label="O star")
ax[0].plot(D_array/D, r_WR_array/D, label="WR star")
ax[0].set_ylabel("$d_{12}/d_{sep}$")
ax[0].legend()

ax[1].plot(D_array/D, vels_O/1000, label="O star", color="C0")
ax[1].plot(D_array/D, vels_WR/1000, label="WR star", color="C1")
# ax[1].plot(D_array/D, np.ones(len(D_array))*vO/1000, color="C0", linestyle="--")
# ax[1].plot(D_array/D, np.ones(len(D_array))*vWR/1000, color="C1", linestyle="--")
ax[1].set_ylabel("$V_{wind}$ $(10^3 km/s)$")
# ax[1].yaxis.set_major_formatter(plt.FormatStrFormatter(''))
ax[1].legend(loc="lower right")


ax[2].plot(D_array[np.where(chis_O<15)]/D, chis_O[np.where(chis_O<15)], label="O star")
ax[2].plot(D_array[np.where(chis_WR<15)]/D, chis_WR[np.where(chis_WR<15)], label="WR star")
ax[2].set_xlabel("$D$ $(d_{sep})$")
ax[2].set_ylabel("$\chi(D)$")
ax[2].legend(loc="lower right")

fig.savefig(os.path.join(image_dir, "wind_velocity.png"), bbox_inches='tight', dpi=500)

# fig.savefig(os.path.join(image_dir, "D_vs_d12.png"), bbox_inches='tight', dpi=300)
# fig2.savefig(os.path.join(image_dir, "D_vs_vwind.png"), bbox_inches='tight', dpi=300)
# fig3.savefig(os.path.join(image_dir, "D_vs_chi.png"), bbox_inches='tight', dpi=300)

# fig, ax = plt.subplots(figsize=(5,3.5))
# ax.plot(D_array, r_O_array, label="O star")
# ax.plot(D_array, r_WR_array, label="WR star")
# ax.set_xlabel("$D$ (cm)")
# ax.set_ylabel("$d_{12}$ (cm)")
# ax.legend(frameon=True)

# fig2, ax2 = plt.subplots(figsize=(5,3.5))
# ax2.plot(D_array, vels_O, label="O star")
# ax2.plot(D_array, vels_WR, label="WR star")
# ax2.set_xlabel("$D$ (cm)")
# ax2.set_ylabel("$V_{wind}$ (km/s)")
# ax2.legend(frameon=True)

# fig3, ax3 = plt.subplots(figsize=(5,3.5))
# ax3.plot(D_array[np.where(chis_O<15)], chis_O[np.where(chis_O<15)], label="O star")
# ax3.plot(D_array[np.where(chis_WR<15)], chis_WR[np.where(chis_WR<15)], label="WR star")
# ax3.set_xlabel("$D$ $(cm)$")
# ax3.set_ylabel("$\chi(D)$")
# ax3.legend(frameon=True)

# fig.savefig(os.path.join(image_dir, "D_vs_d12.png"), bbox_inches='tight', dpi=300)
# fig2.savefig(os.path.join(image_dir, "D_vs_vwind.png"), bbox_inches='tight', dpi=300)
# fig3.savefig(os.path.join(image_dir, "D_vs_chi.png"), bbox_inches='tight', dpi=300)