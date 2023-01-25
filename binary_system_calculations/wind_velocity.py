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

D_array = np.linspace(D, 3*D, 1000) 


######### Accelerating winds case #########
# Calculating the distance from the stars to contact 
# at each distance for accelerating winds case
r_WR_array = 1/(1+phi) * D_array 
r_O_array = D_array - r_WR_array

# Calculating the wind velocity at the contact surface
# for each star at each distance for accelerating winds case
vels_O = velocity(r_O_array, vO, R_O)
vels_WR = velocity(r_WR_array, vWR, R_WR)

# Calculating chi for each star at each distance for accelerating winds case
chis_O = chi(vels_O, r_O_array, MdotO)
chis_WR = chi(vels_WR, r_WR_array, MdotWR)

print("Minimum chi values for accelerating winds case:")
print(f"chi min O = {min(chis_O)} ")
print(f"chi min WR = {min(chis_WR)} ")

######### Constant winds case #########

# Calculating the distance from the stars to contact
# at each distance for constant winds case
r_WR_array_const = 1/(1+phi) * D_array
r_O_array_const = D_array - r_WR_array_const

# Calculating the wind velocity at the contact surface
# for each star at each distance for constant winds case
vels_O_const = vO * np.ones(len(D_array))
vels_WR_const = vWR * np.ones(len(D_array))

# Calculating chi for each star at each distance for constant winds case
chis_O_const = chi(vels_O_const, r_O_array_const, MdotO)
chis_WR_const = chi(vels_WR_const, r_WR_array_const, MdotWR)

x = np.linspace(0, 10*D, 1000)

print("Minimum chi values for constant winds case:")
print(f"chi min O = {min(chis_O_const)} ")
print(f"chi min WR = {min(chis_WR_const)} ")


######### Plotting #########
fig, ax = plt.subplots(3,1, figsize=(5,8), sharex=True)
ax[0].plot(D_array/D, r_O_array/D, label="O star")
ax[0].plot(D_array/D, r_WR_array/D, label="WR star")
ax[0].set_ylabel("$d_{12}/d_{per}$")
# ax[0].legend()

ax[1].plot(D_array/D, vels_O/1000, label="O star", color="C0")
ax[1].plot(D_array/D, vels_WR/1000, label="WR star", color="C1")
ax[1].set_ylim(1, 3.3)
ax[1].set_ylabel("$V_{w}$ $(10^3$ $km$ $s^{-1}$)")

ax[2].plot(D_array/D, chis_O, label="O (AW)", color="C0")
ax[2].plot(D_array/D, chis_WR, label="WR (AW)", color="C1")
ax[2].plot(D_array/D, chis_O_const, label="O (CW)", linestyle="--", color="C12")
ax[2].plot(D_array/D, chis_WR_const, label="WR (CW)", linestyle=":", color="C12")
ax[2].legend(fontsize=8, ncol=2)

ax[2].set_xlabel("$D$ $(d_{per})$")
ax[2].set_ylabel("$\chi(D)$")

fig.subplots_adjust(hspace=0.05)

fig.savefig(os.path.join(image_dir, "wind_velocity.png"), bbox_inches='tight', dpi=500)