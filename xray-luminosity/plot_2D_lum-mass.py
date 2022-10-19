# Author: Sam Green, Created: 14-11-19
# Created to work with the BubbleNebula project. Plots the 3D luminosity data.

import matplotlib.pyplot as plt
import math
import numpy as np
import argparse
from astropy import units as u

parser = argparse.ArgumentParser()
parser.add_argument("sim", type=str)
args = parser.parse_args()
sim = args.sim

xray_file = "xray-lum-"+args.sim+".txt"
mass_file = "wind-mass-"+args.sim+".txt"
img_file  = "lum_"+args.sim+".png"
print(xray_file, mass_file, img_file)

pi = math.pi
m_p = 1.6726219e-24  # g
k = 1.38064852e-16  # erg/K
dist_bn = 3.086e21*1.5 # 1.5 kpc
dist_zeta1 = 3.456e+20 # cm,  Zeta @ 112pc
dist_zeta2 = 4.166e+20 # cm,  Zeta @ 135pc
mu = 0.61  # mH


# Read a line of numbers out of a text file:
data = np.loadtxt(xray_file)
sim_time = np.array(data[:,0])
sim_time = sim_time / 3.156e13
l001 = np.array(data[:,1])
l002 = np.array(data[:,2])
l003 = np.array(data[:,3])
l005 = np.array(data[:,4])
l010 = np.array(data[:,5])
l020 = np.array(data[:,6])
l050 = np.array(data[:,7])
l100 = np.array(data[:,8])

data = np.loadtxt(mass_file)
m1e4 = np.array(data[:,1])
m3e4 = np.array(data[:,2])
m1e5 = np.array(data[:,3])
m3e5 = np.array(data[:,4])
m1e6 = np.array(data[:,5])
m3e6 = np.array(data[:,6])
m1e7 = np.array(data[:,7])
m1e9 = np.array(data[:,8])


ymax = 1e+33
ymin = 1e+27


# plt.rc('font', **{'size': 14})
fig, ax1 = plt.subplots()

plt.title("Unabsorbed Flux vs Time")

ax1.set_xlabel('Time (Myr)')
ax1.set_ylabel('X-ray Luminosity (erg s$^{-1}$)')
ax1.set_yscale('log')
ax1.set_ylim(ymin, ymax)
ax1.grid(True)
ax1.plot(sim_time, (l001 - l003), "c-", label='L(0.1 - 0.3 keV)')
ax1.plot(sim_time, (l003 - l010), "m-", label='L(0.3 - 1 keV)')
ax1.plot(sim_time, (l010 - l020), "y-", label='L(1 - 2 keV)')
ax1.plot(sim_time, (l020 - l100), "k-", label='L(2 - 10 keV)')

ax2 = ax1.twinx()  # create a second axes that shares the same x-axis
ax2.plot(sim_time, m1e4, color='cyan', linestyle='--', label='$M(T<10^4$ K$)$')
ax2.plot(sim_time, m1e5, color='blue', linestyle='--', label='$M(T<10^5$ K$)$')
ax2.plot(sim_time, m1e6, color='red', linestyle='--', label='$M(T<10^6$ K$)$')
ax2.plot(sim_time, m1e7, color='green', linestyle='--', label='$M(T<10^7$ K$)$')
ax2.plot(sim_time, m1e9, color='black', linestyle='--', label='$M(T<10^9$ K$)$')
#ax2.plot(sim_time, l02, color='purple', linestyle='solid', label='L(E>0.2keV)')
#ax2.plot(sim_time, l05, color='red', linestyle='solid', label='L(E>0.5keV)')
#ax2.plot(sim_time, l1, color='yellow', linestyle='solid', label='L(E>1keV)')
#ax2.plot(sim_time, l2, color='cyan', linestyle='solid', label='L(E>2keV)')
#ax2.plot(sim_time, l5, color='lightpink', linestyle='solid', label='L(E>5keV)')
#ax2.plot(sim_time, l10, color='peru', linestyle='solid', label='L(E>10keV)')

ax2.set_ylabel('Mass of wind material ($g$)')
ax2.set_yscale('log')
ax1.legend(loc='upper left', ncol=1, framealpha=1)
ax2.legend(loc='upper right', ncol=1, framealpha=1)
# ax2.set_ylim(bottom=15e27)
ax2.set_ylim(ymin, ymax)
ax2.set_xlim((0, max(sim_time)))

# Save the figures as .png files:
plt.savefig(img_file,
            bbox_inches='tight', dpi=300)
