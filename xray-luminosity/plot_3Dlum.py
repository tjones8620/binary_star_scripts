# Author: Sam Green, Created: 14-11-19
# Created to work with the BubbleNebula project. Plots the 3D luminosity data.

import matplotlib.pyplot as plt
import math
import numpy as np

from astropy import units as u

pi = math.pi
m_p = 1.6726219e-24  # g
k = 1.38064852e-16  # erg/K
#dist_bn = 7.4056e+21  # cm BN
dist_bn = 1.110843929e+22 # old bn, 3.5kpc
dist_zeta1 = 3.456e+20 # cm,  Zeta @ 112pc
dist_zeta2 = 4.166e+20 # cm,  Zeta @ 135pc

mu = 0.61  # mH

l01 = []
l02 = []
l03 = []
l05 = []
l1 = []
l2 = []
l5 = []
l10 = []
sim_time = []

# Read a line of numbers out of a text file:
with open("/home/visitor_ap4/code/project/scripts/xray-luminosity/xray-lum-wr140-hydro-n128.txt") as x:
    for line in x:
        data = line.split()
        l01.append(float(data[1]))
        l02.append(float(data[2]))
        l03.append(float(data[3]))
        l05.append(float(data[4]))
        l1.append(float(data[5]))
        l2.append(float(data[6]))
        l5.append(float(data[7]))
        l10.append(float(data[8]))
        sim_time.append(float(data[0]))


l01 = np.array(l01)
l02 = np.array(l02)
l03 = np.array(l03)
l05 = np.array(l05)
l1 = np.array(l1)
l2 = np.array(l2)
l5 = np.array(l5)
l10 = np.array(l10)
sim_time = np.array(sim_time) / 31536000 / 1e6


# lw = np.full(len(l01), 2.8769e34)
lw = np.full(len(l01), 3.6e36)
#fw = lw / (4.0 * pi * dist ** 2)

xmax = 1e+33
xmin = 1e+28


# plt.rc('font', **{'size': 14})
fig, ax1 = plt.subplots()

plt.title("Unabsorbed Flux vs Time")
# ax1.plot(sim_time, lw, color='black', linestyle='solid', label='Lw')
ax1.set_xlabel('Time (Myr)')
ax1.set_ylabel('Unabsorbed Flux ($erg$ $cm^{-2}$$s^{-1}$)')
# ax1.set_yscale('log')
# ax1.set_ylim((xmin / (4.0 * pi * dist_bn ** 2), xmax / (4.0 * pi * dist_bn ** 2)))
# ax1.set_ylim(bottom=1e-17)
# ax1.set_xlim((0, max(sim_time)))
ax1.grid(True)

ax2 = ax1.twinx()  # create a second axes that shares the same x-axis

# ax2.plot(sim_time, (l03 - l1), color='blue', linestyle='solid', label='L(0.3 - 1keV)')
# ax2.plot(sim_time, (l1 - l2), color='red', linestyle='solid', label='L(1 - 2keV)')
# ax2.plot(sim_time, (l2 - l10), color='black', linestyle='solid', label='L(2 - 10keV)')

ax2.plot(sim_time, l01, color='blue', linestyle='solid', label='L(E>0.1keV)')
ax2.plot(sim_time, l02, color='purple', linestyle='solid', label='L(E>0.2keV)')
ax2.plot(sim_time, l05, color='red', linestyle='solid', label='L(E>0.5keV)')
ax2.plot(sim_time, l1, color='yellow', linestyle='solid', label='L(E>1keV)')
ax2.plot(sim_time, l2, color='cyan', linestyle='solid', label='L(E>2keV)')
ax2.plot(sim_time, l5, color='lightpink', linestyle='solid', label='L(E>5keV)')
ax2.plot(sim_time, l10, color='peru', linestyle='solid', label='L(E>10keV)')


# ax2.plot(sim_time, lw, color='black', linestyle='solid')

ax2.set_ylabel('X-ray Luminosity ($erg$ $s^{-1}$)')
# ax2.set_yscale('log')
ax2.legend(loc='upper left', ncol=2, framealpha=1)
# ax2.set_ylim(bottom=15e27)
# ax2.set_ylim((xmin, xmax))
# ax2.set_xlim((0, max(sim_time)))

# Save the figures as .png files:
plt.savefig("/home/visitor_ap4/code/project/scripts/xray-luminosity/xray-lum-wr140-hydro-n128.png",
            bbox_inches='tight', dpi=300)
