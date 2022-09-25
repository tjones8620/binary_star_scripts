import numpy as np
import matplotlib.pyplot as plt
# from scipy.integrate import odeint
# from scipy.integrate import solve_ivp
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D 
from matplotlib.animation import PillowWriter
from astropy import constants as const
import astropy.units as u


from calculating_velocity import binary_perast_vel

# G = const.G.value
# AU = const.au.value
# M_sun = const.M_sun.value
# daysec = 24.0*60*60

G = const.G.cgs.value
AU = const.au.cgs.value
M_sun = const.M_sun.cgs.value
daysec = 24.0*60*60

P = 2895
e = 0.8993
m1 = 29.27
m2 = 10.31 

v1, v2, x1, x2 = binary_perast_vel(e, P, m1, m2)

# v1 = v1.value
# v2 = v2.to(u.m / u.s).value
# x1 = x1.value/AU
# x2 = x2.to(u.m).value/AU

v1 = v1.value
v2 = v2.value
x1 = x1.value/AU
x2 = x2.value/AU

print(v1, v2, x1, x2)


Ma = m1 * M_sun
Mb = m2 * M_sun

gravconst = G*Ma*Mb

#WR
xa = x1*AU
ya = 0

xva = 0
yva = v1

#O star
xb = x2*AU
yb = 0

xvb = 0
yvb = v2

t = 0.0
dt = 0.1*daysec


xalist = []
yalist = []

xblist = []
yblist = []

tlist = []

while t < P * daysec:
    #Compute force
    rx = xb-xa
    ry = yb-ya

    modr3 = (rx**2 + ry**2)**1.5

    fx = -gravconst * rx/modr3
    fy = -gravconst * ry/modr3

    # Update the quantities

    xvb += fx*dt/Mb
    yvb += fy*dt/Mb
 

    xb += xvb*dt
    yb += yvb*dt


    xva += -fx*dt/Ma
    yva += -fy*dt/Ma

    xa += xva*dt
    ya += yva*dt

    tlist.append(t)
    t += dt

    xalist.append(xa)
    yalist.append(ya)

    xblist.append(xb)
    yblist.append(yb)



# fig, ax = plt.subplots()
# ax.plot(xalist, yalist)
# ax.plot(xblist, yblist)
# ax.set_aspect(aspect=1)
# plt.show()

# # # xa_plot = xalist[::80][0]

# # # print(len(xa_plot), len(xalist))


def animate(i):
    ln1.set_data([xalist[::160][i], xblist[::160][i]], [yalist[::160][i], yblist[::160][i]])
    text.set_text('Time = {:.2f} Years'.format((tlist[::160][i])/(daysec*365)))

fig, ax = plt.subplots(1,1, figsize=(8,8))
ax.set_title("Orbital Motion of WR140")
ax.plot(xalist, yalist, ls = "--", color="k")
ax.plot(xblist, yblist, ls = "--", color="k")
ln1, = plt.plot([], [], 'o', lw=3, markersize=16, color="yellow")
text = plt.text(0, 1e14, '')
ax.set_ylim(-2e14, 2e14)
ax.set_xlim(-1.5e14, 3.5e14)
ani = animation.FuncAnimation(fig, animate, frames=(len(xalist[::160])-1), interval=5)
ani.save('binary.gif',writer='pillow',fps=500)

# def animate(i):
#     ln1.set_data([xalist[::160][i]], [yalist[::160][i]])
#     text.set_text('Time = {:.2f} Years'.format((tlist[::160][i])/(daysec*365)))

# fig, ax = plt.subplots(1,1, figsize=(8,8))
# ax.set_title("Orbital Motion of WR140")
# ax.plot(xalist, yalist, ls = "--", color="k")
# # ax.plot(xblist, yblist, ls = "--", color="k")
# ln1, = plt.plot([], [], 'o', lw=3, markersize=16, color="yellow")
# text = plt.text(0, 1e14, '')
# ax.set_ylim(-2e14, 2e14)
# ax.set_xlim(-1.5e14, 3.5e14)
# ani = animation.FuncAnimation(fig, animate, frames=(len(xalist[::160])-1), interval=5)
# ani.save('binary1.gif',writer='pillow',fps=500)