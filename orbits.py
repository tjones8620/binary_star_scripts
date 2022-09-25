#!/usr/bin/python
# -*- coding: iso-8859-15 -*-

import numpy as np
import astropy.constants as apc
from astropy import units as u

# constants
msun = apc.M_sun.cgs.value
rsun = apc.R_sun.cgs.value
grav = apc.G.cgs.value
day  = (1.0 * u.day).to(u.s).value

#msun = 1.989e33
#rsun = 6.96e10
#grav = 6.6743e-8
#day = 86400.0
AU = 1.495e13


# gamma2-vel
#m1 = 50.0
#m2 = 20.0
#P  = 78.5
#e  = 0.4

# wr22
m1 = 72.0
m2 = 25.7
P = 80.325
e = 0.559

# # wr140
# m1 = 29.27
# m2 = 10.31
# P = 2895
# e = 0.8993

# test
#m1 = 10 
#m2 = 10
#P  = 10
#e  = 0.5


# Kepler's Law
def semimajor(f_m1, f_m2, f_P):
  return (grav * (f_m1 + f_m2) * f_P**2 / (4.*np.pi**2))**(1./3.)

m1 = m1*msun
m2 = m2*msun
P  = P *day

print(m1,m2,P)

a = semimajor(m1,m2,P)
print("a = ","{:e}".format(a),"cm, or","{:e}".format(a/rsun),"Rsun, or","{:e}".format(a/AU),"AU")

# in centre of mass frame, each ellipse has a different size that depends on the mass ratio,
# from a1/a2 = m2/m1
a1 = a / (1.0 + m1/m2)
a2 = a - a1
b1 = np.sqrt(1.0-e**2) * a1
b2 = np.sqrt(1.0-e**2) * a2

print("a1=",a1/rsun,", a2=",a2/rsun,"R_sun")
print("a1=","{:e}".format(a1),", a2=","{:e}".format(a2),"cm")
print("b1=","{:e}".format(b1),", b2=","{:e}".format(b2),"cm")

# periastron distance from the centre of mass:
peri1 = a1 * (1-e)
peri2 = a2 * (1-e)
print("peri1=","{:e}".format(peri1),", peri2=","{:e}".format(peri2),"cm")

mred = m1*m2/(m1+m2)
grav_par = grav * (m1+m2)

# Use Kepler's laws to get the velocity at periastron:

v1 = a1**2 * np.sqrt(grav_par * (1.0-e**2)/ a**3) / peri1
#v1 = 2.0 * np.pi * a1**2 * np.sqrt(1.0-e**2)/(peri1 * P)
print("v1", "{:e}".format(v1))
v2 = a2**2 * np.sqrt(grav_par * (1.0-e**2.0) / a**3) / peri2
#v2 = 2.0 * np.pi * a2**2 * np.sqrt(1.0-e**2)/(peri2 * P)
print("v2", "{:e}".format(v2))

quit()



