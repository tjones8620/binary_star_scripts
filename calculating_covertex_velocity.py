from dataclasses import dataclass
from typing import ClassVar
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D 
from matplotlib.animation import PillowWriter
from astropy import constants as const
import astropy.units as u
from scipy.integrate import odeint, solve_ivp

from calculating_velocity import binary_perast_vel


class BinarySystem():

    G = const.G.cgs
    M_sun = const.M_sun.cgs
    AU = const.au.cgs

    def __init__(self, e, P, m1, m2):
        self.Ecc = e
        self.Period = P
        self.m1 = m1
        self.m2 = m2


    def binary_perast_vel(self):
        """
        Computes the y component of the velocity of
        each star in a binary system at periastron 
        and the distance from each star to the centre
        of mass of the system (cgs units). 

        - e : eccentricity of the binary system
        - P : period of the system in days
        - m1, m2 : masses of stars in units of M_sun

        
        """

        
        P = self.Period * 24 * 60 * 60 * u.s.cgs # Converting period from days to seconds
        m1, m2 = self.m1 * BinarySystem.M_sun, self.m2 * BinarySystem.M_sun # Converting masses to grams

        #Semi-major axis from Keplers third law
        a = ((BinarySystem.G * (m1+m2) * P **2) / (4. * np.pi**2))**(1./3.) 

        # Calculating semi major axis of each individual
        # star's orbit around the centre of mass
        a1 = a / (1 + m1/m2) 
        a2 = a - a1

        # Distance from star to centre of mass at periastron
        self.x1 = a1 * (1 - self.Ecc)
        self.x2 = -a2 * (1 - self.Ecc)

        # Velocity of star at periastron
        self.vy_1 = (2 * np.pi * a1 **2 * np.sqrt(1 - self.Ecc**2)) / (P*self.x1)
        self.vy_2 = (2 * np.pi * a2 **2 * np.sqrt(1 - self.Ecc**2)) / (P*self.x2)

        return self.vy_1, self.vy_2, self.x1, self.x2    

    
    def binary_covert_vel(self):
    
        BinarySystem.binary_perast_vel(self)
        
        m1 = self.m1 * BinarySystem.M_sun.value
        m2 = self.m2 * BinarySystem.M_sun.value

        x1_0 = self.x1.value
        y1_0 = 0
        x2_0 = self.x2.value
        y2_0 = 0
        vx1_0 = 0
        vy1_0 = self.vy_1.value
        vx2_0 = 0
        vy2_0 = self.vy_2.value

        def dSdt(S, t):
            x1, y1, x2, y2, vx1, vy1, vx2, vy2 = S
            r12 = np.sqrt((x2-x1)**2 + (y2-y1)**2)

            return [
                vx1,
                vy1,
                vx2,
                vy2,
                BinarySystem.G.value*m2/r12**3 *(x2-x1),
                BinarySystem.G.value*m2/r12**3 *(y2-y1),
                BinarySystem.G.value*m1/r12**3 *(x1-x2),
                BinarySystem.G.value*m1/r12**3 *(y1-y2)
            ]       


        t = np.linspace(0, self.Period*24*60*60, 10000)

        sol = odeint(dSdt , y0=[x1_0, y1_0, x2_0, y2_0, vx1_0, vy1_0, vx2_0, vy2_0], t=t)
        sol = sol.T

        x1, y1, x2, y2, vx1, vy1, vx2, vy2 = sol                                                                                                                                                                                                                                                                                                                            

        zero_crossings = np.where(np.diff(np.sign(vy1)))[0]
        zero_crossings

        self.vx1_covert = vx1[zero_crossings[0]]*u.cm/u.s
        self.vx2_covert = vx2[zero_crossings[0]]*u.cm/u.s

        self.x1_covert = x1[zero_crossings[0]]*u.cm
        self.x2_covert = x2[zero_crossings[0]]*u.cm

        return self.vx1_covert, self.vx2_covert, self.x1_covert, self.x2_covert

        
WR140 = BinarySystem(0.8993, 2895, 10.31, 29.27)
print(WR140.binary_perast_vel())
WR140.binary_covert_vel()

print(WR140.binary_covert_vel())
# # print(WR140.binary_perast_vel())

# print(WR140.M_sun)

