import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D 
from matplotlib.animation import PillowWriter
from astropy import constants as const
import astropy.units as u
from scipy.integrate import odeint, solve_ivp
import math
import matplotlib
# plt.style.use("science")


class BinarySystem():
    G = const.G.cgs
    M_sun = const.M_sun.cgs
    AU = const.au.cgs

    def __init__(self, e, P, m1, m2, name):
        self.Ecc = e
        self.Period = P
        self.m1 = m1
        self.m2 = m2
        self.name = name


    def get_periast_v_x(self):
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
        self.x1_per = a1 * (1 - self.Ecc)
        self.x2_per = -a2 * (1 - self.Ecc)

        # Velocity of star at periastron
        self.vy1_per = (2 * np.pi * a1 **2 * np.sqrt(1 - self.Ecc**2)) / (P*self.x1_per)
        self.vy2_per = (2 * np.pi * a2 **2 * np.sqrt(1 - self.Ecc**2)) / (P*self.x2_per)

        return self.vy1_per, self.vy2_per, self.x1_per, self.x2_per    

    
    def integrate_orbit(self, k=1):
        """
        Solves for the periastron velocity and position using the get_periast_v_x() method
        and uses this as the initial conditions in order to integrates Newtons equations to
        calculate the x,y positions and velocities for each star at closely spaced discrete 
        points around the stars orbit. 

        - k : Number of periods to plot over. Default set to 1.

        
        """
    
        BinarySystem.get_periast_v_x(self)
        
        m1 = self.m1 * BinarySystem.M_sun.value
        m2 = self.m2 * BinarySystem.M_sun.value

        x1_0 = self.x1_per.value
        y1_0 = 0
        x2_0 = self.x2_per.value
        y2_0 = 0
        vx1_0 = 0
        vy1_0 = self.vy1_per.value
        vx2_0 = 0
        vy2_0 = self.vy2_per.value

        def dSdt(S, t):
            x1, y1, x2, y2, vx1, vy1, vx2, vy2 = S
            r12 = np.sqrt((x2-x1)**2 + (y2-y1)**2)

            return [vx1, vy1, vx2, vy2,
                    BinarySystem.G.value*m2/r12**3 *(x2-x1),
                    BinarySystem.G.value*m2/r12**3 *(y2-y1),
                    BinarySystem.G.value*m1/r12**3 *(x1-x2),
                    BinarySystem.G.value*m1/r12**3 *(y1-y2)]       

        self.t_array = np.linspace(0, k*self.Period*24*60*60, 1000)

        sol = odeint(dSdt , y0=[x1_0, y1_0, x2_0, y2_0, vx1_0, vy1_0, vx2_0, vy2_0], t=self.t_array)
        sol = sol.T

        x1, y1, x2, y2, vx1, vy1, vx2, vy2 = sol

        return sol

    def plot_binary_orbit(self, make_animation=False, length=10, set_fps=30, int_time=50):

        sol = BinarySystem.integrate_orbit(self, k=1)
        x1, y1, x2, y2, vx1, vy1, vx2, vy2 = sol
        fig, ax = plt.subplots()
        ax.plot(x1,y1)
        ax.plot(x2,y2)
        ax.set_title(f"Orbit of {self.name}")
        ax.set_aspect(aspect=1)
        # plt.show()
        #### need to make path to save image ####

        if make_animation==True:
 
            num_frames = int(2.5*length * (1/set_fps + int_time*1e-3)**-1)
            spacing = int(len(x1)/num_frames)
            num_frames = math.ceil(len(x1)/spacing)

            if np.abs(y1).max() >= np.abs(y2).max():
                ylim = np.abs(y1).max()
            else:
                ylim = np.abs(y2).max()
            
            t_array_yr = self.t_array/(365*24*60*60)

            def animate(i):
                ln1.set_data([x1[::spacing][i], x2[::spacing][i]], [y1[::spacing][i], y2[::spacing][i]])
                text.set_text('Time = {:.1f} Years'.format(t_array_yr[::spacing][i]))
            fig, ax = plt.subplots()
            ax.set_title(f"Orbital Motion of {self.name}")
            ax.plot(x1, y1, ls = "--", color="w", lw=1)
            ax.plot(x2, y2, ls = "--", color="w", lw=1)
            ln1, = plt.plot([], [], 'o', lw=3, markersize=10, color="orange")
            text = plt.text(0.6, 0.85, '', color='white', transform=ax.transAxes, fontsize = 10)
            # ax.text(0.7, 0.7, '',transform=ax.transAxes)
            ax.set_xlim(-1.2*np.abs(x1).max(), 1.2*np.abs(x2).max())
            ax.set_ylim(-1.2*((np.abs(x1).max() + np.abs(x2).max())/2), 1.2*((np.abs(x1).max() + np.abs(x2).max())/2))
            ax.set_aspect(aspect=1)
            ax.set_facecolor('k')
            ax.set_xlabel("x [cm]", fontsize=8)
            ax.set_ylabel("y [cm]", fontsize=8)
            ani = animation.FuncAnimation(fig, animate, frames=num_frames, interval=int_time)
            ani.save('WR140_binary.gif',writer='pillow',fps=set_fps, dpi=200)

    def binary_covert_vxy(self): 
        """
        Calculates the velocity of each star at the covertex of their elliptical orbit, 
        as well as the x and y positions of the stars.
        At the covertex the Cartesian y-velocity is 0 
        """  

        #### Fix velocity and position signs ***

        sol = BinarySystem.integrate_orbit(self, k=1)                                                                                                                                                                                                                                                                                                                         

        x1, y1, x2, y2, vx1, vy1, vx2, vy2 = sol

        zero_crossings = np.where(np.diff(np.sign(vy1)))[0]
        zero_crossings

        self.vx1_covert = vx1[zero_crossings[0]]*u.cm/u.s
        self.vx2_covert = vx2[zero_crossings[0]]*u.cm/u.s

        self.x1_covert = x1[zero_crossings[0]]*u.cm
        self.x2_covert = x2[zero_crossings[0]]*u.cm

        self.y1_covert = y1[zero_crossings[0]]*u.cm
        self.y2_covert = y2[zero_crossings[0]]*u.cm

        return self.vx1_covert, self.vx2_covert, self.x1_covert, self.x2_covert, self.y1_covert, self.y2_covert




WR140 = BinarySystem(0.8993, 2895, 10.31, 29.27, "WR140")
WR140.plot_binary_orbit(make_animation=True)
# print(WR140.binary_covert_vel())
# print(WR140.binary_perast_vel())
# WR140.binary_covert_vel()

# print(WR140.binary_covert_vel())
# # # print(WR140.binary_perast_vel())

# # print(WR140.M_sun)

