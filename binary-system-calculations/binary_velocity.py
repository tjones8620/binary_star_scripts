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
import argparse
import os
from tabulate import tabulate
# plt.style.use("science")

class ArgumentInputs:
    def __init__(self):
        parser = argparse.ArgumentParser(description=("Calculate the velocities and positions of masses in a stable binary orbit at periastron, with respect to the centre of mass"))
        parser.add_argument('period', help="Period of binary system in days", type=float)
        parser.add_argument('eccentricity', help='Eccentricity of the system', type=float)
        parser.add_argument('m1', help='Mass of star 1 (M_sun)', type=float)
        parser.add_argument('m2', help='Mass of star 2 (M_sun)', type=float)
        parser.add_argument('name', help='Name of the system', type=str)
        parser.add_argument('--plot', help='Plot the orbit', action='store_true')
        parser.add_argument('--covertex', help='Return the co-vertex (x,y) and Vx', action='store_true')
        parser.add_argument('--arb_phase', help='Arbitrary phase of the system', type=float)
        parser.add_argument('--animate', help='Animate the orbit', action='store_true')
        args=parser.parse_args()

        self.eccentricity = args.eccentricity
        self.period = args.period
        self.mass1 = args.m1
        self.mass2= args.m2
        self.name= args.name

        self.covertex = args.covertex
        self.arb_phase = args.arb_phase
        self.plot = args.plot
        self.ani = args.animate

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

        periastron_v_x = self.get_periast_v_x(e, P, m1, m2)
        self.vy1_per, self.vy2_per, self.x1_per, self.x2_per = periastron_v_x

        init_cond={'m1':m1, 'm2':m2, 
                'x1':self.x1_per.value, 'y1':0, 
                'x2':self.x2_per.value, 'y2':0, 
                'vx1':0, 'vy1':self.vy1_per.value, 
                'vx2':0, 'vy2':self.vy2_per.value}

        self.orb_sol = self.integrate_orbit(k=1, **init_cond)


    @staticmethod
    def get_periast_v_x(e, P, m1, m2):
        """
        Computes the y component of the velocity of
        each star in a binary system at periastron 
        and the distance from each star to the centre
        of mass of the system (cgs units). 

        - e : eccentricity of the binary system
        - P : period of the system in days
        - m1, m2 : masses of stars in units of M_sun

        """
        G = const.G.cgs
        M_sun = const.M_sun.cgs

        Period = P * 24 * 60 * 60 * u.s.cgs # Converting period from days to seconds
        m_1, m_2 = m1 * M_sun, m2 * M_sun # Converting masses to grams

        #Semi-major axis from Keplers third law
        a = ((BinarySystem.G * (m_1+m_2) * Period **2) / (4. * np.pi**2))**(1./3.) 

        # Calculating semi major axis of each individual
        # star's orbit around the centre of mass
        a1 = a / (1 + m_1/m_2) 
        a2 = a - a1

        # Distance from star to centre of mass at periastron
        x1_per = a1 * (1 - e)
        x2_per = -a2 * (1 - e)

        # Velocity of star at periastron
        vy1_per = (2 * np.pi * a1 **2 * np.sqrt(1 - e**2)) / (Period*x1_per)
        vy2_per = (2 * np.pi * a2 **2 * np.sqrt(1 - e**2)) / (Period*x2_per)

        return  vy1_per, vy2_per, x1_per, x2_per

    ########################################################################################
    
    def integrate_orbit(self, k=1, **kwargs):
        """
        Solves for the periastron velocity and position using the get_periast_v_x() method
        and uses this as the initial conditions in order to integrates Newtons equations to
        calculate the x,y positions and velocities for each star at closely spaced discrete 
        points around the stars orbit. 

        - k : Number of periods to plot over. Default set to 1.
        - **kwargs : keyword arguments for the get_periast_v_x() method
        """

        G = const.G.cgs
        M_sun = const.M_sun.cgs
        
        m1 = kwargs['m1'] * M_sun.value
        m2 = kwargs['m2'] * M_sun.value

        x1_0 = kwargs['x1']
        y1_0 = kwargs['y1']
        x2_0 = kwargs['x2']
        y2_0 = kwargs['y2']

        vx1_0 = kwargs['vx1']
        vy1_0 = kwargs['vy1']
        vx2_0 = kwargs['vx2']
        vy2_0 = kwargs['vy2']

        def dSdt(S, t):
            x1, y1, x2, y2, vx1, vy1, vx2, vy2 = S
            r12 = np.sqrt((x2-x1)**2 + (y2-y1)**2)

            return [vx1, vy1, vx2, vy2,
                    G.value*m2/r12**3 *(x2-x1),
                    G.value*m2/r12**3 *(y2-y1),
                    G.value*m1/r12**3 *(x1-x2),
                    G.value*m1/r12**3 *(y1-y2)]       

        self.t_array = np.linspace(0, k*self.Period*24*60*60, 1000)

        sol = odeint(dSdt , y0=[x1_0, y1_0, x2_0, y2_0, vx1_0, vy1_0, vx2_0, vy2_0], t=self.t_array)
        sol = sol.T

        x1_array, y1_array, x2_array, y2_array, vx1_array, vy1_array, vx2_array, vy2_array = sol
        return sol
    
    ########################################################################################

    def plot_binary_orbit(self, sol, make_animation=False, length=10, set_fps=30, int_time=50):

        x1, y1, x2, y2, vx1, vy1, vx2, vy2 = sol
        fig, ax = plt.subplots()
        ax.plot(x1,y1, label='Star 1')
        ax.plot(x2,y2, label='Star 2')
        ax.set_title(f"Orbit of {self.name}")
        ax.set_aspect(aspect=1)
        ax.set_xlabel("x (cm)")
        ax.set_ylabel("y (cm)")
        ax.grid()
        plt.savefig(os.path.join(os.path.abspath(os.path.dirname(__file__)),f"orbit_{self.name}.png"))

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
            ani.save(f'{self.name}_orbit.gif',writer='pillow',fps=set_fps, dpi=200)

    ########################################################################################

    def binary_covert_vxy(self): 
        """
        Calculates the velocity of each star at the covertex of their elliptical orbit, 
        as well as the x and y positions of the stars.
        At the covertex the Cartesian y-velocity is 0 
        """  

        #### Fix velocity and position signs ***                                                                                                                                                                                                                                                                                                                               

        x1, y1, x2, y2, vx1, vy1, vx2, vy2 = self.orb_sol

        zero_crossings = np.where(np.diff(np.sign(vy1)))[0]
        zero_crossings

        self.vx1_covert = -1*vx1[zero_crossings[0]]*u.cm/u.s
        self.vx2_covert = -1*vx2[zero_crossings[0]]*u.cm/u.s

        self.x1_covert = x1[zero_crossings[0]]*u.cm
        self.x2_covert = x2[zero_crossings[0]]*u.cm

        self.y1_covert = -1*y1[zero_crossings[0]]*u.cm
        self.y2_covert = -1*y2[zero_crossings[0]]*u.cm

        return self.vx1_covert, self.vx2_covert, self.x1_covert, self.x2_covert, self.y1_covert, self.y2_covert
    
    ########################################################################################

    def arbitrary_orbital_phase_vxy(self, orbital_solution, phase):
        """
        Calculates the velocity of each star at an arbitrary orbital phase, 
        as well as the x and y positions of the stars.

        0 < phase < 1
        """                                                                                                                                                                                                                                                                                                                        

        x1, y1, x2, y2, vx1, vy1, vx2, vy2 = orbital_solution

        vx1_arb = vx1[int(phase*len(x1))]*u.cm/u.s
        vx2_arb = vx2[int(phase*len(x1))]*u.cm/u.s
        vy1_arb = vy1[int(phase*len(x1))]*u.cm/u.s
        vy2_arb = vy2[int(phase*len(x1))]*u.cm/u.s
        
        x1_arb = x1[int(phase*len(x1))]*u.cm
        x2_arb = x2[int(phase*len(x1))]*u.cm
        y1_arb = y1[int(phase*len(x1))]*u.cm
        y2_arb = y2[int(phase*len(x1))]*u.cm


        t_array_yr = self.t_array/(365*24*60*60)

        fig, ax = plt.subplots()
        text = plt.text(0.6, 0.85, '', color='white', transform=ax.transAxes, fontsize = 10)
        text.set_text('Time = {:.1f} Years'.format(t_array_yr[int(phase*len(x1))]))
        ax.plot(x1, y1, ls = "--", color="w", lw=1)
        ax.plot(x2, y2, ls = "--", color="w", lw=1)
        ax.plot(x1_arb.value, y1_arb.value, 'o', color="orange", markersize=10)
        ax.plot(x2_arb.value, y2_arb.value, 'o', color="orange", markersize=10)
        ax.set_xlim(-1.2*np.abs(x1).max(), 1.2*np.abs(x2).max())
        ax.set_ylim(-1.2*((np.abs(x1).max() + np.abs(x2).max())/2), 1.2*((np.abs(x1).max() + np.abs(x2).max())/2))
        ax.set_aspect(aspect=1)
        ax.set_facecolor('k')
        ax.set_xlabel("x [cm]", fontsize=8)
        ax.set_ylabel("y [cm]", fontsize=8)
        ax.set_title(f"Orbital Motion of {self.name}")
        plt.savefig(f"{self.name}_arbitrary_orbital_phase.png", dpi=200)

        return vx1_arb, vx2_arb, vy1_arb, vy2_arb, x1_arb, x2_arb, y1_arb, y2_arb



def main():
    cmdargs = ArgumentInputs()
    binary = BinarySystem(0.8993, 2895, 10.31, 29.27, "WR140")
    binary = BinarySystem(m1=cmdargs.mass1, m2=cmdargs.mass2, P=cmdargs.period, e=cmdargs.eccentricity, name=cmdargs.name)
    

    headers = ['Quantity', 'Star 1', 'Star 2']

    masses = ['Mass [Msun]', f"{binary.m1}", f"{binary.m2}"]
    period = ['Period [days]', f"{binary.Period}",]
    per_velocities = ['Periastron Vy [cm/s]', f"{binary.vy1_per.value:e}", f"{binary.vy2_per.value:e}"]
    per_positions = ['Periastron X [cm]', f"{binary.x1_per.value:e}", f"{binary.x2_per.value:e}"]

    table_items = [masses, period, per_velocities, per_positions]
    
    if cmdargs.plot:
        if not cmdargs.ani:
            binary.plot_binary_orbit(binary.orb_sol, make_animation=False)
        else:
            binary.plot_binary_orbit(binary.orb_sol, make_animation=True)

    if cmdargs.covertex:
        binary.binary_covert_vxy()
        covert_vx = ['Covertex Vx [cm/s]', f"{binary.vx1_covert.value:e}", f"{binary.vx2_covert.value:e}"]
        covert_x = ['Covertex X [cm]', f"{binary.x1_covert.value:e}", f"{binary.x2_covert.value:e}"]
        covert_y = ['Covertex Y [cm]', f"{binary.y1_covert.value:e}", f"{binary.y2_covert.value:e}"]

        table_items.extend([covert_vx, covert_x, covert_y])

    if cmdargs.arb_phase:
        arb_vxy = binary.arbitrary_orbital_phase_vxy(binary.orb_sol, cmdargs.arb_phase)
        vx1_arb, vx2_arb, vy1_arb, vy2_arb, x1_arb, x2_arb, y1_arb, y2_arb = arb_vxy

        arb_x = ['Arbitrary X [cm]', f"{x1_arb.value:e}", f"{x2_arb.value:e}"]
        arb_y = ['Arbitrary Y [cm]', f"{y1_arb.value:e}", f"{y2_arb.value:e}"]
        arb_vx = ['Arbitrary Vx [cm/s]', f"{vx1_arb.value:e}", f"{vx2_arb.value:e}"]
        arb_vy = ['Arbitrary Vy [cm/s]', f"{vy1_arb.value:e}", f"{vy2_arb.value:e}"]

        table_items.extend([arb_x, arb_y, arb_vx, arb_vy])

    


    table = tabulate(table_items
                    ,
                    headers=headers, 
                    tablefmt='fancy_grid', 
                    numalign="left"
                    )

    print(table)


if __name__ == "__main__":
    main()

