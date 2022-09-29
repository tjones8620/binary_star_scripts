import numpy as np
import astropy.constants as const
import astropy.units as u
import argparse
from tabulate import tabulate

def binary_perast_vel(e, P, m1, m2):
    """
    Function that computes the y component of the velocity of
    each star in a binary system at periastron (in cgs units).

    - e : eccentricity of the binary system
    - P : period of the system in days
    - m1, m2 : masses of stars in units of M_sun
    """

    M_sun = const.M_sun.cgs # M_sun in cgs units
    G = const.G.cgs # Gravitational constant in cgs units
    
    P = P * 24 * 60 * 60 * u.s.cgs # Converting period from days to seconds
    m1, m2 = m1*M_sun, m2*M_sun # Converting masses to grams

    #Semi-major axis from Keplers third law
    a = ((G * (m1+m2) * P **2) / (4. * np.pi**2))**(1./3.) 

    # Calculating semi major axis of each individual
    # star's orbit around the centre of mass
    a1 = a / (1 + m1/m2) 
    a2 = a - a1

    # Distance from star to centre of mass at periastron
    x1 = a1 * (1 - e)
    x2 = -a2 * (1 - e)

    # Velocity of star at periastron
    vy_1 = (2 * np.pi * a1 **2 * np.sqrt(1 - e**2)) / (P*x1)
    vy_2 = (2 * np.pi * a2 **2 * np.sqrt(1 - e**2)) / (P*x2)

    return vy_1, vy_2, x1, x2

class ArgumentInputs:
    def __init__(self):
        parser = argparse.ArgumentParser(description=("Calculate the velocities and positions of masses in a stable binary orbit at periastron, with respect to the centre of mass"))
        parser.add_argument('period', help="Period of binary system in days", type=float)
        parser.add_argument('eccentricity', help='Eccentricity of the system', type=float)
        parser.add_argument('m1', help='Mass of star 1 (M_sun)', type=float)
        parser.add_argument('m2', help='Mass of star 2 (M_sun)', type=float)
        parser.add_argument('-du', '--dist_units', type=str, choices=["cgs", "au", "pc", "si"], default="cgs")
        parser.add_argument('-vu', '--vel_units', type=str, choices=["cgs", "si"], default="cgs")
        parser.add_argument('-pu', '--period_units', type=str, choices=["d", "yr"], default="d")
        parser.add_argument('-mu', '--mass_units', type=str, choices=["M_sun", "cgs", "si"], default="M_sun")
        args=parser.parse_args()

        self.Eccentricity = args.eccentricity
        self.Period = args.period
        self.Mass_1 = args.m1
        self.Mass_2= args.m2

        self.distance_units=args.dist_units
        self.velocity_units = args.vel_units
        self.period_units = args.period_units
        self.mass_units = args.mass_units

        # Changes between selected period units
        if self.period_units=="yr":
            self.Period = self.Period*(365.25)*(u.si.d)
        elif self.period_units=="d":
            self.Period = self.Period*(u.si.d)

        # Changes between inputed mass units
        if self.mass_units=="cgs":
            self.Mass_1 = ((self.Mass_1*u.g)/(const.M_sun.to(u.g)))*u.astrophys.M_sun
            self.Mass_2 = ((self.Mass_2*u.g)/(const.M_sun.to(u.g)))*u.astrophys.M_sun
        elif self.mass_units=="si":
            self.Mass_1 = ((self.Mass_1*u.kg)/(const.M_sun))*u.astrophys.M_sun
            self.Mass_2 = ((self.Mass_2*u.kg)/(const.M_sun))*u.astrophys.M_sun
        elif self.mass_units=="M_sun":
            self.Mass_1 = self.Mass_1*u.astrophys.M_sun
            self.Mass_2 = self.Mass_2*u.astrophys.M_sun
    


if __name__ == "__main__":

    inputs=ArgumentInputs()

    v1, v2, x1, x2 = binary_perast_vel(inputs.Eccentricity, inputs.Period.value, inputs.Mass_1.value, inputs.Mass_2.value)

    if inputs.velocity_units=="si":
        v1, v2 = v1.to(u.m/u.s), v2.to(u.m/u.s)
        
    if inputs.distance_units=="au":
        x1, x2 = x1.to(u.astrophys.au), x2.to(u.astrophys.au)
    elif inputs.distance_units=="pc":
        x1, x2 = x1.to(u.astrophys.pc), x2.to(u.astrophys.pc)
    elif inputs.distance_units=="si":
        x1, x2 = x1.to(u.m), x2.to(u.m)

    table = tabulate([[f"Velocity ({v1.unit})", f"{v1.value:e}", f"{v2.value:e}"],
                    [f"Distance ({x1.unit})", f"{x1.value:e}", f"{x2.value:e}"],
                    [f"Mass (M_sun)", inputs.Mass_1.value, inputs.Mass_2.value],
                    ["e", inputs.Eccentricity],
                    [f"Period ({inputs.Period.unit})", inputs.Period.value]], 
                    headers=['Quantity', 'Star 1', 'Star 2'], tablefmt='fancy_grid',
                    numalign="left"
                    )
    print(table)
    
    
