import argparse
import astropy.constants as const
import astropy.units as u

from calculating_covertex_velocity import BinarySystem

class ArgumentInputs:
    def __init__(self):
        parser = argparse.ArgumentParser(description=("Calculate the velocities and positions of masses in a stable binary orbit at periastron, with respect to the centre of mass"))
        parser.add_argument('period', help="Period of binary system in days", type=float)
        parser.add_argument('eccentricity', help='Eccentricity of the system', type=float)
        parser.add_argument('m1', help='Mass of star 1 (M_sun)', type=float)
        parser.add_argument('m2', help='Mass of star 2 (M_sun)', type=float)
        parser.add_argument('-du', '--dist_units', help="Change distance units for output distances",type=str, choices=["cgs", "au", "pc", "si"], default="cgs")
        parser.add_argument('-vu', '--vel_units', help="Change velocity units for output velocities", type=str, choices=["cgs", "si"], default="cgs")
        parser.add_argument('-pu', '--period_units', help="Change units of time for period input",type=str, choices=["d", "yr"], default="d")
        parser.add_argument('-mu', '--mass_units', type=str, help="Change mass units for mass inputs", choices=["M_sun", "cgs", "si"], default="M_sun")
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

inputs = ArgumentInputs()

WR140 = BinarySystem(inputs.Eccentricity, inputs.Period, inputs.Mass_1, inputs.Mass_2, "WR140")
WR140.plot_binary_orbit(make_animation=False)