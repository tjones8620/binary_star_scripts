# import argparse
# from ast import arg, parse

import astropy.units as u
import astropy.constants as const

q = 25*const.R_sun.cgs
print(q.unit)
print(f"{q:2e}")