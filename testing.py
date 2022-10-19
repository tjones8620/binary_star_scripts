import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as const

from calculating_covertex_velocity import BinarySystem

WR140 = BinarySystem(0.8993, 2895, 10.31, 29.27, "WR140")
print(WR140.vy1_per)

WR140.integrate_orbit()

