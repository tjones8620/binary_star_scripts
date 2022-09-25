import numpy as np
import astropy.constants as const
import astropy.units as u

def binary_perast_vel(e, P, m1, m2):
    """
    Function that computes the y component of the velocity of
    each star in a binary system at periastron (in cgs units).

    - e : eccentricity of the binary system
    - P : period of the system in days
    - m1, m2 : masses of stars in units of M_sun
    """

    M_sun = const.M_sun.cgs # M_sun in cgs units
    AU = const.au.cgs # 1AU in cgs units
    G = const.G.cgs # Gravitational constant in cgs units
    
    P = P * 24 * 60 * 60 * u.s.cgs # Converting period from days to seconds
    m1, m2 = m1*M_sun, m2*M_sun

    a = ((G * (m1+m2) * P **2) / (4. * np.pi**2))**(1./3.)

    a1 = a / (1 + m1/m2)
    a2 = a - a1
    x1 = a1 * (1 - e)
    x2 = -a2 * (1 - e)
    vy_1 = (2 * np.pi * a1 **2 * np.sqrt(1 - e**2)) / (P*x1)
    vy_2 = (2 * np.pi * a2 **2 * np.sqrt(1 - e**2)) / (P*x2)

    # print(f"Semimajor axis of m1 orbit: {a1:e}")
    # print(f"Semimajor axis of m2 orbit: {a2:e}")
    # print(f"Separation between m1 and m2 at periastron: {}")

    return vy_1, vy_2, x1, x2

P = 2895
e = 0.8993
m1 = 29.27
m2 = 10.31 

# # wr22
# m1 = 72.0
# m2 = 25.7
# P = 80.325
# e = 0.559

v1, v2, x1, x2 = binary_perast_vel(e, P, m1, m2)
print(f"{v1:e}", f"{v2:e}")
print(f"{x1:e}", f"{x2:e}")
exit()

