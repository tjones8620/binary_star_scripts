# Author: Sam Green, Created: 12-08-2020
# Analyses 3D data.

# (1) install pypion with `python3 -m pip install pypion`
# (2) install python3-silo from https://git.dias.ie/massive-stars-software/pypion.git
#     - git clone https://git.dias.ie/massive-stars-software/pypion.git
#     - cd pypion/src/silo/
#     - bash install_silo.sh

import sys
sys.path.insert(0,"~/.local/silo/lib")

from pypion.ReadData import ReadData

import numpy as np
import math
from astropy import units as u
import time

pi = math.pi
m_p = 1.6726219e-24  # g
k = 1.38064852e-16  # erg/K
mu = 0.61  # mH


class Analysis(ReadData):
###############################################################################
    def volume3D(self, xmax, xmin, ngrid):  # Calculates the volume of each cell in the image grid

        xmax = xmax
        xmin = xmin
        ngrid = ngrid

        # Calculate the size of each cell in the x, y, and z-direction:
        delta_x = (xmax[0] - xmin[0]) / ngrid[0]
        delta_y = (xmax[1] - xmin[1]) / ngrid[1]
        delta_z = (xmax[2] - xmin[2]) / ngrid[2]

        # Create a 3D array of zeros with the same dimensions as ngrid:
        v = np.zeros((ngrid[2], ngrid[1], ngrid[0]))

        # Loop through each element in the Volume array:
        for zcells in range(ngrid[2]):
            for ycells in range(ngrid[1]):
                for xcells in range(ngrid[0]):
                    v[zcells, ycells, xcells] = delta_x * delta_y * delta_z

        # delete variables:
        del delta_x
        del delta_z
        del delta_y
        del xmax
        del xmin
        del ngrid

        return v

###############################################################################
    def volume2D(self, xmax, xmin, ngrid):  # Calculates the volume of each cell in the image grid

        xmax = xmax
        xmin = xmin
        ngrid = ngrid
        # Calculate the size of each cell in the x, y, and z-direction:
        delta_z = (xmax[0] - xmin[0]) / ngrid[0]
        delta_R = (xmax[1] - xmin[1]) / ngrid[1]
        # Create a 2D array of zeros with the same dimensions as ngrid:
        v = np.zeros((ngrid[1], ngrid[0]))
        # Loop through each element in the Volume array:
        for ycells in range(ngrid[1]):
          rmin = ycells*delta_R
          rmax = (ycells+1)*delta_R
          for xcells in range(ngrid[0]):
                    v[ycells, xcells] = delta_z * np.pi * (rmax**2 - rmin**2)
        del delta_z
        del delta_R
        del xmax
        del xmin
        del ngrid
        return v

###############################################################################
    def volume1D(self, xmax, xmin, ngrid):  # Calculates the volume of each cell in the image grid

        xmax = xmax
        xmin = xmin
        ngrid = ngrid

        # Calculate the size of each cell in the x
        dr = (xmax[0] - xmin[0]) / ngrid[0]

        # Create a 1D array of zeros with the same dimensions as ngrid:
        v = np.zeros((ngrid[0]))

        # Loop through each element in the Volume array:
        for xcells in range(ngrid[0]):
          rmin = xcells*dr
          rmax = (xcells+1)*dr
          v[xcells] = 4.0*np.pi/3.0 * (rmax**3 - rmin**3)

        # delete variables:
        del dr
        del xmax
        del xmin
        del ngrid

        return v
###############################################################################

###############################################################################
    def luminosity_interp3D(self, log_t, l01, l02, l03, l05, l1, l2, l5, l10):
        # This function uses an interpolation method to take known values from
        # a .txt to calculate the x-ray luminosity of each cell

        density = self.get_3Darray("Density")['data']
        temp = self.get_3Darray("Temperature")['data']
        mask = self.get_3Darray("NG_Mask")['data']

        ngrid = self.ngrid()

        lim_max = (self.get_3Darray("Density")['max_extents'] * u.cm).value
        lim_min = (self.get_3Darray("Density")['min_extents'] * u.cm).value

        sim_time = self.get_3Darray("Density")['sim_time']

        l_01 = np.zeros(len(density))
        l_02 = np.zeros(len(density))
        l_03 = np.zeros(len(density))
        l_05 = np.zeros(len(density))
        l_1 = np.zeros(len(density))
        l_2 = np.zeros(len(density))
        l_5 = np.zeros(len(density))
        l_10 = np.zeros(len(density))

        for j in range(len(density)):
            denj = density[j]
            tempj = temp[j]
            maskj = mask[j]

            volj = self.volume3D(lim_max[j], lim_min[j], ngrid)

            li_01 = np.zeros(ngrid[0])
            li_02 = np.zeros(ngrid[0])
            li_03 = np.zeros(ngrid[0])
            li_05 = np.zeros(ngrid[0])
            li_1 = np.zeros(ngrid[0])
            li_2 = np.zeros(ngrid[0])
            li_5 = np.zeros(ngrid[0])
            li_10 = np.zeros(ngrid[0])

            for i in range(ngrid[0]):
                log_masktemp = np.log10(tempj[:, i, :]) * maskj[:, i, :]
                n_ei = np.multiply(denj[:, i, :], np.array(1 / (1.167 * m_p)))
                n_hi = np.multiply(denj[:, i, :], np.array(1 / (1.4 * m_p)))
                vol = volj[:, i, :]

                # Use an interpolation function to calculate the x-ray luminosity at different energy levels:
                li_01[i] = sum(sum(np.multiply(np.interp(log_masktemp, log_t, l01, left=0), vol * n_ei * n_hi)))  # L(E>0.1kev)
                li_02[i] = sum(sum(np.multiply(np.interp(log_masktemp, log_t, l02, left=0), vol * n_ei * n_hi)))  # L(E>0.2kev)
                li_03[i] = sum(sum(np.multiply(np.interp(log_masktemp, log_t, l03, left=0), vol * n_ei * n_hi)))  # L(E>0.3kev)
                li_05[i] = sum(sum(np.multiply(np.interp(log_masktemp, log_t, l05, left=0), vol * n_ei * n_hi)))  # L(E>0.5kev)
                li_1[i] = sum(sum(np.multiply(np.interp(log_masktemp, log_t, l1, left=0), vol * n_ei * n_hi)))  # L(E>1kev)
                li_2[i] = sum(sum(np.multiply(np.interp(log_masktemp, log_t, l2, left=0), vol * n_ei * n_hi)))  # L(E>2kev)
                li_5[i] = sum(sum(np.multiply(np.interp(log_masktemp, log_t, l5, left=0), vol * n_ei * n_hi)))  # L(E>5kev)
                li_10[i] = sum(sum(np.multiply(np.interp(log_masktemp, log_t, l10, left=0), vol * n_ei * n_hi)))  # L(E>10kev)

                del log_masktemp, n_ei, vol

            l_01[j] = sum(li_01)
            l_02[j] = sum(li_02)
            l_03[j] = sum(li_03)
            l_05[j] = sum(li_05)
            l_1[j] = sum(li_1)
            l_2[j] = sum(li_2)
            l_5[j] = sum(li_5)
            l_10[j] = sum(li_10)

        del temp
        del density
        del log_t
        del l05
        del l01
        del l02
        del l03
        del l1
        del l2
        del l5
        del l10

        print(l_01.shape)

        return {'li_01': (l_01), 'li_02': (l_02), 'li_03': (l_03), 'li_05': (l_05), 'li_1': (l_1),
                'li_2': (l_2), 'li_5': (l_5), 'li_10': (l_10), 'sim_time': sim_time}

