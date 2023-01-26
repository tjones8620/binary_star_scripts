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

        return {'li_01': sum(l_01), 'li_02': sum(l_02), 'li_03': sum(l_03), 'li_05': sum(l_05), 'li_1': sum(l_1),
                'li_2': sum(l_2), 'li_5': sum(l_5), 'li_10': sum(l_10), 'sim_time': sim_time}

###############################################################################
    def luminosity_interp2D(self, log_t, l01, l02, l03, l05, l1, l2, l5, l10):
        # This function uses an interpolation method to take known values from
        # a .txt to calculate the x-ray luminosity of each cell
        D = self.get_2Darray("Density")
        density = D['data']
        temp = self.get_2Darray("Temperature")['data']
        mask = self.get_2Darray("NG_Mask")['data']

        ngrid = self.ngrid()
        lim_max = (D['max_extents'] * u.cm).value
        lim_min = (D['min_extents'] * u.cm).value
        sim_time = D['sim_time']

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
            #print(lim_max[j], lim_min[j], ngrid)
            volj = self.volume2D(lim_max[j], lim_min[j], ngrid)

            log_masktemp = np.log10(tempj) * maskj
            n_ei = denj / (1.167 * m_p)
            n_hi = denj / (1.4 * m_p)

            # Use an interpolation function to calculate the x-ray luminosity at different energy levels:
            l_01[j] = sum(sum(np.interp(log_masktemp, log_t, l01, left=0) * volj * n_ei * n_hi))  # L(E>0.1kev)
            l_02[j] = sum(sum(np.interp(log_masktemp, log_t, l02, left=0) * volj * n_ei * n_hi))  # L(E>0.2kev)
            l_03[j] = sum(sum(np.interp(log_masktemp, log_t, l03, left=0) * volj * n_ei * n_hi))  # L(E>0.3kev)
            l_05[j] = sum(sum(np.interp(log_masktemp, log_t, l05, left=0) * volj * n_ei * n_hi))  # L(E>0.5kev)
            l_1[j] = sum(sum(np.interp(log_masktemp, log_t, l1, left=0) * volj * n_ei * n_hi))  # L(E>1kev)
            l_2[j] = sum(sum(np.interp(log_masktemp, log_t, l2, left=0) * volj * n_ei * n_hi))  # L(E>2kev)
            l_5[j] = sum(sum(np.interp(log_masktemp, log_t, l5, left=0) * volj * n_ei * n_hi))  # L(E>5kev)
            l_10[j] = sum(sum(np.interp(log_masktemp, log_t, l10, left=0)*  volj * n_ei * n_hi))  # L(E>10kev)
            del log_masktemp, n_ei, volj

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

        return {'li_01': sum(l_01), 'li_02': sum(l_02), 'li_03': sum(l_03), 'li_05': sum(l_05), 'li_1': sum(l_1),
                'li_2': sum(l_2), 'li_5': sum(l_5), 'li_10': sum(l_10), 'sim_time': sim_time}

###############################################################################
    def luminosity_interp1D(self, log_t, l01, l02, l03, l05, l1, l2, l5, l10):
        # This function uses an interpolation method to take known values from
        # a .txt to calculate the x-ray luminosity of each cell

        data = self.get_1Darray("Density")
        density = data['data']
        temp = self.get_1Darray("Temperature")['data']
        #mask = self.get_1Darray("NG_Mask")['data']

        ngrid = self.ngrid()

        lim_max = (data['max_extents'] * u.cm).value
        lim_min = (data['min_extents'] * u.cm).value

        sim_time = data['sim_time']

        l_001 = np.zeros(len(density))
        l_002 = np.zeros(len(density))
        l_003 = np.zeros(len(density))
        l_005 = np.zeros(len(density))
        l_010 = np.zeros(len(density))
        l_020 = np.zeros(len(density))
        l_050 = np.zeros(len(density))
        l_100 = np.zeros(len(density))

        # loop over grid levels (only 1 level in 1D)
        for j in range(len(density)):
            denj = density[j]
            tempj = temp[j]
            #maskj = mask[j]
            log_masktemp = np.log10(tempj)

            vol = self.volume1D(lim_max[j], lim_min[j], ngrid)

            li_001 = np.zeros(ngrid[0])
            li_002 = np.zeros(ngrid[0])
            li_003 = np.zeros(ngrid[0])
            li_005 = np.zeros(ngrid[0])
            li_010 = np.zeros(ngrid[0])
            li_020 = np.zeros(ngrid[0])
            li_050 = np.zeros(ngrid[0])
            li_100 = np.zeros(ngrid[0])
            
            # electron number density
            n_ei = denj * (1 / (1.167 * m_p))
            # ion number density
            n_hi = denj * (1 / (1.40 * m_p))
            li_001 = np.interp(log_masktemp, log_t, l01, left=0) * vol * n_ei * n_hi
            li_002 = np.interp(log_masktemp, log_t, l02, left=0) * vol * n_ei * n_hi
            li_003 = np.interp(log_masktemp, log_t, l03, left=0) * vol * n_ei * n_hi
            li_005 = np.interp(log_masktemp, log_t, l05, left=0) * vol * n_ei * n_hi
            li_010 = np.interp(log_masktemp, log_t, l1, left=0) * vol * n_ei * n_hi
            li_020 = np.interp(log_masktemp, log_t, l2, left=0) * vol * n_ei * n_hi
            li_050 = np.interp(log_masktemp, log_t, l5, left=0) * vol * n_ei * n_hi
            li_100 = np.interp(log_masktemp, log_t, l10, left=0) * vol * n_ei * n_hi

            l_001[j] = sum(li_001)
            l_002[j] = sum(li_002)
            l_003[j] = sum(li_003)
            l_005[j] = sum(li_005)
            l_010[j] = sum(li_010)
            l_020[j] = sum(li_020)
            l_050[j] = sum(li_050)
            l_100[j] = sum(li_100)

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

        return {'li_01': sum(l_001), 'li_02': sum(l_002), 'li_03': sum(l_003), 'li_05': sum(l_005), 'li_1': sum(l_010),
                'li_2': sum(l_020), 'li_5': sum(l_050), 'li_10': sum(l_100), 'sim_time': sim_time}

###############################################################################
    def wind_mass_2D(self):
        # This function calculates the mass in stellar wind material at
        # temperatures less than certain values.
        D = self.get_2Darray("Density")
        density = D['data']
        temp = self.get_2Darray("Temperature")['data']
        mask = self.get_2Darray("NG_Mask")['data']
        wind = self.get_2Darray("Tr000_WIND")['data']

        ngrid = self.ngrid()
        lim_max = (D['max_extents'] * u.cm).value
        lim_min = (D['min_extents'] * u.cm).value
        sim_time = D['sim_time']

        ml1e4 = np.zeros(len(density))
        ml3e4 = np.zeros(len(density))
        ml1e5 = np.zeros(len(density))
        ml3e5 = np.zeros(len(density))
        ml1e6 = np.zeros(len(density))
        ml3e6 = np.zeros(len(density))
        ml1e7 = np.zeros(len(density))
        ml1e9 = np.zeros(len(density))

        # loop over levels in nested grids
        for j in range(len(density)):
            densj = density[j]
            tempj = temp[j]
            maskj = mask[j]
            windj = wind[j]
            volj = self.volume2D(lim_max[j], lim_min[j], ngrid)

            # This is the mass of wind material in each cell (g)
            maskmass = densj * maskj *windj * volj

            # create mask for cells: 0 if T>T_lim, 1 otherwise, so we only
            # cells with temperature less than T_lim
            t1 = np.where(tempj > 1.0e4, 0, 1);
            ml1e4[j] = np.sum(maskmass * t1)
            t2 = np.where(tempj > 3.0e4, 0, 1);
            ml3e4[j] = np.sum(maskmass * t2)
            t3 = np.where(tempj > 1.0e5, 0, 1);
            ml1e5[j] = np.sum(maskmass * t3)
            t4 = np.where(tempj > 3.0e5, 0, 1);
            ml3e5[j] = np.sum(maskmass * t4)
            t5 = np.where(tempj > 1.0e6, 0, 1);
            ml1e6[j] = np.sum(maskmass * t5)
            t6 = np.where(tempj > 3.0e6, 0, 1);
            ml3e6[j] = np.sum(maskmass * t6)
            t7 = np.where(tempj > 1.0e7, 0, 1);
            ml1e7[j] = np.sum(maskmass * t7)
            t8 = np.where(tempj > 1.0e9, 0, 1);
            ml1e9[j] = np.sum(maskmass * t8)
            del maskmass, t1, t2, t3, t4, t5, t6, t7, t8
        del temp
        del density

        return {'m1e4': ml1e4, 'm3e4': ml3e4, 'm1e5': ml1e5, 'm3e5': ml3e5, 'm1e6': ml1e6, 'm3e6': ml3e6, 'm1e7': ml1e7, 'm1e9': ml1e9, 'sim_time': sim_time}


