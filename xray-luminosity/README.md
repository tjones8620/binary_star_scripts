
Files:
=============

+ `xray-table.txt`: emmissivity of gas in collisional ionization equilibrium.
  This is calculated using XSPEC for gas with Solar System abundances as 
  described in Green et al. (2022) https://ui.adsabs.harvard.edu/abs/2022arXiv220306331G/abstract
+ `Analysis.py`: python functions to calculate an array of cell volumes for 1D,
  2D and 3D grids, assuming 1D = spherical coordinates, 2D = cylindrical (z,R),
  and 3D = Cartesian.  Also calculate X-ray luminosity on these grids, using
  the `xray-table.txt` as input data.  Also calculate the total gas mass in
  temperature bins.
+ `total_lum.py` takes in a list of input Silo files, calculates X-ray
  luminsosity for each of them, and saves result in an ascii text file.
+ `wind-mass.py` takes in a list of input Silo files, calculates the gas mass
  in different temperature ranges for each file, and saves the result in a text
  file.
+ `plot_2D_lum-mass.py`: plot the results from the x-ray luminosity for
  different spatial resolutions.
+ `plot_3Dlum.py`: plot results for 3D MHD sims of Zeta Ophiuchi (Sam Green)
+ `run-xray-2dmhd-analysis.sh`: run all the scripts for 2D hydro and mhd
  simulations, to make plots of X-ray luminosity vs. time for different
  spatial resolutions and different solvers.

