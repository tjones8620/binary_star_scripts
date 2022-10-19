# Author: Sam Green, Created: 12-08-2019
# Analyses 3D data.

# Fixes:
# 13-11-2019: Changed script so that it can be called using the new bash script, cal-lum.
# 13-11-2019: Tidied the script up a bit, initialised arrays (made script 20-30 secs faster).
# 13-11-2019: Data is saved to a txt file after every loop in the bash script.

# (1) install pypion with `python3 -m pip install pypion`
# (2) install python3-silo from https://git.dias.ie/massive-stars-software/pypion.git
#     - git clone https://git.dias.ie/massive-stars-software/pypion.git
#     - cd pypion/src/silo/
#     - bash install_silo.sh


# Arguments to this script:
# <path-to-files>: where the silo files are located
# <base-filename>: filename, excluding `_level00_0000.*.silo`
# <image-path>: not used here, use `./`
# <image-filename>: filename for text file containing output info.
# <dimension>: number of dimensions in PION simulation: 1, 2 or 3

from pathlib import Path
home = str(Path.home())
import sys
sys.path.insert(0,home+"/.local/silo/lib")

import Silo
from Analysis import Analysis
from pypion.argparse_command import InputValues

import numpy as np
import time

cmdargs = InputValues()
time_dicts = cmdargs.time_dicts
# cmdargs.img_file : image base filename
# time_dicts : list of files, somehow.


# opens a txt file to save the data in the following for loop.
dfile="wind-mass-"+cmdargs.img_file+".txt"
data_file = open(dfile, "w")

t0 = time.time()

# loops through all files
i=0
for files in time_dicts:

    t2 = time.time()

    arr = time_dicts[files]
    plotting = Analysis(arr)

    # calls the luminosity function from Analysis_3D.py.
    masses = plotting.wind_mass_2D()

    # saves the data generated from the luminosity funcitons into an array.
    m1e4 = masses['m1e4']
    m3e4 = masses['m3e4']
    m1e5 = masses['m1e5']
    m3e5 = masses['m3e5']
    m1e6 = masses['m1e6']
    m3e6 = masses['m3e6']
    m1e7 = masses['m1e7']
    m1e9 = masses['m1e9']
    sim_time = masses['sim_time']

    # close Analysis_3D.py to clear memory.
    plotting.close()
    # writes the data to the txt file.
    data_file.write("%13.5e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e\n" % (sim_time, np.sum(m1e4), np.sum(m3e4), np.sum(m1e5), np.sum(m3e5), np.sum(m1e6), np.sum(m3e6), np.sum(m1e7), np.sum(m1e9)))

    t3 = time.time()
    print(i,"Running Time per timestep: {0}".format(t3 - t2))
    i = i+1

t1 = time.time()
print("Overall Running Time: {0}".format(t1 - t0))
data_file.close()
