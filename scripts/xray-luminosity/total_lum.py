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

# from pathlib import Path
# home = str(Path.home())
# import sys
# sys.path.insert(0,home+"/.local/silo/lib")

# import Silo
from Analysis import Analysis
from pypion.argparse_command import InputValues

import numpy as np
import time

cmdargs = InputValues()
time_dicts = cmdargs.time_dicts

#print(dir(cmdargs))
#print(cmdargs.__dict__)
dd = cmdargs.dimen
print(dd)

# log10(T/K).0  T(K).1  E(keV).2  j(E>0.1keV).3   j(E>0.2keV).4   j(E>0.5keV).5 /
# j(E>1keV).6   j(E>2keV).7   j(E>5keV).8   j(E>10keV).9
# Create an empty array:
log_t = []
l01 = []
l02 = []
l03 = []
l05 = []
l1 = []
l2 = []
l5 = []
l10 = []

# Read a line of numbers out of a text file:
with open("xray-table.txt") as x:
    for line in x:
        data = line.split()
        if data[0] == "#": 
          continue
        log_t.append(float(data[0]))
        l01.append(float(data[3]))
        l02.append(float(data[4]))
        l03.append(float(data[5]))
        l05.append(float(data[6]))
        l1.append(float(data[7]))
        l2.append(float(data[8]))
        l5.append(float(data[9]))
        l10.append(float(data[10]))

# opens a txt file to save the data in the following for loop.
dfile="xray-lum-"+cmdargs.img_file+".txt"
data_file = open(dfile, "w")

t0 = time.time()

# loops through all files
i=0
for files in time_dicts:

    t2 = time.time()

    arr = time_dicts[files]
    plotting = Analysis(arr)

    # calls the luminosity function from Analysis_3D.py.
    if dd=="2d":
      lum = plotting.luminosity_interp2D(log_t, l01, l02, l03, l05, l1, l2, l5, l10)
    elif dd=="1d":
      lum = plotting.luminosity_interp1D(log_t, l01, l02, l03, l05, l1, l2, l5, l10)
    elif dd=="3d":
      lum = plotting.luminosity_interp3D(log_t, l01, l02, l03, l05, l1, l2, l5, l10)
    else:
      print("bad dimension input on commandline",dd)
      quit()
      

    # saves the data generated from the luminosity funcitons into an array.
    l01i = lum['li_01']
    l02i = lum['li_02']
    l03i = lum['li_03']
    l05i = lum['li_05']
    l1i = lum['li_1']
    l2i = lum['li_2']
    l5i = lum['li_5']
    l10i = lum['li_10']
    sim_time = lum['sim_time']

    # close Analysis_3D.py to clear memory.
    plotting.close()
    # writes the data to the txt file.
    data_file.write("%s %s %s %s %s %s %s %s %s\n" % (sim_time.value, l01i, l02i, l03i, l05i, l1i, l2i, l5i, l10i))

    t3 = time.time()
    print(i,"Running Time per timestep: {0}".format(t3 - t2))
    i = i+1

t1 = time.time()
print("Overall Running Time: {0}".format(t1 - t0))
data_file.close()
