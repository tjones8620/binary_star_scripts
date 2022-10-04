#!/usr/bin/python
__author__ = "Arun Mathew"

# Import the relevant standard packages ###############################
import argparse
import os
import re
import numpy as np
from astropy import units
import matplotlib.pyplot as plt
from matplotlib import ticker
import glob
from pypion.ReadData import ReadData


parser = argparse.ArgumentParser()
parser.add_argument("-p", "--path", help="Path to data folder", default="/mnt/local/thomas/wr140-cool-covertex/wr140-cool-n064")
parser.add_argument("-q", "--fluidquantity", help="Desired fluid quantity to plot", default="Density")
parser.add_argument("-s", "--surface", help ="3D Plot surface - XY, XZ, YZ", default="XY")
parser.add_argument("-imd","--img_dir", help = "Directory for image to be placed")
args = parser.parse_args()


# EDIT THIS BLOCK #######################################################
# Insert the location of simulation data (silo files)
data_path = args.path
if data_path==None:
    print("No path argument given")
    exit()
elif os.path.exists(data_path):
    print(f"Chosen datapath: {data_path}")
else:
    print(f"Directory does not exist: {data_path}")
    exit()

# # Name of your generated silo file
filename=sorted(os.listdir(data_path))[0].replace('_level00_0000.00000000.silo','')
print(filename)

# Fluid Quantity/Parameters and the corresponding tolerance
# maximum and minimum
Quantity = args.fluidquantity
print(f"Fluid quantity: {Quantity}")

# min_tolerance = np.log10(1.0e-18)
# max_tolerance = np.log10(2.0e-13)
min_tolerance = 5
max_tolerance = 9
Tolerance = [min_tolerance, max_tolerance]

# If 3D, specify the coordinate plane
surface_options = ['XY', 'XZ', 'YZ']
Surface = args.surface
if surface_options.count(args.surface) > 0:
    print(f"Chosen surface: {Surface}")
else:
    print(f"Invalid Surface Choice: {Surface}")
    # exit()

# print(Surface)
# # END OF THE BLOCK ********************************************************



############################################################################
# Cataloging silo files
os.chdir(data_path)
file_list = glob.glob('*.silo', recursive=True)
level_list = []
files = []

for file in file_list:
    level = re.search('_level(.*)_', file)
    if level == None:
        pass
    else:
        level = level.group(1)
        if not level in level_list:
            level_list.append(level)

level_list.sort()
#print('level_list =', level_list)

# Display number of levels
# Categorizing data files into levels
if len(level_list) == 1:
    print('Simulation Info: Single level')
    catalog = []
    files = sorted(glob.glob(filename + '_0000.*.silo'))
    catalog.append(files)
else:
    print(f'Simulation Info: {len(level_list)} levels')
    catalog = []
    for i in range(len(level_list)):
        files = sorted(glob.glob(f"{filename}_level{level_list[i]}_0000.*.silo"))
        catalog.append(files)

# End of Cataloging *******************************************************


###########################################################################
# Categorizing data files to time instants
# Bundle silo files of different levels of same time instant into a snapshot.

evolution = np.array(catalog).T.tolist()

# Looking inside the simulation file for dimensions.
info = ReadData(evolution[0])
N_grids = info.ngrid()
N_dims = 0
for j in range(len(N_grids)):
    if N_grids[j] > 1: N_dims += 1
print(f'Simulation Info: {N_dims}D System')
#End of making snapshots ***************************************************


############################################################################
# Making image directory
# If already exist, clean up image directory
os.chdir('../')
ImageDir = f'SimulationPlots/{filename}/{Quantity}/{Surface}'
try:
    # Create target Directory
    os.makedirs(ImageDir)  # changed mkdir to makedirs - creates nested folders
    print("Image directory" , ImageDir,  "created.")
except FileExistsError:
    print("Image directory", ImageDir, "already exists.")

    var = input("Do you want to clean up this directory? (y/n): ")
    if var=="y":
        print( "Continued with code")
    elif var=="n":
        print("Exiting script...")
        exit()
    else:
        print("Invalid input - must be y/n")
        exit()

    print('Cleaning up the',  ImageDir, 'directory ...')
    for f in os.listdir(ImageDir):
        os.remove(os.path.join(ImageDir, f))
os.chdir(data_path)
# End of making image directory *******************************************



###########################################################################
# Plot functions for different dimensions.
class Plot_Functions:
    
    #######################################################################
    # 3D-Surface-Plotter
    def ThreeDSurfacePlotter(self, snaps, fluid_quantity, tolerance, plane):

        print('3D-Surface-Plotter: Plotting', plane, 'plane.')

        if (plane == 'XY'): xlabel = 'x'; x = 2; ylabel = 'y'; y = 1 ; z = 0 
        if (plane == 'XZ'): xlabel = 'x'; x = 0; ylabel = 'z'; y = 2 ; z = 1
        if (plane == 'YZ'): xlabel = 'y'; x = 1; ylabel = 'z'; y = 0 ; z = 2

        for k in range(len(snaps)):
            data = ReadData(snaps[k])
            N_level = data.nlevels()
            N_grids = data.ngrid()
            baseline_data = data.get_3Darray(fluid_quantity)
            fluid_parameter = baseline_data['data']
            dims_min = (baseline_data['min_extents'] * units.cm).to(units.astrophys.au)
            dims_max = (baseline_data['max_extents'] * units.cm).to(units.astrophys.au)
            time = (baseline_data['sim_time'].value * units.s).to(units.yr)

            fig, ax = plt.subplots()
            ax.set_xlim(dims_min[0][x].value, dims_max[0][x].value)
            ax.set_ylim(dims_min[0][y].value, dims_max[0][y].value)
            ax.set_xlabel(xlabel + " (AU)", fontsize=8)
            ax.set_ylabel(ylabel + " (AU)", fontsize=8)

            for l in range(N_level):  # plotting each levels
                level_dims_min = dims_min[l].value
                level_dims_max = dims_max[l].value

                plot_data = np.array(fluid_parameter)
                
                # Slicer ###############################################################
                slice_location = int(0.5 * N_grids[z])
                # print(l, slice_location, ' | ', level_dims_min[1] + slice_location * dy)
                if (plane == 'YZ'): sliced_data = plot_data[l][:, :, slice_location - 1]
                if (plane == 'XZ'): sliced_data = plot_data[l][:, slice_location - 1, :]
                if (plane == 'XY'): sliced_data = plot_data[l][slice_location - 1, :, :]
                # Slicer ends #*********************************************************

            
                extents = [dims_min[l][x].value, dims_max[l][x].value, dims_min[l][y].value, dims_max[l][y].value]

                # image = ax.imshow(np.log10(sliced_data), interpolation="nearest", cmap="plasma",
                #                   extent=extents,
                #                   origin="lower",
                #                   vmin=np.log10(min(sliced_data.flatten().tolist())), vmax=np.log10(max(sliced_data.flatten().tolist()))
                #                   )

                image = ax.imshow(np.log10(sliced_data), interpolation="nearest", cmap="inferno",
                                  extent=extents,
                                  origin="lower",
                                  vmin=tolerance[0], vmax=tolerance[1]
                                  )
                
                
                plt.savefig('../' + ImageDir + '/' + 'image' + str(k).zfill(3) + '.png', bbox_inches="tight", dpi=200)
                
                if (l == 0):
                    cbaxes = fig.add_axes([0.22, 0.95, 0.575, 0.02])
                    cb = fig.colorbar(image, ax=ax, orientation="horizontal", cax=cbaxes, pad=0.0)
                    cb.ax.tick_params(labelsize=8)
                    tick_locator = ticker.MaxNLocator(nbins=5)
                    cb.locator = tick_locator
                    cb.update_ticks()

                
                    tm = str(f"{time.value:.4f}")
                    st = "\nTime = " + tm + str(time.unit) + "\n$\log_{10} \left[" + fluid_quantity + "\, (\mathrm{g\,cm}^{-3}) \\right]$"
                    # ax.text(0.05, 0.9, st, color="white", fontsize=8)

            plt.close(fig)
            print(f'Time: {time:.2e}.',
                  f'Saving snap-{str(k)} to image{str(k).zfill(3)}.png ...')
    # *3D-Plane-Plotter ends ****************************************************

#     # ############################################################################
#     # # 2D-Plotter
#     # def TwoDPlotter(self, snaps, fluid_quantity, tolerance):
#     #     print('2D-Plotter:')

#     #     for k in range(len(snaps)):
#     #         data = ReadData(snaps[k])
#     #         N_level = data.nlevels()
#     #         N_grids = data.ngrid()
#     #         baseline_data = data.get_2Darray(fluid_quantity)
#     #         # print(k, baseline_data)
#     #         fluid_parameter = baseline_data['data']
#     #         dims_min = (baseline_data['min_extents'] * units.cm).to(units.pc)
#     #         dims_max = (baseline_data['max_extents'] * units.cm).to(units.pc)
#     #         time = (baseline_data['sim_time'] * units.s).to(units.Myr)

#     #         fig, ax = plt.subplots()
#     #         ax.set_xlim(dims_min[0][0].value, dims_max[0][0].value)
#     #         ax.set_ylim(dims_min[0][1].value, dims_max[0][1].value)
#     #         ax.set_xlabel("(x pc)", fontsize=8)
#     #         ax.set_ylabel("(y pc)", fontsize=8)

#     #         for l in range(N_level):  # plotting each levels
#     #             plot_data = np.array(fluid_parameter)
#     #             extents = [dims_min[l][0].value, dims_max[l][0].value, dims_min[l][1].value, dims_max[l][1].value]

#     #             image = ax.imshow(np.log10(plot_data[l]), interpolation="nearest", cmap="viridis",
#     #                               extent=extents,
#     #                               origin="lower",
#     #                               vmin=tolerance[0], vmax=tolerance[1]
#     #                               )

#     #             plt.savefig('../' + ImageDir + '/' + str(k).zfill(5) + '.png', bbox_inches="tight", dpi=200)

#     #             if (l == 0):
#     #                 cbaxes = fig.add_axes([0.22, 0.95, 0.575, 0.02])
#     #                 cb = fig.colorbar(image, ax=ax, orientation="horizontal", cax=cbaxes, pad=0.0)
#     #                 cb.ax.tick_params(labelsize=8)
#     #                 tick_locator = ticker.MaxNLocator(nbins=5)
#     #                 cb.locator = tick_locator
#     #                 cb.update_ticks()

#     #                 tm = str("%.4f" % time.value)
#     #                 st = "\nTime = " + tm + " Myr" "\n$\log_{10} \left[" + fluid_quantity + "\, (\mathrm{g\,cm}^{-3}) \\right]$"
#     #                 ax.text(dims_min[l][0].value + 0.1, dims_max[l][1].value - 0.15, st, color="white", fontsize=8)

#     #         print(f'Time: {time:.2e}.',
#     #               'Saving snap-' + str(k) + ' to ' + 'image' + str(k).zfill(3) + '.png ...')
#     #         plt.close(fig)

#     #     return True
#     # # 2D-Plotter ends ****************************************************

#     # ############################################################################
#     # # 1D-Plotter
#     # def OneDPlotter(self, snaps, fluid_quantity, tolerance):

#     #     print('1D-Plotter:')
#     #     pass


plot = Plot_Functions()
# if N_dims == 1: plot.OneDPlotter(evolution, Quantity, Tolerance)
# if N_dims == 2: plot.TwoDPlotter(evolution, Quantity, Tolerance)
if N_dims == 3: plot.ThreeDSurfacePlotter(evolution[::20], Quantity, Tolerance, Surface)









