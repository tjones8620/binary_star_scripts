#!/usr/bin/python

__author__ = "Arun Mathew, Thomas Jones"

# Import the relevant standard packages ###############################
import os
import re
import numpy as np
from astropy import units
import matplotlib.pyplot as plt
from matplotlib import ticker
import glob
import moviepy.video.io.ImageSequenceClip
from pypion.ReadData import ReadData

###########################################################################
# Plot functions for different dimensions.
class Plot_Functions:
    def __init__(self, path, image_dir, fluid_quantity, tolerance, plane):
        ########### Defining path to SILO files #######################
        self.data_path = path
        if os.path.exists(path):
            print(f"Chosen datapath: {path}")
        else:
            print(f"Directory does not exist: {path}")
            exit()
            
        ########### Base name of SILO files ############################
        self.filename=sorted(os.listdir(self.data_path))[0].replace('_level00_0000.00000000.silo','')
        print(self.filename)

        ########### Define desired fluid quantity ######################
        self.Quantity = fluid_quantity
        print(f"Fluid quantity: {self.Quantity}")
        self.Tolerance = tolerance
        print(f"vmin: {self.Tolerance[0]} \nvmax: {self.Tolerance[1]}")

        ########### If 3D, specify the coordinate plane ################
        surface_options = ['XY', 'XZ', 'YZ']
        self.Surface = plane
        if surface_options.count(self.Surface) > 0:
            print(f"Chosen surface: {self.Surface}")
        else:
            print(f"Invalid Surface Choice: {self.Surface}")

        ########### Make snapshots of the simulation ###################
        self.evolution = self.make_snapshots(self.data_path, self.filename)

        print(self.evolution.shape)
        ########### Find the dimensions of the simulation ##############
        self.find_dimensions()
        
        ########### Create image directory #############################
        self.ImageDir = self.create_image_dir(image_dir, self.filename, self.Quantity, self.Surface)

        ######### Defining default fps for mp4 ##########
        self.fps=10

    @staticmethod
    def make_snapshots(data_path, filename):
        ########## Cataloging silo files ###############################
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

        ########## Categorizing data files into levels #################
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
                
        # Bundle silo files of different levels of same time instant into a snapshot.
        evolution = np.array(catalog).T
        print(f"Shape of evolution array: {evolution.shape}")
        return evolution

    def find_dimensions(self):
        # Looking inside the simulation file for dimensions.
        info = ReadData(self.evolution[0])
        N_grids = info.ngrid()
        self.N_dims = 0
        for j in range(len(N_grids)):
            if N_grids[j] > 1: self.N_dims += 1
        print(f'Simulation Info: {self.N_dims}D System')
        #End of making snapshots ***************************************************


    @staticmethod
    def create_image_dir(img_dir, *args):
        # ImageDir = os.path.join(img_dir, f'SimulationPlots/{self.filename}/{self.Quantity}/{self.Surface}')
        ImageDir = os.path.join(img_dir, f'SimulationPlots/{args[0]}/{args[1]}/{args[2]}')
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
            print('Directory cleaned up.')
        return ImageDir




    # #######################################################################
    # # 3D-Surface-Plotter
    # def ThreeDSurfacePlotter(self, colormap='viridis', movie=False):

    #     print('3D-Surface-Plotter: Plotting', self.Surface, 'plane')

    #     if (self.Surface == 'XY'): xlabel = 'x'; x = 0; ylabel = 'y'; y = 1 ; z = 2 
    #     if (self.Surface == 'XZ'): xlabel = 'x'; x = 0; ylabel = 'z'; y = 2 ; z = 1
    #     if (self.Surface == 'YZ'): xlabel = 'y'; x = 1; ylabel = 'z'; y = 2 ; z = 1

    #     for k in range(len(self.evolution)):
    #         data = ReadData(self.evolution[k])
    #         N_level = data.nlevels()
    #         N_grids = data.ngrid()
    #         baseline_data = data.get_3Darray(self.Quantity)
    #         fluid_parameter = baseline_data['data']
    #         dims_min = (baseline_data['min_extents'] * units.cm).to(units.astrophys.au)
    #         dims_max = (baseline_data['max_extents'] * units.cm).to(units.astrophys.au)
    #         time = (baseline_data['sim_time'].value * units.s).to(units.yr)

    #         fig, ax = plt.subplots()
    #         ax.set_xlim(dims_min[0][x].value, dims_max[0][x].value)
    #         ax.set_ylim(dims_min[0][y].value, dims_max[0][y].value)
    #         ax.set_xlabel(xlabel + str(dims_min.unit), fontsize=8)
    #         ax.set_ylabel(ylabel + str(dims_min.unit), fontsize=8)

    #         for l in range(N_level):  # plotting each levels

    #             plot_data = np.array(fluid_parameter)
                
    #             # Slicer ###############################################################
    #             slice_location = int(0.5 * N_grids[z])
    #             if (self.Surface == 'YZ'): sliced_data = plot_data[l][:, :, slice_location - 1]
    #             if (self.Surface == 'XZ'): sliced_data = plot_data[l][:, slice_location - 1, :]
    #             if (self.Surface == 'XY'): sliced_data = plot_data[l][slice_location - 1, :, :]
    #             # Slicer ends #*********************************************************

            
    #             extents = [dims_min[l][x].value, dims_max[l][x].value, dims_min[l][y].value, dims_max[l][y].value]

    #             # if (self.Surface == 'YZ'): extents = [dims_min[l][y].value, dims_max[l][y].value, dims_min[l][x].value, dims_max[l][x].value]

    #             image = ax.imshow(np.log10(sliced_data), interpolation="nearest", cmap=colormap,
    #                               extent=extents,
    #                               origin="lower",
    #                               vmin=self.Tolerance[0], vmax=self.Tolerance[1]
    #                               )
                
                
    #             plt.savefig(f"{self.ImageDir}/image{str(k).zfill(3)}.png", bbox_inches="tight", dpi=200)
                
    #             if (l == 0):
    #                 cbaxes = fig.add_axes([0.22, 0.95, 0.575, 0.02])
    #                 cb = fig.colorbar(image, ax=ax, orientation="horizontal", cax=cbaxes, pad=0.0)
    #                 cb.ax.tick_params(labelsize=8)
    #                 tick_locator = ticker.MaxNLocator(nbins=5)
    #                 cb.locator = tick_locator
    #                 cb.update_ticks()

                
    #                 # tm = str(f"{time.value:.4f}")
    #                 # st = "\nTime = " + tm + str(time.unit) + "\n$\log_{10} \left[" + fluid_quantity + "\, (\mathrm{g\,cm}^{-3}) \\right]$"
    #                 # ax.text(0.3, 0.9, st, color="white", fontsize=8, transform=ax.transAxes)

    #         plt.close(fig)
    #         print(f'Time: {time:.2e}.',
    #               f'Saving snap-{str(k)} to image{str(k).zfill(3)}.png ...')

    
    #     if movie==True:

    #         image_folder=self.ImageDir
    #         image_files = [os.path.join(image_folder,img)
    #                     for img in sorted(os.listdir(image_folder))
    #                     if img.endswith(".png")]
    #         clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=self.fps)
    #         clip.write_videofile(f'{self.ImageDir}/{self.filename}_{self.Quantity}_{self.Surface}.mp4')
        
    #     else:
    #         pass
        

