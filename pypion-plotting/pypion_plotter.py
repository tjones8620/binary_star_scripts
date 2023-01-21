#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
This module contains the Plot_Functions class which is used to plot the simulation data.

The Plot_Functions class contains methods for plotting the simulation data in different ways.

The class is initialized with the path to the simulation data, the desired fluid quantity, the
tolerance for the colorbar, the coordinate plane, and the image directory. 

The class then finds the dimensions of the simulation and creates a directory for the images. 

The class can also be used to make movies of the simulation
data. The class can also be used to plot the orbital phase of the simulation data. 
"""

__author__ = "Thomas Jones"

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
plt.style.use("science")

###########################################################################
# Plot functions for different dimensions.
class Plot_Functions:
    quantity_dict = {'Density':r'$\rho$', 'Temperature':"T", 'Velocity': "v", 'Tr000_WIND': "Wind Tracer"}
    unit_dict = {'Density': "$g \, cm^{-3}$", 'Temperature': "K", 'Velocity': "$cm \, s^{-1}$", "Tr000_WIND": ""}

    def __init__(self, path, image_dir, fluid_quantity="Density", tolerance=[-20,13], plane="XY", start_time=0, period=7.992):
        ########### Defining path to SILO files #######################
        self.data_path = path
        if os.path.exists(path):
            print(f"Chosen datapath: {path}")
        else:
            print(f"Directory does not exist: {path}")
            exit()
            
        ########### Base name of SILO files ############################
        _ = os.listdir(path)
        _= sorted([f for f in _ if f.endswith('.silo')])
        self.filename=(_[0]).replace('_level00_0000.00000000.silo','')
        # self.filename=sorted(os.listdir(self.data_path))[4].replace('_level00_0000.00000000.silo','')
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
            exit()

        ########### Make snapshots of the simulation ###################
        self.evolution = self.make_snapshots(self.data_path, self.filename)

        ########### Find the dimensions of the simulation ##############
        self.find_dimensions()
        
        ########### Create image directory #############################
        self.ImageDir = self.create_image_dir(image_dir, self.filename, self.Quantity, self.Surface)

        ######### Defining default fps for mp4 ##########
        self.fps=10

        self.start_time = start_time
        self.period = period*units.yr

    @staticmethod
    def make_snapshots(data_path, filename):
        ########## Cataloging silo files ###############################
        # os.chdir(data_path)
        file_list = glob.glob(os.path.join(data_path, '*.silo'), recursive=True)
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
                files = sorted(glob.glob(os.path.join(data_path, f"{filename}_level{level_list[i]}_0000.*.silo")))
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
        ImageDir = os.path.join(img_dir, f'SimulationPlots/{args[0]}/{args[1]}/{args[2]}')

        if not os.path.exists(ImageDir):
            os.makedirs(ImageDir)
            print("Image directory" , ImageDir,  "created.")
        else:
            print("Image directory", ImageDir, "already exists.")
        return ImageDir

    #######################################################################
    # 3D-Surface-Plotter
    def ThreeDSurfacePlotter(self, colormap='viridis', movie=False, log=True):

        """
        This function plots the 3D surface of the simulation.

        Parameters
        ----------
        colormap : str, optional
            The colormap to be used for the plot. The default is 'viridis'.
        movie : bool, optional
            If True, a movie of the simulation will be created. The default is False.
        log : bool, optional
            If True, the plot will be in log scale. The default is True.
        
        Returns
        -------
        None.
        """

        print('3D-Surface-Plotter: Plotting', self.Surface, 'plane')

        if (self.Surface == 'XY'): xlabel = 'x'; x = 0; ylabel = 'y'; y = 1 ; z = 2 
        if (self.Surface == 'XZ'): xlabel = 'x'; x = 0; ylabel = 'z'; y = 2 ; z = 1
        if (self.Surface == 'YZ'): xlabel = 'y'; x = 1; ylabel = 'z'; y = 2 ; z = 1

        # Looping over all the snapshots
        for k in range(len(self.evolution)):
            data = ReadData(self.evolution[k])
            N_level = data.nlevels()
            N_grids = data.ngrid()
            baseline_data = data.get_3Darray(self.Quantity)
            fluid_parameter = baseline_data['data']
            dims_min = (baseline_data['min_extents'] * units.cm).to(units.astrophys.au)
            dims_max = (baseline_data['max_extents'] * units.cm).to(units.astrophys.au)
            time = ((baseline_data['sim_time'].value + self.start_time) * units.s).to(units.yr)

            fig, ax = plt.subplots()
            ax.set_xlim(dims_min[0][x].value, dims_max[0][x].value)
            ax.set_ylim(dims_min[0][y].value, dims_max[0][y].value)
            ax.set_xlabel(f"{xlabel} ({str(dims_min.unit)})", fontsize=8)
            ax.set_ylabel(f"{ylabel} ({str(dims_min.unit)})", fontsize=8)

            # Looping over all the levels
            for l in range(N_level): 

                plot_data = np.array(fluid_parameter)
                
                # Slicer ###############################################################
                slice_location = int(0.5 * N_grids[z])
                if (self.Surface == 'YZ'): sliced_data = plot_data[l][:, :, slice_location - 1]
                if (self.Surface == 'XZ'): sliced_data = plot_data[l][:, slice_location - 1, :]
                if (self.Surface == 'XY'): sliced_data = plot_data[l][slice_location - 1, :, :]
                # Slicer ends #*********************************************************
            
                extents = [dims_min[l][x].value, dims_max[l][x].value, dims_min[l][y].value, dims_max[l][y].value]

                if (log == True):
                    image = ax.imshow(np.log10(sliced_data), interpolation="nearest", cmap=colormap,
                                      extent=extents,
                                      origin="lower",
                                      vmin=self.Tolerance[0], vmax=self.Tolerance[1]
                                      )

                elif (log == False):
                    image = ax.imshow(sliced_data, interpolation="nearest", cmap=colormap,
                                      extent=extents,
                                      origin="lower",
                                      vmin=10**self.Tolerance[0], vmax=10**self.Tolerance[1]
                                      )
                else:
                    print('Error: log must be True or False')
                    exit()


            ######### Colorbar ########################################################
            cbaxes = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
            cb = fig.colorbar(image, ax=ax, orientation="vertical", cax=cbaxes, pad=0.0)
            if (self.Quantity == 'Density'): 
                cb.set_label(r"$\log_{10}$ "+r"$(\rho)$" + f" ({self.unit_dict[self.Quantity]})", fontsize=8, labelpad=2)
            else:
                cb.set_label(r"$\log_{10}$ " + f"{self.Quantity} " + f"({self.unit_dict[self.Quantity]})", fontsize=8, labelpad=2)
            cb.ax.tick_params(labelsize=8)
            tick_locator = ticker.MaxNLocator(nbins=5)
            cb.locator = tick_locator
            cb.update_ticks()

            ######### Phase ###########################################################
            phase = str(f"{((self.period + time)/self.period):.2f}")
            st = r"$\phi$ = " + phase 
            ax.text(0.8, 0.9, st, color="black", fontsize=8, transform=ax.transAxes, 
                    bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=1', alpha=0.5))

            plt.savefig(f"{self.ImageDir}/image{str(k).zfill(3)}.png", bbox_inches="tight", dpi=500)
            plt.close(fig)
            print(f'Time: {time:.2e}.',
                  f'Saving snap-{str(k)} to image{str(k).zfill(3)}.png ...')

    
        #### Movie ###############################################################
        if movie==True:

            image_folder=self.ImageDir
            image_files = [os.path.join(image_folder,img)
                        for img in sorted(os.listdir(image_folder))
                        if img.endswith(".png")]
            clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=self.fps)
            clip.write_videofile(f'{self.ImageDir}/{self.filename}_{self.Quantity}_{self.Surface}.mp4')
        
        else:
            pass
    
    def three_time_slice(self, colormap='viridis', d_phase = 0.01, log=True, zoom=1, plot_inset=False):

        """
        Function to plot three time slices of the simulation: pre-periastron, periastron, and post-periastron.
        The time slices are determined by the phase of the simulation. The phase is calculated by dividing the
        simulation time by the orbital period. 

        Parameters
        ----------
        colormap : str, optional
            Colormap to use for the plot. The default is 'viridis'.
        
        d_phase : float, optional
            The phase difference between the pre-periastron, periastron, and post-periastron snapshots.
            The default is 0.01.
        
        log : bool, optional
            If True, the plot will be in log scale. The default is True.
        
        zoom : int, optional
            Zoom factor for the plot. The default is 1.

        plot_inset : bool, optional
            If True, the plot will have an inset of the simulation domain. The default is False.
        
        Returns
        -------
        None.
        """

        if (self.Surface == 'XY'): xlabel = 'x'; x = 0; ylabel = 'y'; y = 1 ; z = 2 
        if (self.Surface == 'XZ'): xlabel = 'x'; x = 0; ylabel = 'z'; y = 2 ; z = 1
        if (self.Surface == 'YZ'): xlabel = 'y'; x = 1; ylabel = 'z'; y = 2 ; z = 1

        pre_periastron_list = []  # list of pre-periastron snapshots
        pre_periastron_phase = [] # list of pre-periastron times
        periastron_list = [] # list of periastron snapshots
        periastron_phase = [] # list of periastron times
        post_periastron_list = [] # list of post-periastron snapshots
        post_periastron_phase = [] # list of post-periastron times

        for k in range(len(self.evolution)):
            data = ReadData(self.evolution[k])
            sim_time = data.sim_time().value
            data.close()
            time = ((sim_time + self.start_time) * units.s).to(units.yr)
            phase = ((self.period + time)/self.period).value
            phase_str = f"{((self.period + time)/self.period):.2f}"

            pre_per_phase = 1 - d_phase
            post_per_phase = 1 + d_phase

            if round(phase, len(str(d_phase))-2) == round(1, len(str(d_phase))-2):
                periastron_list.append(k)
                periastron_phase.append(abs(phase - 1.0))
            elif round(phase, len(str(d_phase))-2) == round(pre_per_phase, len(str(d_phase))-2):
                pre_periastron_list.append(k)
                pre_periastron_phase.append(abs(phase - pre_per_phase))
            elif round(phase, len(str(d_phase))-2) == round(post_per_phase, len(str(d_phase))-2):
                post_periastron_list.append(k)
                post_periastron_phase.append(abs(phase - post_per_phase))
                        
        
        plot_indices = [pre_periastron_list[np.argmin(pre_periastron_phase)], periastron_list[np.argmin(periastron_phase)], post_periastron_list[np.argmin(post_periastron_phase)]]

        threeplotfig, ax = plt.subplots(1, 3, figsize=(12, 4), sharey=True)
        num = 0
        for index in plot_indices:
            data = ReadData(self.evolution[index])
            N_level = data.nlevels()
            N_grids = data.ngrid()
            baseline_data = data.get_3Darray(self.Quantity)
            data.close()
            fluid_parameter = baseline_data['data']
            dims_min = (baseline_data['min_extents'] * units.cm).to(units.astrophys.au)
            dims_max = (baseline_data['max_extents'] * units.cm).to(units.astrophys.au)
            time = ((baseline_data['sim_time'].value + self.start_time) * units.s).to(units.yr)

            ax[num].set_xlim(dims_min[0][x].value/zoom, dims_max[0][x].value/zoom)
            ax[num].set_ylim(dims_min[0][y].value/zoom, dims_max[0][y].value/zoom)
            ax[num].set_xlabel(f"{xlabel} ({str(dims_min.unit)})", fontsize=8)

            if plot_inset==True:
                axins = ax[num].inset_axes([0.55, 0, 0.45, 0.45], transform=ax[num].transAxes)
                axins.set_xlim(dims_min[0][x].value/(zoom*16), dims_max[0][x].value/(zoom*16))
                axins.set_ylim(dims_min[0][y].value/(zoom*16), dims_max[0][y].value/(zoom*16))
                axins.set_xticks([])
                axins.set_yticks([])
                axins.set_xlabel("")
                axins.set_ylabel("")
        
            # ax.indicate_inset_zoom(axins, edgecolor="black")


            for l in range(N_level):  # plotting each levels

                plot_data = np.array(fluid_parameter)
                
                # Slicer ###############################################################
                slice_location = int(0.5 * N_grids[z])
                if (self.Surface == 'YZ'): sliced_data = plot_data[l][:, :, slice_location - 1]
                if (self.Surface == 'XZ'): sliced_data = plot_data[l][:, slice_location - 1, :]
                if (self.Surface == 'XY'): sliced_data = plot_data[l][slice_location - 1, :, :]
                # Slicer ends #*********************************************************
            
                extents = [dims_min[l][x].value, dims_max[l][x].value, dims_min[l][y].value, dims_max[l][y].value]

                if log == True:
                    image = ax[num].imshow(np.log10(sliced_data), interpolation="nearest", cmap=colormap,
                                    extent=extents,
                                    origin="lower",
                                    vmin=self.Tolerance[0], vmax=self.Tolerance[1]
                                    )

                    label = r"$\log_{10}$" + f"({self.quantity_dict[self.Quantity]}) "   

                    if plot_inset==True:
                        axins.imshow(np.log10(sliced_data), interpolation="nearest", cmap=colormap,
                                    extent=extents,
                                    origin="lower",
                                    vmin=self.Tolerance[0], vmax=self.Tolerance[1]
                                    )

                else:
                    image = ax[num].imshow(sliced_data, interpolation="nearest", cmap=colormap,
                                    extent=extents,
                                    origin="lower",
                                    vmin=self.Tolerance[0], vmax=self.Tolerance[1]
                                    )                  

                    label = fr"{self.quantity_dict[self.Quantity]} " 

            if not self.Quantity == 'Tr000_WIND':
                label = label + f"({self.unit_dict[self.Quantity]})" 
            
            ax[num].indicate_inset_zoom(axins, edgecolor="black")

            phase = str(f"{((self.period + time)/self.period):.{len(str(d_phase))-2}f}")
            st = r"$\phi$ = " + phase 
            ax[num].text(0.7, 0.9, st, color="black", fontsize=8, transform=ax[num].transAxes, 
                    bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=1', alpha=0.5))
    
            num+=1
            print(num)
        
        ax[0].set_ylabel(f"{ylabel} ({str(dims_min.unit)})", fontsize=8)
        cbaxes = threeplotfig.add_axes([ax[num-1].get_position().x1+0.01,ax[num-1].get_position().y0,0.02,ax[num-1].get_position().height])
        cb = threeplotfig.colorbar(image, ax=ax[num-1], orientation="vertical", cax=cbaxes, pad=3.0)
        if (self.Quantity == 'Density'): 
            cb.set_label(label, fontsize=8, labelpad=4)
        else:
            cb.set_label(label, fontsize=8, labelpad=2)
        tick_locator = ticker.MaxNLocator(nbins=5)
        cb.locator = tick_locator
        tick_font_size = 8
        cb.ax.tick_params(labelsize=tick_font_size)
        cb.update_ticks()

        plt.savefig(f"{self.ImageDir}/3sliceimage_{self.Quantity}_{self.Surface}.png", bbox_inches="tight", dpi=900)
 
    def plot_orbital_phase(self, phase_choice, ax, colormap='viridis', zoom=1, log=False, plot_inset=False):
        """
        Plots a 2D slice of the simulation at a given orbital phase.

        Parameters
        ----------
        phase_choice : float
            The orbital phase at which to plot the slice.

        ax : matplotlib.axes._subplots.AxesSubplot
            The axes on which to plot the slice.

        colormap : str, optional
            The colormap to use for the slice. The default is 'viridis'.

        zoom : float, optional
            The factor by which to zoom in on the slice. The default is 1.

        log : bool, optional
            Whether to plot the slice on a log scale. The default is False.

        plot_inset : bool, optional
            Whether to plot an inset of the slice. The default is False.
        
        Returns
        -------
        None.

        """
        if (self.Surface == 'XY'): xlabel = 'x'; x = 0; ylabel = 'y'; y = 1 ; z = 2 
        if (self.Surface == 'XZ'): xlabel = 'x'; x = 0; ylabel = 'z'; y = 2 ; z = 1
        if (self.Surface == 'YZ'): xlabel = 'y'; x = 1; ylabel = 'z'; y = 2 ; z = 1

        snapshot_list = []  # list of pre-periastron snapshots
        snapshot_phase = [] # list of pre-periastron times

        for k in range(len(self.evolution)):
            data = ReadData(self.evolution[k])
            sim_time = data.sim_time().value
            data.close()
            time = ((sim_time + self.start_time) * units.s).to(units.yr)
            phase = ((self.period + time)/self.period).value

            if round(phase, len(str(phase_choice))-2) == round(phase_choice, len(str(phase_choice))-2):
                snapshot_list.append(k)
                snapshot_phase.append(phase-phase_choice)
            
        index = snapshot_list[np.argmin(np.abs(snapshot_phase))]

        
        data = ReadData(self.evolution[index])
        N_level = data.nlevels()
        N_grids = data.ngrid()
        baseline_data = data.get_3Darray(self.Quantity)
        data.close()
        fluid_parameter = baseline_data['data']
        dims_min = (baseline_data['min_extents'] * units.cm).to(units.astrophys.au)
        dims_max = (baseline_data['max_extents'] * units.cm).to(units.astrophys.au)
        time = ((baseline_data['sim_time'].value + self.start_time) * units.s).to(units.yr)

        ax.set_xlim(dims_min[0][x].value/zoom, dims_max[0][x].value/zoom)
        ax.set_ylim(dims_min[0][y].value/zoom, dims_max[0][y].value/zoom)
        ax.set_xlabel(f"{xlabel} ({str(dims_min.unit)})", fontsize=8)
        ax.set_ylabel(f"{ylabel} ({str(dims_min.unit)})", fontsize=8)

        if plot_inset==True:
            axins = ax.inset_axes([0.6, 0, 0.4, 0.4], transform=ax.transAxes)
            axins.set_xlim(dims_min[0][x].value/(zoom*8), dims_max[0][x].value/(zoom*8))
            axins.set_ylim(dims_min[0][y].value/(zoom*8), dims_max[0][y].value/(zoom*8))
            axins.set_xticks([])
            axins.set_yticks([])
            axins.set_xlabel("")
            axins.set_ylabel("")
            ax.indicate_inset_zoom(axins, edgecolor="black")

        for l in range(N_level):  # plotting each levels

            plot_data = np.array(fluid_parameter)
            
            # Slicer ###############################################################
            slice_location = int(0.5 * N_grids[z])
            if (self.Surface == 'YZ'): sliced_data = plot_data[l][:, :, slice_location - 1]
            if (self.Surface == 'XZ'): sliced_data = plot_data[l][:, slice_location - 1, :]
            if (self.Surface == 'XY'): sliced_data = plot_data[l][slice_location - 1, :, :]
            # Slicer ends #*********************************************************
        
            extents = [dims_min[l][x].value, dims_max[l][x].value, dims_min[l][y].value, dims_max[l][y].value]

            if log == True:
                image = ax.imshow(np.log10(sliced_data), interpolation="nearest", cmap=colormap,
                                extent=extents,
                                origin="lower",
                                vmin=self.Tolerance[0], vmax=self.Tolerance[1]
                                )

                label = r"$\log_{10}$" + f"({self.quantity_dict[self.Quantity]}) "   

                if plot_inset==True:
                    axins.imshow(np.log10(sliced_data), interpolation="nearest", cmap=colormap,
                                extent=extents,
                                origin="lower",
                                vmin=self.Tolerance[0], vmax=self.Tolerance[1]
                                )                

            else:
                image = ax.imshow(sliced_data, interpolation="nearest", cmap=colormap,
                                extent=extents,
                                origin="lower",
                                vmin=self.Tolerance[0], vmax=self.Tolerance[1]
                                )                  

                label = fr"{self.quantity_dict[self.Quantity]} " 
                
                if plot_inset==True:
                    axins.imshow(sliced_data, interpolation="nearest", cmap=colormap,
                                extent=extents,
                                origin="lower",
                                vmin=self.Tolerance[0], vmax=self.Tolerance[1]
                                )    

        if not self.Quantity == 'Tr000_WIND':
            label = label + f"({self.unit_dict[self.Quantity]})" 
        

