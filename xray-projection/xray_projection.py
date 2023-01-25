__author__ = "Thomas Jones"

import numpy as np
import matplotlib.pyplot as plt
plt.style.use('science')
from astropy.io import fits
import os
import re
import glob
import argparse
from pypion.ReadData import ReadData
import astropy.units as u
import pandas as pd
import sigfig
from matplotlib import ticker
from palettable.mycarta import LinearL_5
from palettable.cmocean.sequential import Ice_12
from math import log10, floor



import pathlib
home = pathlib.Path.home()

import sys
sys.path.insert(0, os.path.join(home, 'code/project/scripts/'))

from yt_folder.silo_to_yt import make_snapshots


class XrayProj:
    def __init__(self, data_path, fits_dir, png_dir, start_time, period, name = "binary", **kwargs):
        self.start_time = start_time*u.s
        self.evolution = make_snapshots(data_path)
        self.fits_dir = fits_dir
        self.png_dir = png_dir
        self.period = period*u.yr
        self.kwargs = kwargs
        self.plane = kwargs.get('plane', 'xy')
        self.name = name

        self.trajectory_file = kwargs.get('trajectory_file', None)
        self.create_png_dir()

        # Getting paths to all fits files in fits_dir
        self.fits_files = sorted(os.listdir(fits_dir))
        self.fits_files.sort(key=lambda test_string : list(
            map(int, re.findall(r'\d+', test_string))))

        self.get_extents()

    def create_png_dir(self):
        """
        Creates a directory to store the png files in
        """
        if not os.path.exists(self.png_dir):
            os.makedirs(self.png_dir)
        else:
            print(f"Directory {self.png_dir} already exists")

    @staticmethod
    def get_star_position(trajectory_file, times):
        """
        Returns the x, y, z positions of the stars at the given times
        from the trajectory file saved by the simulation code.

        Parameters
        ----------
        trajectory_file : str
            Path to the trajectory file
        times : list
            List of times to get the positions of the stars at

        Returns
        -------
        star1_x : list
            List of x positions of star 1
        star1_y : list
            List of y positions of star 1
        star2_x : list
            List of x positions of star 2
        star2_y : list  
            List of y positions of star 2
        star1_z : list
            List of z positions of star 1
        star2_z : list
            List of z positions of star 2
        """

        def find_neighbours(value, df, colname):
            exactmatch = df[df[colname] == value]
            if not exactmatch.empty:
                return exactmatch.index
            else:
                lowerneighbour_ind = df[df[colname] < value][colname].idxmax()
                upperneighbour_ind = df[df[colname] > value][colname].idxmin()
                return lowerneighbour_ind

        headers = ['time', 'unknown','star1_x', 'star1_y', 'star1_z', 
                'star1_vx', 'star1_vy', 'star1_vz', 'star2_x', 'star2_y', 
                'star2_z', 'star2_vx', 'star2_vy', 'star2_vz']

        times = [sigfig.round(time, sigfigs=7) for time in times]

        df = pd.read_csv(trajectory_file, delim_whitespace=True, names=headers)

        index_closest_time = find_neighbours(times[0], df, 'time')
        df = df.loc[index_closest_time]

        star1_x = (df['star1_x'].tolist()*u.cm).to(u.au).value
        star1_y = (df['star1_y'].tolist()*u.cm).to(u.au).value
        star2_x = (df['star2_x'].tolist()*u.cm).to(u.au).value
        star2_y = (df['star2_y'].tolist()*u.cm).to(u.au).value
        star1_z = (df['star1_z'].tolist()*u.cm).to(u.au).value
        star2_z = (df['star2_z'].tolist()*u.cm).to(u.au).value

        return star1_x, star1_y, star2_x, star2_y, star1_z, star2_z

    def get_extents(self):
        """
        Gets the extents of the simulation box from the first file in the
        directory with the simulation SILO files
        """
        baseline = ReadData(self.evolution[0]).get_3Darray('Density')
        min_extents = baseline['min_extents'][1]
        max_extents = baseline['max_extents'][1]
        del baseline

        x_min = ((min_extents[0]*u.cm).to(u.au).value)/(2**2)
        x_max = ((max_extents[0]*u.cm).to(u.au).value)/(2**2)
        y_min = ((min_extents[1]*u.cm).to(u.au).value)/(2**2)
        y_max = ((max_extents[1]*u.cm).to(u.au).value)/(2**2)
        z_min = ((min_extents[2]*u.cm).to(u.au).value)/(2**2)
        z_max = ((max_extents[2]*u.cm).to(u.au).value)/(2**2)

        if self.plane == 'xy':
            self.extents = [x_min, x_max, y_min, y_max]
        elif self.plane == 'xz':
            self.extents = [z_min, z_max, x_min, x_max]

    def create_imgs(self, i, ax, band="hard", fig=None):
        """
        Creates X-ray images for a given time simulation
        snapshot

        Parameters
        ----------
        i : int
            The index of the simulation snapshot to be processed
        ax : matplotlib.axes._subplots.AxesSubplot
            The axes to plot the image on
        band : str, optional
            The X-ray band to plot, by default "hard"
        fig : matplotlib.figure.Figure, optional
            The figure to plot the image on, by default None
        
        Returns
        -------
        matplotlib.image.AxesImage
            The image of the X-ray data
        """

        file = self.fits_files[i]
        data = ReadData(self.evolution[i])
        sim_time = data.sim_time().value
        time = ((sim_time + self.start_time.value) * u.s).to(u.yr)
        phase = ((self.period + time)/self.period).value

        data.close()
        print(f"Processing {file} at orbital phase {phase:.3f}")

        hdul = fits.open(os.path.join(self.fits_dir, file))
        hard_xray_data = hdul[9].data - hdul[11].data
        soft_xray_data = hdul[4].data - hdul[9].data
        if  band == "hard":
            data = hard_xray_data
        elif band == "soft":
            data = soft_xray_data
        hdul.close()

        order_of_mag = floor(log10(np.max(data)))
        image = ax.imshow(data/10**(order_of_mag), cmap=self.kwargs.get("cmap", "gist_heat"), origin="lower", 
                            vmin=self.kwargs.get("vmin", 0), vmax=self.kwargs.get("vmax", 0.00225), 
                            extent=self.extents)
        

        zoom_factor = self.kwargs.get('zoom', 1)
        ax.set_xlim(self.extents[0]/zoom_factor, self.extents[1]/zoom_factor)
        ax.set_ylim(self.extents[2]/zoom_factor, self.extents[3]/zoom_factor)

        if self.trajectory_file:
            star1_x, star1_y, star2_x, star2_y, star1_z, star2_z = self.get_star_position(self.trajectory_file, [sim_time])

            if self.plane == "xy":
                x1, y1 = star1_x*-1, star1_y
                x2, y2 = star2_x*-1, star2_y
                xlabel = "x (AU)"
                ylabel = "y (AU)"

            if self.plane == "xz":
                x1, y1 = star1_z*-1, star1_x
                x2, y2 = star2_z*-1, star2_x
                xlabel = "z (AU)"
                ylabel = "x (AU)"

            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.scatter(x1, y1, marker="*", color="white", s=150)
            ax.scatter(x2, y2, marker="*", color="green", s=150)
        
        else:
            xlabel = "x (AU)"
            ylabel = "y (AU)"
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)

        if fig:
            cbaxes = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.03, ax.get_position().width-0.02])
            cb = fig.colorbar(image, ax=ax, orientation="vertical", cax=cbaxes, pad=10.0)
            tick_locator = ticker.MaxNLocator(nbins=5)
            cb.locator = tick_locator
            tick_font_size = 8
            cb.set_label(label="Flux Density ($10^{" + f"{order_of_mag}" +"}$ erg cm$^{-2}$ s$^{-1}$ arcsec$^{-2}$)", fontsize=18, labelpad=10, )
            # Option to put label above colorbar
            cb.ax.xaxis.set_label_position('top')
            cb.ax.xaxis.set_ticks_position('top')
            cb.formatter.set_powerlimits((0, 0))
            cb.formatter.set_useMathText(True)
            cb.ax.tick_params(labelsize=tick_font_size)
            cb.update_ticks()
            fig.patch.set_facecolor('white')
            fig.patch.set_alpha(0.0)
            fig.savefig(os.path.join(home, "code/project/scripts/images/xray-proj/wr140-compton-mhd-n256/XY/soft/test.png"), pad_inches=0.1, bbox_inches="tight", dpi=900, transparent=True)

        return image

    def three_time_slice(self, colormap='viridis', d_phase = 0.005):

        """
        Creates a 3x1 plot of the x-ray image at the pre-periastron, periastron, and post-periastron snapshots.
        The snapshots are chosen to be the closest to the pre-periastron, periastron, and post-periastron times.
        
        Parameters
        ----------
        colormap : str
            The colormap to use for the image.
        d_phase : float
            The phase difference between the pre-periastron, periastron, and post-periastron snapshots.
        
        Returns
        -------
        fig : matplotlib.pyplot.figure
            The figure containing the 3x1 plot of the x-ray image at the pre-periastron, periastron, and post-periastron snapshots.
        axes : matplotlib.pyplot.axes
            The axes of the figure.
        """
        
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
            time = ((sim_time + self.start_time.value) * u.s).to(u.yr)
            phase = ((self.period + time)/self.period).value

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

        fig, axes = plt.subplots(1, 3, figsize=(12, 4), sharey=True)
        num = 0

        for index in plot_indices:
            # print(index)
            image = self.create_imgs(index, ax=axes[num], band="hard")
            num += 1

        axes[1].set_ylabel(None)
        axes[2].set_ylabel(None)

        cbaxes = fig.add_axes([axes[num-1].get_position().x1+0.01,axes[num-1].get_position().y0,0.02,axes[num-1].get_position().height])
        cb = fig.colorbar(image, ax=axes[num-1], orientation="vertical", cax=cbaxes, pad=3.0)
        tick_locator = ticker.MaxNLocator(nbins=5)
        cb.locator = tick_locator
        tick_font_size = 8
        cb.set_label(label=r"Flux Density (erg cm$^{-2}$ s$^{-1}$ arcsec$^{-2}$)", fontsize=10, labelpad=4)
        cb.formatter.set_powerlimits((0, 0))
        cb.formatter.set_useMathText(True)
        cb.ax.tick_params(labelsize=tick_font_size)
        cb.update_ticks()

        print(f"Saving 3-time plot to {os.path.join(self.png_dir, 'three_time_slice.png')} ...")
        plt.savefig(os.path.join(self.png_dir, "three_time_slice.png"), pad_inches=0.1, bbox_inches="tight", dpi=300)

    def five_time_slice(self, colormap='viridis', d_phase = 0.005, band="hard"):
        """
        Save a 5-time slice plot of the x-ray emission, with the phases marked.

        Parameters
        ----------
        colormap : str, optional
            The colormap to use for the image. Default is 'viridis'.
        d_phase : float, optional
            The phase difference between each time slice. Default is 0.005.
        band : str, optional
            The band to use for the image. Default is 'hard'.
        
        Returns
        -------
        None
        """
        
        phase1_list = [] ;  phase1_phase = [] 
        phase2_list = [] ;  phase2_phase = []
        phase3_list = [] ;  phase3_phase = [] 
        phase4_list = [] ;  phase4_phase = []
        phase5_list = [] ;  phase5_phase = [] 

        for k in range(len(self.evolution)):
            data = ReadData(self.evolution[k])
            sim_time = data.sim_time().value
            data.close()
            time = ((sim_time + self.start_time.value) * u.s).to(u.yr)
            phase = ((self.period + time)/self.period).value

            pre_per_phase1 = 1 - d_phase
            pre_per_phase2 = 1 - d_phase/2
            post_per_phase1 = 1 + d_phase
            post_per_phase2 = 1 + d_phase/2

            if round(phase, len(str(d_phase))-2) == round(pre_per_phase1, len(str(d_phase))-2):
                phase1_list.append(k)
                phase1_phase.append(abs(phase - pre_per_phase1))
            elif round(phase, len(str(d_phase))-2) == round(pre_per_phase2, len(str(d_phase))-2):
                phase2_list.append(k)
                phase2_phase.append(abs(phase - pre_per_phase2))
            elif round(phase, len(str(d_phase))-2) == round(1, len(str(d_phase))-2):
                phase3_list.append(k)
                phase3_phase.append(abs(phase - 1.0))
            elif round(phase, len(str(d_phase))-2) == round(post_per_phase2, len(str(d_phase))-2):
                phase4_list.append(k)
                phase4_phase.append(abs(phase - post_per_phase2))
            elif round(phase, len(str(d_phase))-2) == round(post_per_phase1, len(str(d_phase))-2):
                phase5_list.append(k)
                phase5_phase.append(abs(phase - post_per_phase1))

        try:
            plot_indices = [phase1_list[np.argmin(phase1_phase)], phase2_list[np.argmin(phase2_phase)], phase3_list[np.argmin(phase3_phase)], phase4_list[np.argmin(phase4_phase)], phase5_list[np.argmin(phase5_phase)]]

        except:
            if ValueError:
                print("The time slice is not found in the evolution data. Please check the d_phase value.")
                return None

        fig, axes = plt.subplots(1, 5, figsize=(12*(5/3), 4*(5/3)), sharey=True)
        num = 0

        for index in plot_indices:
            image = self.create_imgs(index, ax=axes[num], band=band)
            axes[num].tick_params(axis='both', which='major', labelsize=8*5/3)
            axes[num].set_xlabel(r"x (AU)", fontsize=8*5/3)
            if num == 0:
                axes[num].set_ylabel(r"y (AU)", fontsize=8*5/3)

            else:
                axes[num].set_ylabel(None)
            num += 1
            

        plt.subplots_adjust(wspace=0)

        cbaxes = fig.add_axes([axes[num-1].get_position().x1+0.01,axes[num-1].get_position().y0,0.015,axes[num-1].get_position().height])
        # cbaxes = fig.add_axes([axes[0].get_position().x0, axes[0].get_position().y1+0.01,5*(axes[num-1].get_position().width), 0.02])
        cb = fig.colorbar(image, ax=axes[0], orientation="vertical", cax=cbaxes, pad=3.0)
        # cb.ax.xaxis.set_label_position('top')
        tick_locator = ticker.MaxNLocator(nbins=5)
        cb.locator = tick_locator
        tick_font_size = 8
        cb.set_label(label=r"Flux Density ($erg$ $cm^{-2}$ $s^{-1}$ $arcsec^{-2}$)", fontsize=14, labelpad=10)
        cb.formatter.set_powerlimits((0, 0))
        cb.formatter.set_useMathText(True)
        cb.ax.tick_params(labelsize=tick_font_size)
        cb.update_ticks()

        print(f"Saving 5-time plot to {os.path.join(self.png_dir,  f'{self.name}_five_time.png')} ...")
        plt.savefig(os.path.join(self.png_dir, f"{self.name}_five_time.png"), pad_inches=0.1, bbox_inches="tight", dpi=300)




# def main(): 
#     # Defining time before periastron where simulation was started
#     # start_time = -2.679e7
#     start_time = -2.671e7 # start time for simulations starting at covertex
#     # start_time = -7.25e6 # start time for simulations starting at orbital phase 0.98


#     # Path to directory containing silo files that were used to create the fits files
#     data_dir = "/mnt/massive-stars/data/thomas_simulations/wr140-sims/covertex_start/red_z_res/wr140-mhd-l7n256/"
#     # data_dir = "/mnt/massive-stars/data/thomas_simulations/wr140-sims/late_start/wr140-hydro-n256"

#     # Path to directory containing fits files
#     fits_dir = os.path.join(data_dir, "proj/XY_proj2")
#     # fits_dir = os.path.join(data_dir, "proj/proj_XZ")

#     # Path to directory where png files will be saved
#     png_dir = os.path.join(home, "code/project/scripts/images/xray-proj/wr140-mhd-n256/XY/soft/")
#     # png_dir = os.path.join(home, "code/project/images/xray-proj/wr140-hydro-n256/XZ/")

#     # Path to file containing trajectory data
#     trajectory_file = os.path.join(data_dir, "trajectory.txt")

#     plot = XrayProj(data_dir, fits_dir, png_dir, start_time, trajectory_file=trajectory_file, vmin=0, vmax=0.002, cmap=LinearL_5.mpl_colormap, tolerance=0.03, plane="xy", zoom = 6, period=7.926, name = "mhd_n256_hardxrayproj")
#     plot.five_time_slice(d_phase=0.003, band="hard")

def main(): 
    # Defining time before periastron where simulation was started
    start_time = -1.23e7 # start time for simulations starting at covertex

    # Path to directory containing silo files that were used to create the fits files
    data_dir = "/mnt/massive-stars/data/thomas_simulations/wr140-sims/compton_tests/compton_mhd_n256/"

    # Path to directory containing fits files
    fits_dir = os.path.join(data_dir, "proj/XY")

    # Path to directory where png files will be saved
    png_dir = os.path.join(home, "code/project/scripts/images/xray-proj/wr140-compton-mhd-n256/XY/soft/")

    # Path to file containing trajectory data
    trajectory_file = "/home/visitor_ap4/code/project/scripts/binary_system_calculations/trajectory_phasept95.txt"

    plot = XrayProj(data_dir, fits_dir, png_dir, start_time, trajectory_file=trajectory_file, vmin=0, vmax=0.002, 
                    cmap=LinearL_5.mpl_colormap, plane="xy", zoom = 6, period=7.926, name="compton_mhd_n256_hardxrayproj")
    plot.five_time_slice(d_phase=0.003, band="hard")

def main(): 
    # Defining time before periastron where simulation was started
    start_time = -1.23e7 # start time for simulations starting at covertex

    # Path to directory containing silo files that were used to create the fits files
    data_dir = "/mnt/massive-stars/data/thomas_simulations/wr140-sims/compton_tests/compton_mhd_n256/"

    # Path to directory containing fits files
    fits_dir = os.path.join(data_dir, "proj/XY")

    # Path to directory where png files will be saved
    png_dir = os.path.join(home, "code/project/scripts/images/xray-proj/wr140-compton-mhd-n256/XY/soft/")

    # Path to file containing trajectory data
    trajectory_file = "/home/visitor_ap4/code/project/scripts/binary_system_calculations/trajectory_phasept95.txt"

    fig, ax = plt.subplots(figsize=(6, 6))

    plot = XrayProj(data_dir, fits_dir, png_dir, start_time, trajectory_file=trajectory_file, vmin=0, vmax=2, 
                    cmap=LinearL_5.mpl_colormap, plane="xy", zoom = 6, period=7.926, name="compton_mhd_n256_hardxrayproj")
    fig = plot.create_imgs(ax=ax, i=115, band="hard", fig=fig)
    
    ax.set_ylabel("")
    ax.set_xlabel("")
    ax.set_xticks([])
    ax.set_yticks([])

if __name__ == "__main__":
    main()
