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

import pathlib
home = pathlib.Path.home()

import sys
sys.path.insert(0, os.path.join(home, 'code/project/scripts/'))

from yt_folder.silo_to_yt import make_snapshots


class XrayProj:
    def __init__(self, data_path, fits_dir, png_dir, start_time, trajectory_file, **kwargs):
        self.start_time = start_time*u.s
        self.evolution = make_snapshots(data_path)
        self.fits_dir = fits_dir
        self.png_dir = png_dir
        self.kwargs = kwargs
        self.plane = kwargs.get('plane', 'xy')
        self.create_png_dir()

        # Getting paths to all fits files in fits_dir
        self.fits_files = sorted(os.listdir(fits_dir))
        self.fits_files.sort(key=lambda test_string : list(
            map(int, re.findall(r'\d+', test_string))))

        self.get_extents()

        self.args, self.arg_times = self.snapshots_before_after_periastron(self.evolution, self.kwargs.get("tolerance", 0.03), self.start_time)

        # self.arg_times = np.round(self.arg_times, 6)
        self.arg_times = [sigfig.round(time, sigfigs=7) for time in self.arg_times]

        self.get_star_position(trajectory_file)

    def create_png_dir(self):
        if not os.path.exists(self.png_dir):
            os.makedirs(self.png_dir)
        else:
            print("Directory already exists")
            input("Press enter to continue and clear directory")
            # Clean directory
            files = glob.glob(os.path.join(self.png_dir, "*.png"))
            for f in files:
                os.remove(f)

    # Function returning index of snapshots within a given time range of periastron
    @staticmethod
    def snapshots_before_after_periastron(evolution, tolerance, start_time):
        sim_time_array = []
        sim_time_sec_array = []
        for file in evolution:
            data = ReadData(file)

            sim_time_sec = data.sim_time().to(u.s).value
            sim_time_sec_array.append(sim_time_sec)

            sim_time = data.sim_time().to(u.yr) + start_time.to(u.yr)
            data.close()
            sim_time_array.append(np.abs(sim_time.value))
        
        sim_time_array = np.array(sim_time_array)
        args = np.where(sim_time_array < tolerance)
        args = list(args[0])

        sim_time_sec_array = np.array(sim_time_sec_array)

        return args, sim_time_sec_array[args]

    def get_star_position(self, trajectory_file):

        headers = ['time', 'unknown','star1_x', 'star1_y', 'star1_z', 
                'star1_vx', 'star1_vy', 'star1_vz', 'star2_x', 'star2_y', 
                'star2_z', 'star2_vx', 'star2_vy', 'star2_vz']

        df = pd.read_csv(trajectory_file, delim_whitespace=True, names=headers)
        times = self.arg_times
        df = df.loc[df['time'].isin(times)]

        self.star1_x = (df['star1_x'].tolist()*u.cm).to(u.au).value
        self.star1_y = (df['star1_y'].tolist()*u.cm).to(u.au).value
        self.star2_x = (df['star2_x'].tolist()*u.cm).to(u.au).value
        self.star2_y = (df['star2_y'].tolist()*u.cm).to(u.au).value
        self.star1_z = (df['star1_z'].tolist()*u.cm).to(u.au).value
        self.star2_z = (df['star2_z'].tolist()*u.cm).to(u.au).value

    def get_extents(self):
        # Getting extents of the simulation box
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

    def create_imgs(self, i, band="hard"):
        # Creating xray projection images
        file = self.fits_files[i]
        data = ReadData(self.evolution[i])
        sim_time = data.sim_time().to(u.yr) + self.start_time.to(u.yr)

        data.close()
        print(f"Processing {file} at {sim_time:3e}")

        hdul = fits.open(os.path.join(self.fits_dir, file))

        hard_xray_data = hdul[9].data - hdul[11].data
        soft_xray_data = hdul[4].data - hdul[9].data

        if  band == "hard":
            data = hard_xray_data
        elif band == "soft":
            data = soft_xray_data
        # data = kwargs.get('data', hard_xray_data)
        hdul.close()

        fig, ax = plt.subplots(figsize=(10,10))
        image = ax.imshow(data, cmap=self.kwargs.get("cmap", "inferno"), origin="lower", vmin=self.kwargs.get("vmin", 0), vmax=self.kwargs.get("vmax", 0.00225), extent=self.extents)
        
        cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
        
        fig.colorbar(image, ax=ax, cax=cax,label="Flux Density (erg cm$^{-2}$ s$^{-1}$ arcsec$^{-2}$)")

        zoom_factor = self.kwargs.get('zoom', 1)
        ax.set_xlim(self.extents[0]/zoom_factor, self.extents[1]/zoom_factor)
        ax.set_ylim(self.extents[2]/zoom_factor, self.extents[3]/zoom_factor)

        if self.plane == "xy":
            x1, y1 = self.star1_x[i-self.args[0]]*-1, self.star1_y[i-self.args[0]]
            x2, y2 = self.star2_x[i-self.args[0]]*-1, self.star2_y[i-self.args[0]]
            xlabel = "x [AU]"
            ylabel = "y [AU]"

        if self.plane == "xz":
            x1, y1 = self.star1_z[i-self.args[0]]*-1, self.star1_x[i-self.args[0]]
            x2, y2 = self.star2_z[i-self.args[0]]*-1, self.star2_x[i-self.args[0]]
            xlabel = "z [AU]"
            ylabel = "x [AU]"

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        ax.scatter(x1, y1, marker="*", color="yellow", s=100, label="O-Star")
        ax.scatter(x2, y2, marker="*", color="blue", s=50, label="WR-Star")
        ax.legend(loc='upper right', framealpha=0.5, frameon=True)

        if sim_time < 0:
            ax.text(0.05, 0.95, f"{abs(sim_time):.3f} before periastron", 
            bbox=dict(facecolor='white', alpha=0.5), transform=ax.transAxes)
        else:
            ax.text(0.05, 0.95, f"{sim_time:.3f} after periastron", 
            bbox=dict(facecolor='white', alpha=0.5), transform=ax.transAxes)

        plt.savefig(os.path.join(self.png_dir, file[:-5] + ".png"), pad_inches=0.1, bbox_inches="tight")
        plt.close(fig)



def main(): 
    # Defining time before periastron where simulation was started
    # start_time = -2.679e7
    # start_time = -2.671e7 # start time for simulations starting at covertex
    start_time = -7.25e6 # start time for simulations starting at orbital phase 0.98


    # Path to directory containing silo files that were used to create the fits files
    # data_dir = "/mnt/massive-stars/data/thomas_simulations/wr140-sims/covertex_start/red_z_res/wr140-mhd-l7n256/"
    data_dir = "/mnt/massive-stars/data/thomas_simulations/wr140-sims/late_start/wr140-hydro-n256"

    # Path to directory containing fits files
    # fits_dir = os.path.join(data_dir, "proj/XY_proj2")
    fits_dir = os.path.join(data_dir, "proj/proj_XZ")

    # Path to directory where png files will be saved
    # png_dir = os.path.join(home, "code/project/images/xray-proj/wr140-mhd-n256/XY/")
    # png_dir = os.path.join(home, "code/project/images/xray-proj/wr140-mhd-n256/XY/")
    png_dir = os.path.join(home, "code/project/images/xray-proj/wr140-hydro-n256/XZ/")



    # Path to file containing trajectory data
    trajectory_file = os.path.join(data_dir, "trajectory.txt")

    plot = XrayProj(data_dir, fits_dir, png_dir, start_time, trajectory_file=trajectory_file, vmin=0, vmax=0.00225, cmap="inferno", tolerance=0.03, plane="xz", zoom = 4)

    
    for i in plot.args:
        plot.create_imgs(i, band="hard")

if __name__ == "__main__":
    main()
