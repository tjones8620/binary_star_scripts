import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
import re
import glob
import argparse
from pypion.ReadData import ReadData
import astropy.units as u

import pathlib
home = str(pathlib.Path.home())
import sys
sys.path.insert(0, os.path.join(home, 'code/project/scripts/'))

from yt_folder.silo_to_yt import make_snapshots

# class ArgparseInputs:
#     def __init__(self):
#         parser = argparse.ArgumentParser()
#         parser.add_argument("Path", help="Path to data folder")
#         parser.add_argument("fits_dir", help="Path to fits directory")
#         parser.add_argument("png_dir", help="Path to png directory")
#         parser.add_argument('--tol', nargs='+', type=float, default=[np.log10(1.0e-18), np.log10(2.0e-13)])
#         parser.add_argument('--cmap', default='viridis')
#         args = parser.parse_args()

#         self.path = args.Path
#         self.tolerance = args.tol
#         self.img_dir = args.img_dir





# Defining function to create xray projection images

def create_xray_proj_img(data_path, fits_dir, png_dir, start_time, **kwargs):
    start_time = start_time*u.s

    evolution = make_snapshots(data_path) # Creating snapshots of data

    fits_dir = os.path.join(data_path, fits_dir) # Path to directory containing fits files
    png_dir = os.path.join(home, png_dir) # Path from home dir to directory where png files will be saved

    if not os.path.exists(png_dir):
        os.makedirs(png_dir)
    else:
        print("Directory already exists")
        # Clean directory
        files = glob.glob(os.path.join(png_dir, "*.png"))
        for f in files:
            os.remove(f)
            

    # Getting paths to all fits files in fits_dir
    fits_files = sorted(os.listdir(fits_dir))
    fits_files.sort(key=lambda test_string : list(
        map(int, re.findall(r'\d+', test_string))))

    # Function returning index of snapshots within a given time range of periastron
    def snapshots_before_after_periastron(evolution, tolerance, start_time):
        sim_time_array = []
        for file in evolution:
            data = ReadData(file)
            sim_time = data.sim_time().to(u.yr) + start_time.to(u.yr)
            data.close()
            sim_time_array.append(np.abs(sim_time.value))
        
        sim_time_array = np.array(sim_time_array)
        args = np.where(sim_time_array < tolerance)
        return list(args[0])
    
    # Getting indices of snapshots within a given time range of periastron
    args = snapshots_before_after_periastron(evolution, kwargs.get("timerange", 0.01), start_time)

    # Creating xray projection images
    for i in args:
        file = fits_files[i]
        data = ReadData(evolution[i])
        sim_time = data.sim_time().to(u.yr) + start_time.to(u.yr)

        data.close()
        print(f"Processing {file} at {sim_time:3e}")

        hdul = fits.open(os.path.join(fits_dir, file))
        data = hdul[9].data - hdul[11].data
        hdul.close()

        fig, ax = plt.subplots(figsize=(10,10))
        image = ax.imshow(data, cmap="viridis", origin="lower", vmin=0, vmax=0.00225)
        fig.colorbar

        cb = fig.colorbar(image, ax=ax, orientation="vertical", pad=0.1)

        if sim_time < 0:
            ax.text(0.05, 0.95, f"t = {abs(sim_time):.4f} before periastron", 
            bbox=dict(facecolor='white', alpha=0.5), transform=ax.transAxes)
        else:
            ax.text(0.05, 0.95, f"t = {sim_time:.4f} after periastron", 
            bbox=dict(facecolor='white', alpha=0.5), transform=ax.transAxes)

        plt.tick_params(axis='both', bottom=False, left=False, labelbottom=False, labelleft=False)
        plt.savefig(os.path.join(png_dir, file[:-5] + ".png"), bbox_inches="tight", pad_inches=0)
        plt.close(fig)
    
def main():
    # Defining time before periastron where simulation was started
    start_time = -2.679e7 * u.s
    print(start_time.to(u.yr))

    # Path to directory containing silo files that were used to create the fits files
    data_dir = "/mnt/massive-stars/data/thomas_simulations/wr140-sims/covertex_start/red_z_res/wr140-mhd-l7n256/"

    # Path to directory containing fits files
    fits_dir = os.path.join(data_dir, "proj/XY_proj2")

    # Path to directory where png files will be saved
    png_dir = os.path.join(home, "code/project/images/xray-proj/wr140-mhd-n256/XY/")


    create_xray_proj_img(data_dir, "proj/XY_proj2", "code/project/images/xray-proj/wr140-mhd-n256/XY/", start_time.value, timerange=0.04)


if __name__ == "__main__":
    main()
