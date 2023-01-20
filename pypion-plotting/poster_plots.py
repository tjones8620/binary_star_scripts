from pypion_plotter import Plot_Functions
import numpy as np
import matplotlib.pyplot as plt
import os
from palettable.cmocean.sequential import Ice_10
import pathlib
plt.style.use('science')
home = str(pathlib.Path.home())

def main():

    fluidquantity = "Density"
    surface = "XY"
    tolerance = [-20,-13]
    img_dir = os.path.join(home, "code/project/scripts/images/SimulationPlots/windtracer/")

    mhd1 = {
        "path":"/mnt/massive-stars/data/thomas_simulations/wr140-sims/covertex_start/red_z_res/wr140-mhd-l7n256",
        "start_time":-2.671e7,
        "label": "mhd-1",
        "color": None,
        "inset": False

    }

    mhd2 = {
        "path":"/mnt/massive-stars/data/thomas_simulations/wr140-sims/compton_tests/compton_mhd_n256",
        "start_time":-1.239e7,
        "label": "mhd-2",
        "inset": True
    }

    sim_list = [mhd1, mhd2]


    for sim in sim_list:
        fig, ax = plt.subplots(1, 1, figsize=(6, 6))
        plot = Plot_Functions(sim['path'], img_dir, fluidquantity, tolerance, surface, period=7.992, start_time=sim['start_time'])
        plot.plot_orbital_phase(ax=ax, log=True, colormap="viridis", phase_choice=1.001, zoom=1, plot_inset=sim['inset'])

        plt.savefig(os.path.join(img_dir, "density_{}_orbital_phase.png".format(sim['label'])), dpi=900, bbox_inches="tight")

if __name__ == "__main__":
    main()