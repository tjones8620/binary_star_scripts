from pypion_plotter import Plot_Functions
import numpy as np
import matplotlib.pyplot as plt
import os
from palettable.cmocean.sequential import Ice_10
import pathlib
plt.style.use('science')
home = str(pathlib.Path.home())

def main():

    fluidquantity = "Tr000_WIND"
    surface = "XY"
    tolerance = [0, 1]
    img_dir = os.path.join(home, "code/project/scripts/images/SimulationPlots/windtracer/")

    hd1 = {
        "path":"/mnt/massive-stars/data/thomas_simulations/wr140-sims/late_start/wr140-hydro-n256",
        "start_time":-7.25e6,
        "label": "hd-1",
        "color": None,

    }

    hd2 = {
        "path":"/mnt/massive-stars/data/thomas_simulations/wr140-sims/accel_wind/wr140-hydro-l7n256-accel",
        "start_time":-7.25e6,
        "label": "hd-2",
        "color": None,
    }

    mhd1 = {
        "path":"/mnt/massive-stars/data/thomas_simulations/wr140-sims/covertex_start/red_z_res/wr140-mhd-l7n256",
        "start_time":-2.671e7,
        "label": "mhd-1",
        "color": None,

    }

    mhd2 = {
        "path":"/mnt/massive-stars/data/thomas_simulations/wr140-sims/compton_tests/compton_mhd_n256",
        "start_time":-1.239e7,
        "label": "mhd-2",
    }

    sim_list = [hd1, hd2, mhd1, mhd2]

    fig, ax = plt.subplots(ncols=4, figsize=(10, 6), sharey=True)
    i = 0
    for sim in sim_list:
        plot = Plot_Functions(sim['path'], img_dir, fluidquantity, tolerance, surface, period=7.992, start_time=sim['start_time'])
        plot.plot_orbital_phase(ax=ax[i], log=False, colormap="RdBu", phase_choice=1.002, zoom=4, plot_inset=False)
        ax[i].text(0.7, 0.95, sim['label'], transform=ax[i].transAxes, va='top', ha='left', fontsize=14, color="white", bbox=dict(facecolor='black', alpha=0.5))

        if not i==0:
            ax[i].set_ylabel("")

        i += 1
        
    fig.subplots_adjust(hspace=0.01, wspace=0.05)
    plt.savefig(os.path.join(img_dir, f"orbital_phase_windtracer2-combined.png"), dpi=300)

if __name__=="__main__":
    main()
