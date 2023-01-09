from silo_to_yt import *
import sys
import pickle
yt.set_log_level("ERROR")
import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np
plt.style.use('science')
from matplotlib import ticker
import sys
import pathlib
home = str(pathlib.Path.home())
sys.path.insert(0, os.path.join(home, "code/project/scripts/"))
from binary_system_calculations.binary_velocity import BinarySystem

class YTPlotFunction():

    def __init__(self, data_path, img_dir, quantities, tolerance, period=7.992, start_time=0):

        self.evolution = make_snapshots(data_path) # list of all snapshots
        self.data_path = data_path # path to data
        self.img_dir = img_dir # path to save images
        self.quantities = quantities # list of quantities to include in dataset
        self.tolerance = tolerance 
        self.period = period*u.yr # period of binary system
        self.start_time = start_time*u.s # time of first snapshot wrt periastron

    def add_synchrotron_emission(self, ds):
        
        def _Isync(field, data):
            Isync = data["magnetic_field_magnitude"]**(3/2) * data["pressure"] * data["NG_Mask"]
            return Isync

        # add synchrotron emission field to dataset
        ds.add_field(("gas", "Isync"), function=_Isync, units=None, sampling_type="cell", force_override=True) 

    def three_slice_indices(self, d_phase=0.01): 
        
        pre_periastron_list = []  # list of pre-periastron snapshots
        pre_periastron_phase = [] # list of pre-periastron times
        periastron_list = [] # list of periastron snapshots
        periastron_phase = [] # list of periastron times
        post_periastron_list = [] # list of post-periastron snapshots
        post_periastron_phase = [] # list of post-periastron times

        for k in range(len(self.evolution)):
            data = ReadData(self.evolution[k])
            sim_time = data.sim_time()
            data.close()
            time = ((sim_time + self.start_time)).to(u.yr)
            phase = ((self.period + time)/self.period).value
            phase_str = f"{((self.period + time)/self.period):.2f}"

            pre_per_phase = 1 - d_phase
            post_per_phase = 1 + d_phase

            if phase_str == "1.00":
                periastron_list.append(k)
                periastron_phase.append(abs(phase - 1.0))
            elif phase_str == f"{pre_per_phase}":
                pre_periastron_list.append(k)
                pre_periastron_phase.append(abs(phase - pre_per_phase))
            elif phase_str == f"{post_per_phase}":
                post_periastron_list.append(k)
                post_periastron_phase.append(abs(phase - post_per_phase))
        
        plot_indices = [pre_periastron_list[np.argmin(pre_periastron_phase)], periastron_list[np.argmin(periastron_phase)], post_periastron_list[np.argmin(post_periastron_phase)]]
        return plot_indices

    def plot_synchrotron_emission(self, i, north_vec, norm):

        print(f"Loading dataset: Sim Snapshot {i}")
        ds = get_ds(self.evolution[i], quantities=["magnetic_field", "pressure", "NG_Mask"], start_time=self.start_time)
        time = ds.current_time.value
        print(time)
        print(f"Successfully Loaded dataset: {str(ds)}")
        self.add_synchrotron_emission(ds)
        print(f"Added Synchrotron emission field")
        print(f"Plotting Isync for snapshot {i}...")
        prj = yt.OffAxisProjectionPlot(ds, normal=norm, north_vector=north_vec, fields=("gas", "Isync")) # create projection plot
        prj.set_cmap(("gas", "Isync"), "plasma")
        prj.hide_colorbar("Isync")
        prj.set_figure_size(5)
        prj.set_zlim(("gas", "Isync"), 1e10, 1e16)
        print(f"Saving image Isync_{i}.png...")

        phase = ((self.period + (time*u.s).to(u.yr))/self.period).value

        fig = prj.export_to_mpl_figure((1,1), cbar_mode="None") # export to matplotlib figure
        ax = fig.axes[0] 
        ax.set_xlabel("x (AU)") 
        ax.set_ylabel("y (AU)")
        st = r"$\phi$ = " + f"{phase:.2f}" 
        ax.text(0.1, 0.9, st, color="black", fontsize=8, transform=ax.transAxes, 
                    bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=1', alpha=0.5))
        # fig.savefig(os.path.join(self.img_dir, f"Isync_{i}.png"), dpi=300, bbox_inches="tight")
        return fig 

    def temp_plotter(self, norm, north_vec):
        self.plot_indices = self.three_slice_indices() # list of indices for pre-periastron, periastron, and post-periastron snapshots
        num = 0
        for index in self.plot_indices:
            fig =  self.plot_synchrotron_emission(index, norm, north_vec)
            ax = fig.axes[0]
            fig.savefig(os.path.join(self.img_dir, f"Isync_{index}.png"), dpi=300, bbox_inches="tight")

            num += 1

    def plot_projected_quantity(self, i, plot_quantity, ds_quantities, north_vec, norm, **kwargs):
        print(f"Loading dataset: Sim Snapshot {i}")
        ds = get_ds(self.evolution[i], ds_quantities, start_time=self.start_time)
        time = ds.current_time.value
        print(time)
        print(f"Successfully Loaded dataset: {str(ds)}")
        print(f"Plotting {plot_quantity} for snapshot {i}...")
        prj = yt.OffAxisProjectionPlot(ds, normal=norm, north_vector=north_vec, fields=("gas", plot_quantity))
        prj.set_cmap(("gas", plot_quantity), kwargs.get("cmap", "viridis"))
        prj.set_figure_size(kwargs.get("figsize", 5))
        # prj.set_zlim(("gas", plot_quantity), kwargs.get("zlim", None))
        print(f"Saving image {plot_quantity}_{i}.png...")
        phase = ((self.period + (time*u.s).to(u.yr))/self.period).value

        fig = prj.export_to_mpl_figure((1,1), cbar_mode="None") # export to matplotlib figure
        ax = fig.axes[0]
        ax.set_xlabel("x (AU)")
        ax.set_ylabel("y (AU)")
        st = r"$\phi$ = " + f"{phase:.2f}"
        ax.text(0.1, 0.9, st, color="black", fontsize=8,
                transform=ax.transAxes,
                bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=1', alpha=0.5))
        
        return fig
    
    def temp_quantity_plotter(self, plot_quantity, ds_quantities, north_vec, norm, **kwargs):
        self.plot_indices = self.three_slice_indices() # list of indices for pre-periastron, periastron, and post-periastron snapshots
        num = 0
        for index in self.plot_indices:
            fig = self.plot_projected_quantity(index, plot_quantity, ds_quantities, north_vec, norm, **kwargs)
            fig.savefig(os.path.join(self.img_dir, f"{plot_quantity}_{index}.png"), dpi=300, bbox_inches="tight")
            num += 1

    def threetimeslice(self, plot_quantity, ds_quantities, north_vec, norm, **kwargs):
        self.plot_indices = self.three_slice_indices() # list of indices for pre-periastron, periastron, and post-periastron snapshots

        fig1, axes = plt.subplots(1,3, figsize=(15,5))
        num=0

        for index in self.plot_indices:
            fig = self.plot_projected_quantity(index, plot_quantity, ds_quantities, north_vec, norm, **kwargs)
            axes[num].imshow(fig.axes[0].images[0].get_array(), origin="lower", 
                            extent=fig.axes[0].images[0].get_extent(), cmap="plasma",
                            vmin=1e3, vmax=1e9)
            num += 1
        
        cbaxes = fig.add_axes([axes[-1].get_position().x1+0.01,axes[-1].get_position().y0,0.02,axes[-1].get_position().height])
        cb = fig.colorbar(axes[-1].images[0], ax=axes[-1], orientation="vertical", cax=cbaxes, pad=3.0)
        tick_locator = ticker.MaxNLocator(nbins=5)
        cb.locator = tick_locator
        tick_font_size = 8
        # cb.set_label(label=r"Flux Density (erg cm$^{-2}$ s$^{-1}$ arcsec$^{-2}$)", fontsize=10, labelpad=4)
        cb.formatter.set_powerlimits((0, 0))
        cb.formatter.set_useMathText(True)
        cb.ax.tick_params(labelsize=tick_font_size)
        cb.update_ticks()
        
        fig1.savefig(os.path.join(self.img_dir, f"{plot_quantity}_wr140.png"), dpi=300)

    def plot_smr_plot(self, i):
        print(f"Loading dataset: Sim Snapshot {i}")
        ds = get_ds(self.evolution[i], quantities=["windtracer"], start_time=self.start_time)
        time = ds.current_time.value
        print(time)
        print(f"Successfully Loaded dataset: {str(ds)}")
        print(f"Plotting SMR for snapshot {i}...")
        slc = yt.SlicePlot(ds, "z", "windtracer")
        slc.set_cmap("windtracer", "gray")
        slc.set_figure_size(5)
        slc.zoom(2)
        slc.set_log("windtracer", False)
        slc.set_zlim("windtracer", -10, 0)
        slc.annotate_grids(linewidth=1, edgecolors="black")
        slc.annotate_cell_edges(line_width=0.00001, alpha=0.5, color="black")

        fig = slc.export_to_mpl_figure((1,1), cbar_mode="None") # export to matplotlib figure
        ax = fig.axes[0]
        ax.set_xlabel("x (AU)")
        ax.set_ylabel("y (AU)")

        wr140 = BinarySystem(0.8993, 2895, 10.31, 29.27, "wr140")
        orb_sol = wr140.orb_sol

        x1, y1, x2, y2, vx1, vy1, vx2, vy2 = orb_sol
        x1_au, y1_au = (x1*u.cm).to(u.au).value, (y1*u.cm).to(u.au).value
        x2_au, y2_au = (x2*u.cm).to(u.au).value, (y2*u.cm).to(u.au).value
        ax.plot(x1_au, y1_au, label="WR star", lw=2)
        ax.plot(x2_au, y2_au, label="O star", lw=2)
        ax.legend(loc="upper right", frameon=True, framealpha=1, facecolor="white")
        fig.savefig(os.path.join(self.img_dir, f"smr_{i}.png"), dpi=300, bbox_inches="tight")

    def plot_Bfield(self, i):
        print(f"Loading dataset: Sim Snapshot {i}")
        ds = get_ds(self.evolution[i], quantities=["windtracer", "magnetic_field"], start_time=self.start_time)
        time = ds.current_time.value
        print(time)
        print(f"Successfully Loaded dataset: {str(ds)}")
        print(f"Plotting SMR for snapshot {i}...")
        slc = yt.SlicePlot(ds, "z", "windtracer")
        slc.set_cmap("windtracer", "RdBu")
        slc.set_figure_size(5)
        slc.zoom(2)
        slc.set_log("windtracer", False)
        # slc.set_zlim("windtracer", -10, 0)
        slc.annotate_magnetic_field(normalize=True, cmap="RdBu")

        slc.annotate_quiver("magnetic_field_x", "magnetic_field_y", factor=16, cmap="inferno")


        fig = slc.export_to_mpl_figure((1,1), cbar_mode="None") # export to matplotlib figure
        ax = fig.axes[0]
        ax.set_xlabel("x (AU)")
        ax.set_ylabel("y (AU)")

        fig.savefig(os.path.join(self.img_dir, f"Bfield_contours_{i}.png"), dpi=300, bbox_inches="tight")


    

# def main():
#     # data_path = '/mnt/massive-stars/data/thomas_simulations/wr140-sims/covertex_start/red_z_res/wr140-mhd-l7n256/'
#     data_path = "/mnt/massive-stars/data/thomas_simulations/wr140-sims/compton_cool/compton_mhd_l7n128_v2"
#     # img_dir = "/home/visitor_ap4/code/project/images/sync-emission/"
#     img_dir = "/home/visitor_ap4/code/project/images/"
#     # quantities = ["density", "pressure", "magnetic_field_magnitude"]
#     quantities = ["windtracer"]
#     tolerance = 0.0001
#     # start = -7.25e6
#     start=-2.671e7

#     yt_plot = YTPlotFunction(data_path, img_dir, quantities, tolerance, start_time=start)
#     # plot_indices = yt_plot.three_slice_indices()
#     # print(plot_indices)

#     # yt_plot.plot_synchrotron_emission(plot_indices[1], norm=[-1,0,0], north_vec=[-0.703,0.452,-0.549])
#     # yt_plot.threetimeslice()
#     # yt_plot.temp_plotter(norm=[-1,0,0], north_vec=[-0.703,0.452,-0.549])

#     yt_plot.plot_smr_plot(50)

#     # # yt_plot.temp_quantity_plotter("temperature", ["temperature"], norm=[-1,0,0], north_vec=[0,0,1])
#     # yt_plot.threetimeslice("temperature", ["temperature"], norm=[-1,0,0], north_vec=[0,0,1])

# def main():

#     data_path = "/mnt/massive-stars/data/thomas_simulations/wr140-sims/compton_cool/compton_mhd_l7n128_v2"
#     img_dir = "/home/visitor_ap4/code/project/images/"
#     quantities = ["windtracer"]
#     tolerance = 0.0001
#     start=-2.671e7
#     yt_plot = YTPlotFunction(data_path, img_dir, quantities, tolerance, start_time=start)
#     yt_plot.plot_smr_plot(50)

def main():

    data_path = "/mnt/massive-stars/data/thomas_simulations/wr140-sims/compton_cool/compton_mhd_l7n128_v2"
    img_dir = "/home/visitor_ap4/code/project/images/"
    quantities = ["windtracer", "magnetic_field"]
    tolerance = 0.0001
    start = -7.25e6
    yt_plot = YTPlotFunction(data_path, img_dir, quantities, tolerance, start_time=start)
    yt_plot.plot_smr_plot(50)

if __name__ == "__main__":
    main()