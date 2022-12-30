from silo_to_yt import *
import sys
import pickle
yt.set_log_level("ERROR")
import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np
plt.style.use('science')



class YTPlotFunction():

    def __init__(self, data_path, img_dir, quantities, tolerance, period=7.992, start_time=-2.671e7):

        self.evolution = make_snapshots(data_path)
        self.data_path = data_path
        self.img_dir = img_dir
        self.quantities = quantities
        self.tolerance = tolerance
        self.period = period*u.yr
        self.start_time = start_time*u.s

        self.plot_indices = self.three_slice_indices()

    
    def add_synchrotron_emission(self, ds):
        def _Isync(field, data):
            Isync = data["magnetic_field_magnitude"]**(3/2) * data["pressure"] * data["NG_Mask"]
            return Isync

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
        prj = yt.OffAxisProjectionPlot(ds, normal=norm, north_vector=north_vec, fields=("gas", "Isync"))
        prj.set_cmap(("gas", "Isync"), "plasma")
        prj.hide_colorbar("Isync")
        prj.set_figure_size(5)
        prj.set_zlim(("gas", "Isync"), 1e10, 1e16)
        print(f"Saving image Isync_{i}.png...")

        phase = ((self.period + (time*u.s).to(u.yr))/self.period).value

        fig = prj.export_to_mpl_figure((1,1), cbar_mode="None")
        ax = fig.axes[0]
        ax.set_xlabel("x (AU)")
        ax.set_ylabel("y (AU)")
        st = r"$\phi$ = " + f"{phase:.2f}" 
        ax.text(0.1, 0.9, st, color="black", fontsize=8, transform=ax.transAxes, 
                    bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=1', alpha=0.5))
        # fig.savefig(os.path.join(self.img_dir, f"Isync_{i}.png"), dpi=300, bbox_inches="tight")
        return fig


    def temp_plotter(self):
        num = 0
        for index in self.plot_indices:
            fig =  self.plot_synchrotron_emission(index, norm=[-1,0,0], north_vec=[-0.703,0.452,-0.549])
            ax = fig.axes[0]
            
            # if num == 0:
            #     pass
            # elif num == 1:
            #     ax.set_ylabel("")
            #     ax.set_yticklabels([])
            # elif num == 2:
            #     ax.set_ylabel("")
            #     ax.set_yticklabels([])
            
            fig.savefig(os.path.join(self.img_dir, f"Isync_{index}.png"), dpi=300, bbox_inches="tight")

            num += 1

    # def threetimeslice(self):

    #     fig1, axes = plt.subplots(1,3, figsize=(15,5))
    #     num=0

    #     for index in self.plot_indices:
    #         fig = self.plot_synchrotron_emission(index, norm=[-1,0,0], north_vec=[-0.703,0.452,-0.549])

    #         axes[num].imshow(fig.axes[0].images[0].get_array(), origin="lower", 
    #                         extent=fig.axes[0].images[0].get_extent(), cmap="plasma",
    #                         vmin=-1e5, vmax=1e15)
    #         num += 1

    #     fig1.savefig(os.path.join(self.img_dir, f"Isync_wr140.png"), dpi=300)


def main():
    data_path = '/mnt/massive-stars/data/thomas_simulations/wr140-sims/covertex_start/red_z_res/wr140-mhd-l7n256/'
    img_dir = "/home/visitor_ap4/code/project/images/sync-emission/"
    quantities = ["density", "pressure", "magnetic_field_magnitude", "Isync"]
    tolerance = 0.01

    yt_plot = YTPlotFunction(data_path, img_dir, quantities, tolerance)
    plot_indices = yt_plot.three_slice_indices()
    print(plot_indices)

    # yt_plot.plot_synchrotron_emission(plot_indices[1], norm=[-1,0,0], north_vec=[-0.703,0.452,-0.549])
    # yt_plot.threetimeslice()
    yt_plot.temp_plotter()

if __name__ == "__main__":
    main()