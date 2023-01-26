from silo_to_yt import *
import sys
import pickle
yt.set_log_level("ERROR")
import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np
import pandas as pd
import sigfig
from yt.visualization.volume_rendering.api import PointSource
# print(plt.style.available)
# plt.style.use('nature')
# plt.rcParams['font.family'] = 'sans-serif'
# plt.rcParams['font.sans-serif'] = 'Helvetica'
# import matplotlib as mpl
# mpl.rcParams['font.family'] = 'serif'

from matplotlib import ticker
from matplotlib import cm
import sys
import pathlib
home = str(pathlib.Path.home())
sys.path.insert(0, os.path.join(home, "code/project/scripts/"))
from binary_system_calculations.binary_velocity import BinarySystem
from misc.make_movies import make_movies

class YTPlotFunction():

    def __init__(self, data_path, img_dir, quantities=["density"], period=7.992, start_time=0):

        self.evolution = make_snapshots(data_path) # list of all snapshots
        self.data_path = data_path # path to data
        self.img_dir = img_dir # path to save images
        if not os.path.exists(self.img_dir):
            os.makedirs(self.img_dir)
            print("Created images directory:{} ".format(self.img_dir))
        else:
            print("Images directory already exists")
        self.quantities = quantities # list of quantities to include in dataset
        self.period = period*u.yr # period of binary system
        self.start_time = start_time*u.s # time of first snapshot wrt periastron

    def add_synchrotron_emission(self, ds):
        from unyt import K, g, cm, s
        def _Isync(field, data):
            Isync = data["magnetic_field_magnitude"]**(3/2) * data["pressure"] * data["NG_Mask"]
            return Isync * (K*g**(7/4)/(cm**(3/4)*s**(7/2)))**-1 * cm

        # add synchrotron emission field to dataset
        ds.add_field(("gas", "Isync"), function=_Isync, units="auto", sampling_type="cell", force_override=True) 

    def three_slice_indices(self, d_phase=0.01): 
        
        """
        Function that finds the index of SILO snapshots closest to a phase of
        1-d_phase, 1, and 1+d_phase. This is used to plot the three slices
        of the binary system at the pre-periastron, periastron, and post-periastron

        Parameters
        ----------

        d_phase : float

            phase difference between the three slices
        
        Returns
        -------
        plot_indices : list
        """

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
        return plot_indices

    def plot_synchrotron_emission(self, i, north_vec, norm, **kwargs):

        """
        Function that creates synchrotron emission maps using the
        yt.OffAxisProjectionPlot function. Converts the yt figure
        into a mpl figure and returns the mpl fig object.

        Parameters
        ----------
        i : int

            index of snapshot to plot

        north_vec : list

            list of three floats that define the north vector of the
            projection plot
        
        norm : list

            list of three floats that define the normal vector to the
            projection plot

        **kwargs : dict
        
            dictionary of keyword arguments to pass to the yt.OffAxisProjectionPlot
        
        Returns
        -------
        fig : matplotlib.figure.Figure

        """

        print(f"Loading dataset: Sim Snapshot {i}")
        ds = get_ds(self.evolution[i], quantities=["magnetic_field", "pressure", "NG_Mask"], start_time=self.start_time)
        time = ds.current_time.value
        phase = ((self.period + (time*u.s).to(u.yr))/self.period).value
        print(f"Successfully Loaded dataset: {str(ds)} at orbital phase {phase:.2f}")
        self.add_synchrotron_emission(ds)
        print(f"Added Synchrotron emission field")
        print(f"Plotting Isync for snapshot {i}...")
        prj = yt.OffAxisProjectionPlot(ds, normal=norm, north_vector=north_vec, fields=("gas", "Isync")) # create projection plot
        prj.set_cmap(("gas", "Isync"), "gist_heat")
        prj.hide_colorbar("Isync")
        prj.set_figure_size(5)
        prj.zoom(kwargs.get("zoom", 1))
        prj.set_zlim(("gas", "Isync"), 1e13, 5e15)
        # prj.set_font({"family": "sans-serif", "size": 8})
        

        # phase = ((self.period + (time*u.s).to(u.yr))/self.period).value

        fig = prj.export_to_mpl_figure((1,1), cbar_mode=kwargs.get("cbarmode", "None")) # export to matplotlib figure
        ax = fig.axes[0] 
        ax.set_xlabel("x (AU)") 
        ax.set_ylabel("y (AU)")
        st = r"$\phi$ = " + f"{phase:.2f}" 
        ax.text(0.8, 0.9, st, color="black", fontsize=12, transform=ax.transAxes, 
                    bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=1', alpha=0.5))
        return fig 

    def temp_plotter(self, norm, north_vec, zoom=1, d_phase=0.01):
        self.plot_indices = self.three_slice_indices(d_phase = d_phase) # list of indices for pre-periastron, periastron, and post-periastron snapshots
        num = 0
        for index in self.plot_indices:
            fig =  self.plot_synchrotron_emission(index, norm, north_vec, zoom=zoom)
            ax = fig.axes[0]
            if num == 1 or num == 2:
                ax.set_ylabel("")
                ax.set_yticklabels([])

            if num == 2:
                cbaxes = fig.add_axes([ax.get_position().x1-0.04,ax.get_position().y0-0.034,0.06,ax.get_position().height-0.037])
                cb = fig.colorbar(cm.ScalarMappable(cmap="gist_heat"), ax=ax, orientation="vertical", cax=cbaxes, pad=0.0)
                cb.set_label(r"$I_{sync}$ (Normalised Units)", fontsize=18, labelpad=8)
                cb.ax.tick_params(labelsize=12)
                tick_locator = ticker.MaxNLocator(nbins=5)
                cb.locator = tick_locator
                cb.update_ticks()
            print(f"Saving image Isync_{index}.png...")
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
        # slc.set_font({"family": "times new roman"})
        slc.set_font({"family": "mpl-default", "size": 11, "weight": "bold"})
        # slc.set_xlabel("x (AU)")
        fig = slc.export_to_mpl_figure((1,1), cbar_mode="None") # export to matplotlib figure
        ax = fig.axes[0]
        ax.set_xlabel("$\mathrm{x}$ (AU)")
        ax.set_ylabel("y (AU)")
        wr140 = BinarySystem(0.8993, 2895, 10.31, 29.27, "wr140")
        orb_sol = wr140.orb_sol

        x1, y1, x2, y2, vx1, vy1, vx2, vy2 = orb_sol
        x1_au, y1_au = (x1*u.cm).to(u.au).value, (y1*u.cm).to(u.au).value
        x2_au, y2_au = (x2*u.cm).to(u.au).value, (y2*u.cm).to(u.au).value
        ax.plot(x1_au, y1_au, label="WR star", lw=2)
        ax.plot(x2_au, y2_au, label="O star", lw=2)
        ax.legend(loc="upper right", frameon=True, framealpha=1, facecolor="white", fontsize=14)
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

        print("Getting star positions from trajectory file...")
        # print(trajectory_file, times)

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
        print(times)

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

    def volume_rendering(self, **kwargs):
        from yt.units import cm
        # Creating volume renderings at each time step in ts
        i=0
        for file in self.evolution[40:]:
            ds = get_ds(file)
            sc = yt.create_scene(ds)
            # Print the time of the current scene
            print(f"Time: {ds.current_time.to('s')}")
            time = np.float64(ds.current_time.to('s').value)
            print(type(time))
            # if kwargs['trajectory_file'] is not None:
                # star1_x, star1_y, star2_x, star2_y, star1_z, star2_z = self.get_star_position(kwargs['trajectory_file'], [time])
            # identifying the source
            source = sc[0]
            source.set_field('density')
            source.set_log(True)


            # building transfer function
            bounds = (1e-18, 1e-11)
            tf = yt.ColorTransferFunction(x_bounds=np.log10(bounds), nbins=500)
            # Automatically add a number of layers
            tf.add_layers(8, w=0.002, colormap='viridis')
            
            source.tfh.tf = tf
            source.tfh.bounds = bounds
            
            cam = sc.camera
            cam
            cam.zoom(2.2)
            cam.resolution = (4096, 4096)
            cam.switch_orientation(normal_vector=[-1,0,0], north_vector=[0.2,1,0])

            # colors = np.random.random([1, 4])
            # colors[:, 3] = 1.0
        

            # points = PointSource(np.array([star1_x*cm, star1_y*cm, star1_z*cm]), colors=colors, radii=5e12)
            # sc.add_source(points)

            sc.save(os.path.join(self.img_dir, f"wr140_mhd_compton_n256_{i+1}"), sigma_clip=6.0)
            print(f"Saving image {os.path.join(self.img_dir, f'wr140_mhd_compton_n256_{i+1}.png')}")
            del sc, cam
            i+=1

        make_movies(self.img_dir, self.img_dir, "wr140_mhd_compton_n256_density.mp4")

# def main():
#     data_path = '/mnt/massive-stars/data/thomas_simulations/wr140-sims/covertex_start/red_z_res/wr140-mhd-l7n256/'
#     # data_path = "/mnt/massive-stars/data/thomas_simulations/wr140-sims/compton_cool/compton_mhd_l7n128_v2"
#     img_dir = "/home/visitor_ap4/code/project/images/sync-emission/"
#     # img_dir = "/home/visitor_ap4/code/project/images/"
#     quantities = ["density", "pressure", "magnetic_field_magnitude"]
#     # quantities = ["windtracer"]
#     # start = -7.25e6
#     start=-2.671e7

#     yt_plot = YTPlotFunction(data_path, img_dir, quantities, start_time=start)


#     # yt_plot.plot_synchrotron_emission(plot_indices[1], norm=[-1,0,0], north_vec=[-0.703,0.452,-0.549])
#     # yt_plot.threetimeslice()
#     yt_plot.temp_plotter(norm=[-1,0,0], north_vec=[-0.703,0.452,-0.549], zoom=4)

#     # yt_plot.plot_smr_plot(50)

#     # # yt_plot.temp_quantity_plotter("temperature", ["temperature"], norm=[-1,0,0], north_vec=[0,0,1])
#     # yt_plot.threetimeslice("temperature", ["temperature"], norm=[-1,0,0], north_vec=[0,0,1])

def main():
    data_path = "/mnt/massive-stars/data/thomas_simulations/wr140-sims/compton_tests/compton_mhd_n256"
    # data_path = '/mnt/massive-stars/data/thomas_simulations/wr140-sims/covertex_start/red_z_res/wr140-mhd-l7n256/'
    # img_dir = "/home/visitor_ap4/code/project/scripts/images/sync-emission/wr140_mhd_n256"
    img_dir = "/home/visitor_ap4/code/project/scripts/images/sync-emission/compton_mhd_n256"
    quantities = ["density", "pressure", "magnetic_field_magnitude"]
    start=-1.24e7
    # start=-2.671e7

    yt_plot = YTPlotFunction(data_path, img_dir, quantities, start_time=start)
    yt_plot.temp_plotter(norm=[-1,0,0], north_vec=[0.703,0.452,-0.549], zoom=6, d_phase = 0.01)


# def main():

#     data_path = "/mnt/massive-stars/data/thomas_simulations/wr140-sims/compton_cool/compton_mhd_l7n128_large_eta/"
#     img_dir = "/home/visitor_ap4/code/project/scripts/images/"
#     quantities = ["windtracer"]
#     start=-2.671e7
#     yt_plot = YTPlotFunction(data_path, img_dir, quantities, start_time=start)
#     yt_plot.plot_smr_plot(0)

# def main():

#     data_path = "/mnt/massive-stars/data/thomas_simulations/wr140-sims/compton_cool/compton_mhd_l7n128_v2"
#     img_dir = "/home/visitor_ap4/code/project/images/"
#     quantities = ["windtracer", "magnetic_field"]
#     start = -7.25e6
#     yt_plot = YTPlotFunction(data_path, img_dir, quantities, start_time=start)
#     yt_plot.plot_smr_plot(50)

# def main():
#     # data_path = "/mnt/massive-stars/data/thomas_simulations/wr140-sims/covertex_start/orig_res/wr140-hydro-cool-n064/"
#     data_path = "/mnt/massive-stars/data/thomas_simulations/wr140-sims/covertex_start/orig_res/wr140-hydro-cool-n128/"
#     img_dir = "/home/visitor_ap4/code/project/images/volume_renderings/wr140-hydro-cool-n128/density/"

#     yt_plot = YTPlotFunction(data_path, img_dir)
#     yt_plot.volume_rendering(start=50, end=100, step=10)


if __name__ == "__main__":
    main()