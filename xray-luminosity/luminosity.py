__author__ = "Thomas Jones"

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.style.use("science")
import os
import astropy.units as u
import pathlib
home = pathlib.Path.home()
import palettable
from palettable.mycarta import LinearL_5
import argparse
from math import log10, floor


class ArgparseInputs:
    def __init__(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("data_path", help="Path to x-ray luminosity data file (.txt)")
        parser.add_argument("img_dir", help="Path to image directory")
        parser.add_argument('--tol', nargs='+', type=float, default=[0, 1e-3])
        parser.add_argument('--cmap', default=LinearL_5.mpl_colormap)
        parser.add_argument('--start_time', type=float, default=-1.25e7, help="Start time of simulation before periastron in seconds")
        parser.add_argument('--name', default='X-ray-Luminosity')
        parser.add_argument('--plot_xlims', nargs='+', type=float, default=[-100, 100], help="X-axis limits for plot")
        parser.add_argument('--plot_ylims', nargs='+', type=float, default=[0, None], help="Y-axis limits for plot")
        parser.add_argument('--data_xlims', nargs='+', type=float, default=[None, None], help="X-axis limits for data")
        parser.add_argument('--plot_hardsoftband', action='store_true', help="Plot the hard and soft X-ray bands")
        parser.add_argument('--plot_band', default="hard", choices=["hard", "soft"], help="Plot the hard or soft X-ray band")

        args = parser.parse_args()
        data_path = args.data_path
        self.data_path = data_path
        self.img_dir = args.img_dir
        self.tolerance = args.tol
        self.cmap = args.cmap
        self.start_time = args.start_time
        self.name = args.name
        self.plot_xlims = args.plot_xlims
        self.plot_ylims = args.plot_ylims
        self.data_xlims = args.data_xlims
        self.plot_hardsoftband = args.plot_hardsoftband
        self.plot_band = args.plot_band
    

############################################################################################################
# X-ray Luminosity
############################################################################################################
# This class is used to plot the total X-ray luminosity of a given source in different energy bands as a
# function of time. The luminosity and time values are take from a table in a .txt file, which contains 
# the luminosity values in different energy bands and the time values in seconds. This table is produced 
# with the total_luminosity.py script. The time values are converted to days and the luminosity values are
# in erg/s. The plot_xray_lum() method plots the total X-ray luminosity in different energy bands


class XrayLum:
    def __init__(self, path, image_folder, start_time, name):
        self.path = path
        self.image_folder = image_folder

        if not os.path.exists(self.image_folder):
            print(f"Creating folder {self.image_folder}...")
            input("Press any key to continue...")
            os.makedirs(self.image_folder)

        self.start_time = start_time
        self.name = name
        self.df = self.get_df(self.path, self.start_time)

    @staticmethod
    def get_df(path, start_time):
    
        headers = ['time','g0.1', 'g0.2', 'g0.3', 'g0.5', 'g1', 'g2', 'g5', 'g10']
        df = pd.read_csv(path, delim_whitespace=True, names=headers, index_col=False)
        df['2-10'] = df['g2'] - df['g10']
        df['0.3-2'] = df['g0.3'] - df['g2']
        df['time'] = (df['time'] + start_time) / (86400)
        df.rename(columns={'time': 'time [d]'}, inplace=True)
        return df

    def plot_xray_lum(self, **kwargs):

        fig, ax = plt.subplots()
        ax.set_xlabel('Time (d)')
        ax.grid(True)
        ax.plot(self.df['time [d]'], self.df['g0.1'], color='blue', linestyle='solid', label='$L(E>0.1keV)$')
        ax.plot(self.df['time [d]'], self.df['g0.2'], color='purple', linestyle='solid', label='$L(E>0.2keV)$')
        ax.plot(self.df['time [d]'], self.df['g0.5'], color='red', linestyle='solid', label='$L(E>0.5keV)$')
        ax.plot(self.df['time [d]'], self.df['g1'], color='yellow', linestyle='solid', label='$L(E>1keV)$')
        ax.plot(self.df['time [d]'], self.df['g2'], color='cyan', linestyle='solid', label='$L(E>2keV)$')
        ax.plot(self.df['time [d]'], self.df['g5'], color='lightpink', linestyle='solid', label='$L(E>5keV)$')
        ax.plot(self.df['time [d]'], self.df['g10'], color='peru', linestyle='solid', label='$L(E>10keV)$')
        ax.set_ylabel('X-ray Luminosity ($erg$ $s^{-1}$)')
        ax.legend(loc='upper left', ncol=2, framealpha=1, fontsize=5)
        ax.set_xlim(kwargs.get('xmin', None), kwargs.get('xmax', None))
        ax.set_ylim(kwargs.get('ymin', None), kwargs.get('ymax', None))
        plt.savefig(os.path.join(self.image_folder, f'xray-lum-{self.name}.png'), bbox_inches='tight', dpi=300)

    def plot_hardsoftband(self, ax, **kwargs):
        time = self.df['time [d]']
        lum1 = self.df['2-10']
        lum2 = self.df['0.3-2']
        dataxmin = kwargs.get('dataxmin', None)
        dataxmax = kwargs.get('dataxmax', None)
        indices = np.where(np.logical_and(dataxmin <= time, time <= dataxmax))[0]
        
        if dataxmin is not None and dataxmax is not None:
            time = time[indices]
            lum1 = lum1[indices]
            lum2 = lum2[indices]

        # order_of_mag = np.log10(np.max(lum1))
        order_of_mag = floor(np.log10(np.max(lum2)))

        ax.plot(time, lum1/10**order_of_mag, linestyle='solid', label='L(2keV - 10keV)')
        ax.plot(time, lum2/10**order_of_mag, linestyle='solid', label='L(0.3keV - 2keV)')
        # ax.legend(loc='upper right', ncol=1)
        ax.set_xlim(kwargs.get('xmin', None), kwargs.get('xmax', None))
        ax.set_ylim(kwargs.get('ymin', None), kwargs.get('ymax', None))
        ax.set_xlabel('Days relative to periastron', fontsize=14, labelpad=10)
        ax.set_ylabel('Unabsorbed Luminosity ' + '($10^{' + str(order_of_mag) + '}$' + ' $erg$ $s^{-1}$)', fontsize=14, labelpad=10)
        # plt.savefig(os.path.join(self.image_folder, f'xray-lum-{self.name}-hard-soft.png'), bbox_inches='tight', dpi=kwargs.get('dpi', 300))
        # return fig, ax

    def plot_band(self, ax, band = "hard", **kwargs):
        # fig, ax = plt.subplots(figsize=(5, 3.5))
        ax.set_xlabel('Days relative to periastron')
        ax.set_ylabel('Unabsorbed Luminosity ($10^{34}$ $erg$ $s^{-1}$)')

        if band == "hard":
            lum = self.df['2-10']
            bandlabel = 'L(2keV - 10keV)'
        elif band == "soft":
            lum = self.df['0.3-2']
            bandlabel = 'L(0.3keV - 2keV)'

        time = self.df['time [d]']
        dataxmin = kwargs.get('dataxmin', None)
        dataxmax = kwargs.get('dataxmax', None)
        indices = np.where(np.logical_and(dataxmin <= time, time <= dataxmax))[0]
        
        if dataxmin is not None and dataxmax is not None:
            time = time[indices]
            lum = lum[indices]

        ax.plot(time, lum/1e34, ls=kwargs.get('ls', 'solid'), label=kwargs.get("label", bandlabel), color=kwargs.get("color", None), lw=kwargs.get("lw", 1))
        ax.set_xlim(kwargs.get('xmin', None), kwargs.get('xmax', None))
        ax.set_ylim(kwargs.get('ymin', None), kwargs.get('ymax', None))
        # plt.savefig(os.path.join(self.image_folder, f'xray-lum-{self.name}-{band}.png'), bbox_inches='tight', dpi=kwargs.get('dpi', 300))
        # return fig, ax   
    
    def plot_obs_data(self, ax, **kwargs):
        pass

# def main():
    # path = os.path.join(home, "code/project/scripts/xray-luminosity/xray-tables/xray-lum-wr140-mhd-n256.txt")
    # start_time = -2.671e7
    # image_folder = os.path.join(home, "code/project/scripts/images/xray-lum")
    # plot = XrayLum(path=path, start_time=start_time, image_folder=image_folder, name='wr140-mhd-n256')
    # fig, ax = plot.plot_hardsoftband(xmin=-100,xmax=100, ymin=0, ymax = 1.5e35, dataxmin=-100, dataxmax=100)
    # # fig4, ax4 = plot.plot_band(xmin=-180,xmax=90, band="hard")

    # path = os.path.join(home, "code/project/scripts/xray-luminosity/xray-tables/xray-lum-wr140-hydro-n256.txt")
    # start_time = -7.25e6
    # image_folder = os.path.join(home, "code/project/scripts/images/xray-lum")
    # plot = XrayLum(path=path, start_time=start_time, image_folder=image_folder, name='wr140-hydro-n256')
    # fig, ax = plot.plot_hardsoftband(xmin=-100,xmax=100, ymin=0, ymax = 1.5e35, dataxmin=-65, dataxmax=100)

    # path = os.path.join(home, "code/project/scripts/xray-luminosity/xray-tables/xray-lum-wr140-hydro-accel-n256.txt")
    # start_time = -7.25e6
    # image_folder = os.path.join(home, "code/project/scripts/images/xray-lum")
    # plot = XrayLum(path=path, start_time=start_time, image_folder=image_folder, name='wr140-hydro-accel-n256')
    # fig, ax = plot.plot_hardsoftband(xmin=-100,xmax=100, ymin=0, ymax = 1.5e35, dataxmin=-50, dataxmax=50)
    # # fig, ax = plot.plot_band(xmin=-100,xmax=100, ymin=0, ymax = 1.5e35, dataxmin=-50, dataxmax=50, band="hard")

    # fig, axes = plt.subplots(figsize=(5, 3.5))
    # path = os.path.join(home, "code/project/scripts/xray-luminosity/xray-tables/xray-lum-wr140-mhd-compton-n256.txt")
    # start_time = -1.25e7
    # image_folder = os.path.join(home, "code/project/scripts/images/xray-lum")
    # name = 'wr140-mhd-compton-n256'
    # plot = XrayLum(path=path, start_time=start_time, image_folder=image_folder, name=name)
    # plot_xlims = [-100, 100]
    # plot_ylims = [0, 1.5e35]
    # data_xlims = [-100, 100]
    # plot.plot_hardsoftband(ax=axes, xmin=plot_xlims[0],xmax=plot_xlims[1], ymin=plot_ylims[0], ymax = plot_ylims[1], dataxmin=data_xlims[0], dataxmax=data_xlims[1])
    # plt.savefig(os.path.join(image_folder, f'xray-lum-{name}-hard-soft.png'), bbox_inches='tight', dpi=300)
    # fig, ax = plot.plot_band(xmin=-100,xmax=100, ymin=0, ymax = 1.5e35, dataxmin=-50, dataxmax=50, band="hard")


# def main():

#     cmdargs = ArgparseInputs()

#     fig, axes = plt.subplots(figsize=(5, 3.5))
#     path = cmdargs.data_path
#     start_time = cmdargs.start_time
#     img_dir = cmdargs.img_dir
#     name = cmdargs.name
#     plot = XrayLum(path=path, start_time=start_time, image_folder=img_dir, name=name)
#     plot_xlims = cmdargs.plot_xlims
#     plot_ylims = cmdargs.plot_ylims
#     data_xlims = cmdargs.data_xlims

#     if cmdargs.plot_hardsoftband:
#         plot.plot_hardsoftband(ax=axes, xmin=plot_xlims[0],xmax=plot_xlims[1], ymin=plot_ylims[0], ymax = plot_ylims[1], dataxmin=data_xlims[0], dataxmax=data_xlims[1])
#         plt.savefig(os.path.join(img_dir, f'xray-lum-{name}-hard-soft.png'), bbox_inches='tight', dpi=300)
#     if cmdargs.plot_band==True:
#         plot.plot_band(ax=axes, xmin=plot_xlims[0],xmax=plot_xlims[1], ymin=plot_ylims[0], ymax = plot_ylims[1], dataxmin=data_xlims[0], dataxmax=data_xlims[1], band=cmdargs.band)
#         plt.savefig(os.path.join(img_dir, f'xray-lum-{name}-{cmdargs.band}.png'), bbox_inches='tight', dpi=300)
    


def main():
    image_folder = os.path.join(home, "code/project/scripts/images/xray-lum")

    # axes.set_prop_cycle('color', LinearL_5.mpl_colors)

    mhd_1 = {"path":os.path.join(home, "code/project/scripts/xray-luminosity/xray-tables/xray-lum-wr140-mhd-n256.txt"), 
            "start_time":-2.671e7,
            "label": "mhd-1",
            "color": None,
            "dataxmin":-150,
            "dataxmax":120
            }
    hd_1 = {"path":os.path.join(home, "code/project/scripts/xray-luminosity/xray-tables/xray-lum-wr140-hydro-n256.txt"), 
            "start_time":-7.25e6,
            "label": "hd-1",
            "color": None,
            "dataxmin":-70,
            "dataxmax":120
            }
    hd_2 = {"path":os.path.join(home, "code/project/scripts/xray-luminosity/xray-tables/xray-lum-wr140-hydro-accel-n256.txt"), 
            "start_time":-7.25e6,
            "label": "hd-2",
            "color": None,
            "dataxmin":-70,
            "dataxmax":100
            }

    mhd_2 = {"path":os.path.join(home, "code/project/scripts/xray-luminosity/xray-tables/xray-lum-wr140-mhd-compton-n256.txt"), 
            "start_time":-1.25e7,
            "label": "mhd-2",
            "color": None,
            "dataxmin":-130,
            "dataxmax":135
            }

    sim_list = [mhd_1, hd_1, hd_2, mhd_2]

    fig, axes = plt.subplots(figsize=(15, 5), ncols=4, nrows=1, sharey=True, sharex=True) 
    i = 0 
    for sim in sim_list:
        plot = XrayLum(path=sim['path'], start_time=sim['start_time'], image_folder=image_folder, name=sim['label'])
        plot.plot_hardsoftband(ax=axes[i], xmin=-70,xmax=90, ymin=0, ymax = 1.4, dataxmin=sim['dataxmin'], dataxmax=sim['dataxmax'])
        
        axes[i].text(0.05, 0.95, sim['label'], transform=axes[i].transAxes, fontsize=14, verticalalignment='top', weight='bold')

        if not i == 0:
            axes[i].set_ylabel("")
        i += 1
    
    fig.subplots_adjust(wspace=0.05)

    fig.savefig(os.path.join(image_folder, f"xray-lum-soft-hard.png"), bbox_inches='tight', dpi=300)


if __name__ == "__main__":
    main()
