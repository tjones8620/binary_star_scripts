import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.style.use("science")
import os
import astropy.units as u
import pathlib
home = pathlib.Path.home()
from scipy.optimize import curve_fit

import sys
sys.path.insert(0, os.path.join(home, 'code/project/scripts/'))


class XrayLum:
    def __init__(self, path, image_folder, start_time, name):
        self.path = path
        self.image_folder = image_folder
        self.start_time = start_time
        self.name = name
        self.df = self.get_df(self.path, self.start_time)

    @staticmethod
    def get_df(path, start_time):
        headers = ['time','g0.1', 'g0.2', 'g0.3', 'g0.5', 'g1', 'g2', 'g5', 'g10']
        df = pd.read_csv(path, delim_whitespace=True, names=headers, index_col=False)
        df['2-10'] = df['g2'] - df['g10']
        df['time'] = (df['time'] + start_time) / (86400)
        df.rename(columns={'time': 'time [d]'}, inplace=True)
        return df

    def plot_xray_lum(self, **kwargs):
        fig, ax1 = plt.subplots()
        ax1.set_xlabel('Time (d)')
        ax1.set_ylabel('Unabsorbed Flux ($erg$ $cm^{-2}$$s^{-1}$)')
        ax1.grid(True)

        ax2 = ax1.twinx()  # create a second axes that shares the same x-axis
        ax2.plot(self.df['time [d]'], self.df['g0.1'], color='blue', linestyle='solid', label='$L(E>0.1keV)$')
        ax2.plot(self.df['time [d]'], self.df['g0.2'], color='purple', linestyle='solid', label='$L(E>0.2keV)$')
        ax2.plot(self.df['time [d]'], self.df['g0.5'], color='red', linestyle='solid', label='$L(E>0.5keV)$')
        ax2.plot(self.df['time [d]'], self.df['g1'], color='yellow', linestyle='solid', label='$L(E>1keV)$')
        ax2.plot(self.df['time [d]'], self.df['g2'], color='cyan', linestyle='solid', label='$L(E>2keV)$')
        ax2.plot(self.df['time [d]'], self.df['g5'], color='lightpink', linestyle='solid', label='$L(E>5keV)$')
        ax2.plot(self.df['time [d]'], self.df['g10'], color='peru', linestyle='solid', label='$L(E>10keV)$')
        ax2.set_ylabel('X-ray Luminosity ($erg$ $s^{-1}$)')
        ax2.legend(loc='upper left', ncol=2, framealpha=1, frameon=True, fontsize=5)
        ax2.set_xlim(kwargs.get('xmin', None), kwargs.get('xmax', None))
        ax2.set_ylim(kwargs.get('ymin', None), kwargs.get('ymax', None))
        plt.savefig(os.path.join(self.image_folder, f'xray-lum-{self.name}.png'), bbox_inches='tight', dpi=300)

    def plot_xray_lum_2_10(self, **kwargs):
        fig, ax1 = plt.subplots()
        ax1.set_xlabel('Time (d)')
        ax1.set_ylabel('Unabsorbed Flux ($erg$ $cm^{-2}$$s^{-1}$)')
        ax1.grid(True)

        ax2 = ax1.twinx()  # create a second axes that shares the same x-axis
        ax2.plot(self.df['time [d]'], self.df['2-10'], color='blue', linestyle='solid', label='$L(2keV - 10keV)$')
        ax2.set_ylabel('X-ray Luminosity ($erg$ $s^{-1}$)')
        ax2.legend(loc='upper left', ncol=1, framealpha=1)
        ax2.set_xlim(kwargs.get('xmin', None), kwargs.get('xmax', None))
        ax2.set_ylim(kwargs.get('ymin', None), kwargs.get('ymax', None))
        plt.savefig(os.path.join(self.image_folder, f'xray-lum-{self.name}-2-10keV.png'), bbox_inches='tight', dpi=300)

def main():
    # path = os.path.join(home, "code/project/scripts/xray-luminosity/xray-lum-wr140-mhd-n256.txt")
    # start_time = -2.671e7
    # image_folder = os.path.join(home, "code/project/images/xray-luminosity/wr140-mhd-n256")

    path = os.path.join(home, "code/project/scripts/xray-luminosity/xray-lum-wr140-hydro-n256.txt")
    start_time = -2.671e7
    image_folder = os.path.join(home, "code/project/images/xray-luminosity/wr140-hydro-n256")

    plot = XrayLum(path=path, start_time=start_time, image_folder=image_folder, name='wr140-hydro-n256')
    plot.plot_xray_lum()
    plot.plot_xray_lum_2_10()

if __name__ == "__main__":
    main()
