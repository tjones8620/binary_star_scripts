import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.style.use("science")
# plt.rcParams['mathtext.default'] = 'rm'
import os
import astropy.units as u
import pathlib
home = pathlib.Path.home()


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
        ax.legend(loc='upper left', ncol=2, framealpha=1, frameon=True, fontsize=5)
        ax.set_xlim(kwargs.get('xmin', None), kwargs.get('xmax', None))
        ax.set_ylim(kwargs.get('ymin', None), kwargs.get('ymax', None))
        plt.savefig(os.path.join(self.image_folder, f'xray-lum-{self.name}.png'), bbox_inches='tight', dpi=300)

    def plot_hardsoftband(self, log=False, **kwargs):
        fig, ax = plt.subplots(figsize=(5, 3.5))
        ax.set_xlabel('Days relative to periastron')
        ax.set_ylabel('Unabsorbed Luminosity ($erg$ $s^{-1}$)')

        if log==True:
            ax.set_yscale('log')
            ax.set_ylim(1e27, 1e37)
        else:
            ax.set_yscale('linear')

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

        ax.plot(time, lum1, linestyle='solid', label='L(2keV - 10keV)')
        ax.plot(time, lum2, linestyle='solid', label='L(0.3keV - 2keV)')
        ax.legend(loc='upper left', ncol=1)
        ax.set_xlim(kwargs.get('xmin', None), kwargs.get('xmax', None))
        ax.set_ylim(kwargs.get('ymin', None), kwargs.get('ymax', None))
        plt.savefig(os.path.join(self.image_folder, f'xray-lum-{self.name}-hard-soft.png'), bbox_inches='tight', dpi=kwargs.get('dpi', 300))
        return fig, ax

    def plot_band(self, band = "hard", **kwargs):
        fig, ax = plt.subplots(figsize=(5, 3.5))
        ax.set_xlabel('Days relative to periastron')
        ax.set_ylabel('Unabsorbed Luminosity ($erg$ $s^{-1}$)')

        if band == "hard":
            lum = self.df['2-10']
            label = 'L(2keV - 10keV)'
        elif band == "soft":
            lum = self.df['0.3-2']
            label = 'L(0.3keV - 2keV)'

        time = self.df['time [d]']
        dataxmin = kwargs.get('dataxmin', None)
        dataxmax = kwargs.get('dataxmax', None)
        indices = np.where(np.logical_and(dataxmin <= time, time <= dataxmax))[0]
        
        if dataxmin is not None and dataxmax is not None:
            time = time[indices]
            lum = lum[indices]

        ax.plot(time, lum, linestyle='solid', label=label)
        ax.set_xlim(kwargs.get('xmin', None), kwargs.get('xmax', None))
        ax.set_ylim(kwargs.get('ymin', None), kwargs.get('ymax', None))
        plt.savefig(os.path.join(self.image_folder, f'xray-lum-{self.name}-{band}.png'), bbox_inches='tight', dpi=kwargs.get('dpi', 300))
        return fig, ax   
    
    def plot_hard_vs_obs():
        pass

def main():
    path = os.path.join(home, "code/project/scripts/xray-luminosity/xray-tables/xray-lum-wr140-mhd-n256.txt")
    start_time = -2.671e7
    image_folder = os.path.join(home, "code/project/images/xray-luminosity/wr140-mhd-n256")
    plot = XrayLum(path=path, start_time=start_time, image_folder=image_folder, name='wr140-mhd-n256')
    fig, ax = plot.plot_hardsoftband(xmin=-100,xmax=100, ymin=0, ymax = 1.5e35, dataxmin=-100, dataxmax=100)
    # fig4, ax4 = plot.plot_band(xmin=-180,xmax=90, band="hard")

    path = os.path.join(home, "code/project/scripts/xray-luminosity/xray-tables/xray-lum-wr140-hydro-n256.txt")
    start_time = -7.25e6
    image_folder = os.path.join(home, "code/project/images/xray-luminosity/wr140-hydro-n256")
    plot = XrayLum(path=path, start_time=start_time, image_folder=image_folder, name='wr140-hydro-n256')
    fig, ax = plot.plot_hardsoftband(xmin=-100,xmax=100, ymin=0, ymax = 1.5e35, dataxmin=-65, dataxmax=100)

    path = os.path.join(home, "code/project/scripts/xray-luminosity/xray-tables/xray-lum-wr140-hydro-accel-n256.txt")
    start_time = -7.25e6
    image_folder = os.path.join(home, "code/project/images/xray-luminosity/wr140-hydro-accel-n256")
    plot = XrayLum(path=path, start_time=start_time, image_folder=image_folder, name='wr140-hydro-accel-n256')
    fig, ax = plot.plot_hardsoftband(xmin=-100,xmax=100, ymin=0, ymax = 1.5e35, dataxmin=-50, dataxmax=50)
    # fig, ax = plot.plot_band(xmin=-100,xmax=100, ymin=0, ymax = 1.5e35, dataxmin=-50, dataxmax=50, band="hard")

    path = os.path.join(home, "code/project/scripts/xray-luminosity/xray-tables/xray-lum-wr140-compton-n128-new-updated.txt")
    start_time = -1.25e7
    image_folder = os.path.join(home, "code/project/images/xray-luminosity/wr140-compton-n128-new")
    plot = XrayLum(path=path, start_time=start_time, image_folder=image_folder, name='wr140-mhd-compton-n128')
    fig, ax = plot.plot_hardsoftband(xmin=-100,xmax=100, ymin=0, ymax = 1.5e35, dataxmin=-100, dataxmax=100)
    # fig, ax = plot.plot_band(xmin=-100,xmax=100, ymin=0, ymax = 1.5e35, dataxmin=-50, dataxmax=50, band="hard")


if __name__ == "__main__":
    main()
