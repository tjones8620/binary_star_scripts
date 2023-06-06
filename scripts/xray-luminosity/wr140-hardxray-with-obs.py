from luminosity import XrayLum
import matplotlib.pyplot as plt
import os
from palettable.mycarta import LinearL_5
import pathlib
plt.style.use('default')
plt.rcParams['font.family'] = 'sans-serif'
home = str(pathlib.Path.home())
import pandas as pd

def main():
    image_folder = os.path.join(home, "code/project/scripts/images/xray-lum")
    fig, axes = plt.subplots(figsize=(2.5, 3.5))
#     axes.set_prop_cycle('color', LinearL_5.mpl_colors)

    mhd_1 = {"path":os.path.join(home, "code/project/scripts/xray-luminosity/xray-tables/xray-lum-wr140-mhd-n256.txt"), 
            "start_time":-2.671e7,
            "label": "mhd-1",
            "color": "magenta",
            "dataxmin":-150,
            "dataxmax":135,
            "lw": 1,
            "ls": "-"
            }
    hd_1 = {"path":os.path.join(home, "code/project/scripts/xray-luminosity/xray-tables/xray-lum-wr140-hydro-n256.txt"), 
            "start_time":-7.25e6,
            "label": "hd-1",
            "color": "C0",
            "dataxmin":-75,
            "dataxmax":120,
            "lw": 1,
            "ls": "-."
            }
    hd_2 = {"path":os.path.join(home, "code/project/scripts/xray-luminosity/xray-tables/xray-lum-wr140-hydro-accel-n256.txt"), 
            "start_time":-7.24e6,
            "label": "hd-2",
            "color": "C6",
            "dataxmin":-50,
            "dataxmax":100,
            "lw": 1,
                "ls": "-"
            }

    mhd_2 = {"path":os.path.join(home, "code/project/scripts/xray-luminosity/xray-tables/xray-lum-wr140-mhd-compton-n256.txt"), 
            "start_time":-1.25e7,
            "label": "mhd-2",
            "color": "C1",
            "dataxmin":-135,
            "dataxmax":140,
            "lw": 1,
            "ls": "-"
            }

    sim_list = [mhd_1, hd_1, hd_2, mhd_2]

    for sim in sim_list:
        plot = XrayLum(path=sim['path'], start_time=sim['start_time'], image_folder=image_folder, name=sim['label'])
        plot.plot_band(xmin=-130,xmax=130, ymax = 5,dataxmin=sim['dataxmin'], dataxmax=sim['dataxmax'], band="hard", ax=axes, label=sim['label'], color=sim['color'], lw=sim['lw'], ls=sim['ls'])
        # axes.legend()
    axes.legend(loc="upper right", fontsize=8)

    axes.set_ylim(0, None)

    # Make sure scattered points are in front of the line
    # axes.scatter([-100, -95, -90], [1, 1.02, 1.04], color='k', label='HD-1', zorder=10, s=3)
    # df_path = os.path.join(home, "code/project/scripts/xray-luminosity/wr140-phasefolded-unabslum.csv")
    # df = pd.read_csv(df_path)
    # t_rxte = df['t_rxte'].values
    # lum_rxte = df['l_rxte'].values
    # t_nicer = df['t_nicer'].values
    # lum_nicer = df['l_nicer'].values
    # t_swiftpc = df['t_swiftpc'].values
    # lum_swiftpc = df['l_swiftpc'].values
    # t_swiftwt = df['t_swiftwt'].values
    # lum_swiftwt = df['l_swiftwt'].values

    
    # axes.plot(t_rxte, lum_rxte, 'ko', label='RXTE', markersize=3)
    # axes.plot(t_nicer, lum_nicer, 'go', label='NICER', markersize=3)
    # axes.plot(t_swiftpc, lum_swiftpc, 'bo',label='Swift PC',markersize=3)
    # axes.plot(t_swiftwt, lum_swiftwt, 'bo',label='Swift WT', markerfacecolor='none', markersize=3)
    # axes.legend(loc="upper left")



    fig.savefig(os.path.join(image_folder, f'xray-lum-combined-hard-narrow.png'), bbox_inches='tight', dpi=300)


def main():
    image_folder = os.path.join(home, "code/project/scripts/images/xray-lum")
    fig, axes = plt.subplots(figsize=(5, 3.5))
#     axes.set_prop_cycle('color', LinearL_5.mpl_colors)


    # mhd_2_n128 = {"path":os.path.join(home, "code/project/scripts/xray-luminosity/xray-tables/xray-lum-wr140-compton-n128-new-updated.txt"), 
    #         "start_time":-1.25e7,
    #         "label": "n128",
    #         "color": None,
    #         "dataxmin":-130,
    #         "dataxmax":120,
    #         "lw": 1.5,
    #         "ls": "-"
    #         }
    mhd_1 = {"path":os.path.join(home, "code/project/scripts/xray-luminosity/xray-tables/xray-lum-wr140-mhd-n256.txt"), 
            "start_time":-2.671e7,
            "label": "mhd-1 (No IC)",
            "color": "magenta",
            "dataxmin":-150,
            "dataxmax":150,
            "lw": 1.5,
            "ls": "-"
            }

    mhd_2 = {"path":os.path.join(home, "code/project/scripts/xray-luminosity/xray-tables/xray-lum-wr140-mhd-compton-n256.txt"), 
            "start_time":-1.25e7,
            "label": "mhd-2 (With IC)",
            "color": None,
            "dataxmin":-130,
            "dataxmax":160,
            "lw": 1.5,
            "ls": "-"
            }

    sim_list = [mhd_1, mhd_2]

    for sim in sim_list:
        plot = XrayLum(path=sim['path'], start_time=sim['start_time'], image_folder=image_folder, name=sim['label'])
        plot.plot_band(xmin=-150, xmax=160, dataxmin=sim['dataxmin'], dataxmax=sim['dataxmax'], band="hard", ax=axes, label=sim['label'], color=sim['color'], lw=sim['lw'], ls=sim['ls'])
        # axes.legend()

    axes.legend()
    axes.set_ylim(0, None)

    df_path = os.path.join(home, "code/project/scripts/xray-luminosity/wr140-phasefolded-unabslum.csv")
    df = pd.read_csv(df_path)
    t_rxte = df['t_rxte'].values
    lum_rxte = df['l_rxte'].values
    t_nicer = df['t_nicer'].values
    lum_nicer = df['l_nicer'].values
    t_swiftpc = df['t_swiftpc'].values
    lum_swiftpc = df['l_swiftpc'].values
    t_swiftwt = df['t_swiftwt'].values
    lum_swiftwt = df['l_swiftwt'].values

    
    axes.plot(t_rxte, lum_rxte, 'ko', label='RXTE ($2-10$ keV)', markersize=1.5)
    # axes.plot(t_nicer, lum_nicer, 'go', label='NICER', markersize=3)
    # axes.plot(t_swiftpc, lum_swiftpc, 'bo',label='Swift PC',markersize=1)
    # axes.plot(t_swiftwt, lum_swiftwt, 'bo',label='Swift WT', markerfacecolor='none', markersize=3)
    axes.legend(loc="upper left", fontsize=8)
    
    # Make facecolor white with background transparent
    fig.patch.set_facecolor('white')
    fig.patch.set_alpha(0.0)

    # fig.savefig(os.path.join(image_folder, f'xray-mhd2-withobslum.png'), bbox_inches='tight', dpi=300)
    fig.savefig(os.path.join(image_folder, f'xray-mhd2-withobslum-poster.png'), bbox_inches='tight', dpi=900)



def main():
    image_folder = os.path.join(home, "code/project/scripts/images/xray-lum")
    fig, axes = plt.subplots(figsize=(5, 3.5))
#     axes.set_prop_cycle('color', LinearL_5.mpl_colors)


    mhd_1 = {"path":os.path.join(home, "code/project/scripts/scripts/xray-luminosity/xray-tables/xray-lum-wr140-mhd-n256.txt"), 
            "start_time":-2.671e7,
            "label": "mhd-1",
            "color": "magenta",
            "dataxmin":-150,
            "dataxmax":150,
            "lw": 1.3,
            "ls": "-.",
            "zorder": 3
            }

    mhd_2 = {"path":os.path.join(home, "code/project/scripts/scripts/xray-luminosity/xray-tables/xray-lum-wr140-mide-compton.txt"), 
            "start_time":-1.25e7,
            "label": "mhd-2",
            "color": "lime",
            "dataxmin":-130,
            "dataxmax":160,
            "lw": 1.3,
            "ls": "-",
            "zorder": 2
            }

    mhd_3 = {"path":os.path.join(home, "code/project/scripts/scripts/xray-luminosity/xray-tables/xray-lum-wr140-mide-control.txt"), 
            "start_time":-1.25e7,
            "label": "mhd-3",
            "color": "blue",
            "dataxmin":-130,
            "dataxmax":160,
            "lw": 1.3,
            "ls": "--",
            "zorder": 1
            }
        
    

    sim_list = [mhd_1, mhd_2, mhd_3]

    for sim in sim_list:
        plot = XrayLum(path=sim['path'], start_time=sim['start_time'], image_folder=image_folder, name=sim['label'])
        plot.plot_band(xmin=-125, xmax=115, dataxmin=sim['dataxmin'], dataxmax=sim['dataxmax'], band="hard", ax=axes, label=sim['label'], color=sim['color'], lw=sim['lw'], ls=sim['ls'])
        # Change plot zorder
        axes.lines[-1].set_zorder(sim['zorder'])
        # axes.legend()

    axes.legend()
    axes.set_ylim(0, None)

    df_path = os.path.join(home, "code/project/scripts/scripts/xray-luminosity/wr140-phasefolded-unabslum.csv")
    df = pd.read_csv(df_path)
    t_rxte = df['t_rxte'].values
    lum_rxte = df['l_rxte'].values
 
    axes.plot(t_rxte, lum_rxte, 'ko', label='RXTE', markersize=1)
    axes.legend(loc="upper left", fontsize=8)
    axes.grid()
    
    # Turn off box in figure legend
    # axes.get_legend().get_frame().set_linewidth(0.0)
    fig.savefig(os.path.join(image_folder, f'xray-lum-withobs-papernew.png'), bbox_inches='tight', dpi=300)


def main():
    image_folder = os.path.join(home, "code/project/scripts/images/xray-lum")
    fig, axes = plt.subplots(figsize=(5, 3.5))
#     axes.set_prop_cycle('color', LinearL_5.mpl_colors)


    mhd_2 = {"path":os.path.join(home, "code/project/scripts/scripts/xray-luminosity/xray-tables/xray-lum-wr140-mide-compton.txt"), 
            "start_time":-1.25e7,
            "label": r"$256\times256\times64$",
            "color": "lime",
            "dataxmin":-130,
            "dataxmax":160,
            "lw": 1.3,
            "ls": "-",
            "zorder": 2
            }

    mhd_3 = {"path":os.path.join(home, "code/project/scripts/scripts/xray-luminosity/xray-tables/xray-lum-wr140-mide-compton-128.txt"), 
            "start_time":-1.25e7,
            "label": r"$128\times128\times32$",
            "color": "blue",
            "dataxmin":-130,
            "dataxmax":160,
            "lw": 1.3,
            "ls": "--",
            "zorder": 1
            }
        
    

    sim_list = [mhd_2, mhd_3]

    for sim in sim_list:
        plot = XrayLum(path=sim['path'], start_time=sim['start_time'], image_folder=image_folder, name=sim['label'])
        plot.plot_band(xmin=-125, xmax=115, ymin=0, ymax=3, dataxmin=sim['dataxmin'], dataxmax=sim['dataxmax'], band="hard", ax=axes, label=sim['label'], color=sim['color'], lw=sim['lw'], ls=sim['ls'])
        # Change plot zorder
        axes.lines[-1].set_zorder(sim['zorder'])
        # axes.legend()

    axes.legend()
    axes.set_ylim(0, None)

    axes.legend(loc="upper left", fontsize=8)
    axes.grid()
    
    # Turn off box in figure legend
    # axes.get_legend().get_frame().set_linewidth(0.0)
    fig.savefig(os.path.join(image_folder, f'xray-lum-mhd2-resolution-study2.png'), bbox_inches='tight', dpi=300)


if __name__ == "__main__":
    main()