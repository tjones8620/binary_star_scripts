import time
import numpy as np
import os
from make_movies import make_movies
from pypion.ReadData import ReadData
import pickle
import json
import argparse
import numpy as np
import yt
import glob
import re


class ArgparseInputs:
    def __init__(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("data_path", help="Path to data folder")
        parser.add_argument("pickle_file", help = "Full path to pickle file")
        args = parser.parse_args()

        self.data_path = args.data_path
        self.pickle_file = args.pickle_file

cwd = os.getcwd()
    
def make_snapshots(data_path):
        ########## Cataloging silo files ###############################
        # Get the list of silo files
        silo_files = os.listdir(data_path)
        silo_files = sorted([f for f in silo_files if f.endswith('.silo')])
        filename=(silo_files[0]).replace('_level00_0000.00000000.silo','')
        
        print("Info from silo files:")
        print(f"Basename of silo files: {filename}")

        file_list = glob.glob(os.path.join(data_path, '*.silo'), recursive=False)

        level_list = []
        files = []

        for file in file_list:
            level = re.search('_level(.*)_', file)
            if level == None:
                pass
            else:
                level = level.group(1)
                if not level in level_list:
                    level_list.append(level)
        level_list.sort()

        ########## Categorizing data files into levels #################
        if len(level_list) == 1: 
            print('Simulation Info: Single level')
            catalog = []
            files = sorted(glob.glob(filename + '_0000.*.silo'))
            catalog.append(files)
        else:
            print(f'Simulation Info: {len(level_list)} levels')
            catalog = []
            for i in range(len(level_list)):
                files = sorted(glob.glob(os.path.join(data_path, f"{filename}_level{level_list[i]}_0000.*.silo")))
                catalog.append(files)
                

        # Bundle silo files of different levels of same time instant into a snapshot.
        evolution = np.array(catalog).T
        print(f"Number of snapshots: {evolution.shape[0]}")

        return evolution



def save_grid_data_to_pickle(files, pickle_file):

    grid_data_list = []
    sim_time_list = []
    N_grids_list = []
    Dom_size_list = []

    #### Creating yt compatable grid data for each snapshot ###
    for file in files:
        data = ReadData(file)

        sim_time = data.sim_time()
        N_levels = data.nlevels()
        N_grids = data.ngrid()
        Dom_size = max(data.level_max()) - max(data.level_min())

        data_den = data.get_3Darray("Density")['data']
        data_temp = data.get_3Darray("Temperature")['data']
        data_velx = data.get_3Darray("VelocityX")['data']
        data_vely = data.get_3Darray("VelocityY")['data']
        data_velz = data.get_3Darray("VelocityZ")['data']
        # data_ngmask = data.get_3Darray("NG_Mask")['data']

        grid_data = [dict(left_edge=[0.5-0.5**(n+1)]*len(N_grids), 
                        right_edge=[0.5+0.5**(n+1)]*len(N_grids), 
                        level=n, dimensions=N_grids) for n in range(N_levels)]

        i = 0
        for g in grid_data:
            g[("gas", "density")] = (data_den[i], "g/cm**3")
            g[("gas", "temperature")] = (data_temp[i], "K")
            g[("gas", "velocity_x")] = (data_velx[i], "cm/s")
            g[("gas", "velocity_y")] = (data_vely[i], "cm/s")
            g[("gas", "velocity_z")] = (data_velz[i], "cm/s")
            # g[("gas", "NG_Mask")] = (data_ngmask[i], "cm/s")
            i += 1

        grid_data_list.append(grid_data)
        sim_time_list.append(sim_time)
        N_grids_list.append(N_grids)
        Dom_size_list.append(Dom_size)

        print(f"\nCreating Grid Data at sim-time: {sim_time} ...")

    print("\n")
    print("#"*30)
    print("Finished creating Grid Data")
    print("\n")
    print("#"*30)

    # Creating dictionary containing grid data and sim-time for each snapshot, and N_grids and Domain Size for the simulation
    grid_dict = {'grid_data': grid_data_list, 'sim_time': sim_time_list, 'N_grids': N_grids_list[0], 'Dom_Size': Dom_size_list[0]}


    print(f"\nSaving array grid information to {pickle_file} ...")
    
    os.chdir(cwd)

    if not pickle_file.endswith('.pickle'):
       print("Pickle file name should end with '.pickle'")
       print("Exiting ...")
       exit()

    #### Saving dictionary to a pickle file ####
    with open(pickle_file, 'xb') as f:
        pickle.dump(grid_dict, f)
    
    print(f"Finished saving data to {pickle_file}")


def yt_ts_from_pickle(pickle_file):
    with open(pickle_file, 'rb') as f:
        grid_dict = pickle.load(f)
    
    grid_data_list = grid_dict['grid_data']
    sim_time_list = grid_dict['sim_time']
    N_grids = grid_dict['N_grids']
    Dom_size = grid_dict['Dom_Size']

    ds_list = []
    for grid_data, sim_time in zip(grid_data_list, sim_time_list):
        ds = yt.load_amr_grids(grid_data, N_grids, length_unit=f"{Dom_size} * cm", geometry=("cartesian", ("z","y","x")), sim_time=sim_time)
        ds_list.append(ds)
        del ds
    
    ts = yt.DatasetSeries(ds_list)
    return ts
    


def main():

    cmdargs = ArgparseInputs()

    data_path = cmdargs.data_path
    pickle_file = cmdargs.pickle_file

    # Get the list of silo files
    evolution = make_snapshots(data_path)

    save_grid_data_to_pickle(evolution, pickle_file)

    # ts = yt_ts_from_pickle(pickle_file)

if __name__ == "__main__":
    main()
