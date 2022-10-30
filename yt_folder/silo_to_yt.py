from pypion.ReadData import ReadData
import yt
import os
import glob
import numpy as np
import re

def make_snapshots(data_path):
        ########## Cataloging silo files ###############################
        # Get the list of silo files
        silo_files = os.listdir(data_path)
        silo_files = sorted([f for f in silo_files if f.endswith('.silo')])
        filename=(silo_files[0]).replace('_level00_0000.00000000.silo','')
        print(filename)

        cwd = os.getcwd()

        os.chdir(data_path)
        file_list = glob.glob('*.silo', recursive=True)

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
                files = sorted(glob.glob(f"{filename}_level{level_list[i]}_0000.*.silo"))
                catalog.append(files)
                
        
        # catalog = [os.path.join(data_path, file) for file in catalog[i] for i in range(len(catalog))]

        # Bundle silo files of different levels of same time instant into a snapshot.
        evolution = np.array(catalog).T
        print(f"Shape of evolution array: {evolution.shape}")

        os.chdir(cwd)
        return evolution

def get_ds(file) -> yt.data_objects.static_output.Dataset:
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

    grid_data = [dict(left_edge=[0.5-0.5**(n+1)]*len(N_grids), right_edge=[0.5+0.5**(n+1)]*len(N_grids), level=n, dimensions=N_grids) for n in range(N_levels)]

    i = 0
    for g in grid_data:
        g[("gas", "density")] = (data_den[i], "g/cm**3")
        g[("gas", "temperature")] = (data_temp[i], "K")
        g[("gas", "velocity_x")] = (data_velx[i], "cm/s")
        g[("gas", "velocity_y")] = (data_vely[i], "cm/s")
        g[("gas", "velocity_z")] = (data_velz[i], "cm/s")
        # g[("gas", "NG_Mask")] = (data_ngmask[i], "cm/s")
        i += 1

    ds = yt.load_amr_grids(grid_data, N_grids, length_unit=f"{Dom_size} * cm", geometry=("cartesian", ("z","y","x")), sim_time=sim_time)
    return ds