__author__ = "Thomas Jones"

from pypion.ReadData import ReadData
import yt
import os
import glob
import numpy as np
import re

def make_snapshots(data_path):
    """
    Function to make snapshots of silo files.

    Input: path to silo files
    Output: list of snapshots - each snapshot is a list of silo files of different levels at the same time instant.

    Example:
    evolution = make_snapshots(data_path)
    
    """
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

def get_ds(file, quantities=["density"], **kwargs) -> yt.data_objects.static_output.Dataset:

    """
    Function to create a yt dataset from a silo file.
    Input: snapshot of silo file (list of filepaths for a single snapshot), list of quantities to be loaded
    Output: yt dataset

    Example:
    evolution = make_snapshots(data_path)
    file = evolution[0]
    ds = get_ds(file, quantities=["density", "temperature", "velocity", "magnetic_field"])

    """

    data = ReadData(file)
    sim_time = data.sim_time()
    N_levels = data.nlevels()
    N_grids = data.ngrid()[::-1]
    Dom_size = max(data.level_max()) - min(data.level_min())
    print(f"Simulation Info: {N_levels} levels, {N_grids} grids, {Dom_size} cm")

    # Arrays for quantities
    if quantities.__contains__("density"):
        data_den = data.get_3Darray("Density")['data']
    if quantities.__contains__("temperature"):
        data_temp = data.get_3Darray("Temperature")['data']
    if quantities.__contains__("pressure"):
        data_pressure = data.get_3Darray("Pressure")['data']
    if quantities.__contains__("NG_Mask"):
        data_ngmask = data.get_3Darray("NG_Mask")['data']
    if quantities.__contains__("windtracer"):
        data_windtr = data.get_3Darray("Tr000_WIND")['data']
    if quantities.__contains__("velocity"):
        data_velx = data.get_3Darray("VelocityX")['data']
        data_vely = data.get_3Darray("VelocityY")['data']
        data_velz = data.get_3Darray("VelocityZ")['data']
    if quantities.__contains__("magnetic_field"):
        data_Bx = data.get_3Darray("MagneticFieldX")['data']
        data_By = data.get_3Darray("MagneticFieldY")['data']
        data_Bz = data.get_3Darray("MagneticFieldZ")['data']

    grid_size = np.array(N_grids)/max(N_grids) # Normalized grid size

    grid_data = np.array([dict(left_edge=np.array([0.5-0.5**(n+1)]*len(N_grids))*grid_size, 
                            right_edge=np.array([0.5+0.5**(n+1)]*len(N_grids))*grid_size, 
                            level=n, dimensions=N_grids) for n in range(N_levels)])

    bbox = np.array([np.zeros(3), grid_size]).T

    i = 0
    for g in grid_data:
        if quantities.__contains__("density"):
            g[("gas", "density")] = (data_den[i], "g/cm**3")
        if quantities.__contains__("temperature"):
            g[("gas", "temperature")] = (data_temp[i], "K")
        if quantities.__contains__("pressure"):
            g[("gas", "pressure")] = (data_pressure[i], "dyne/cm**2")
        if quantities.__contains__("NG_Mask"):
            g[("gas", "NG_Mask")] = (data_ngmask[i], "K")
        if quantities.__contains__("windtracer"):
            g[("gas", "windtracer")] = (data_windtr[i], "")
        if quantities.__contains__("velocity"):
            g[("gas", "velocity_x")] = (data_velx[i], "cm/s")
            g[("gas", "velocity_y")] = (data_vely[i], "cm/s")
            g[("gas", "velocity_z")] = (data_velz[i], "cm/s")
        if quantities.__contains__("magnetic_field"):
            g[("gas", "magnetic_field_x")] = (data_Bx[i], "G")
            g[("gas", "magnetic_field_y")] = (data_By[i], "G")
            g[("gas", "magnetic_field_z")] = (data_Bz[i], "G")

        i += 1

    ds = yt.load_amr_grids(grid_data, N_grids, length_unit=f"{Dom_size} * cm", geometry=("cartesian", ("z","y","x")), sim_time=sim_time+kwargs.get("start_time", 0), bbox=bbox, time_unit="s", unit_system="cgs")
    return ds



# def get_ts(files, quantities=["density"]):
    
#     ds_list = []
#     for file in files:
#         ds = get_ds(file, quantities)
#         ds_list.append(ds)
#         del ds

#     global ts    
#     ts = yt.DatasetSeries(ds_list)

#     return ts

def get_ts(evolution, **kwargs):

    ds_list = []

    start = kwargs.get('start', 1)
    end = kwargs.get('end', 10)
    step = kwargs.get('step', 1)

    # Load the desired snapshots
    for i in range(start, end, step):
        ds_list.append(get_ds(evolution[i], quantities=kwargs.get("quantities", ["density"]), start_time=kwargs.get("start_time", 0)))

    print("Number of datasets: ", len(ds_list))
    # Create time series object
    return yt.DatasetSeries(ds_list)
    del ds_list


