import os
import glob
import numpy as np
import re
from pypion.ReadData import ReadData

def make_snapshots(data_path):
        ########## Cataloging silo files ###############################
        cwd = os.getcwd()
        # Get the list of silo files
        silo_files = os.listdir(data_path)
        silo_files = sorted([f for f in silo_files if f.endswith('.silo')])
        filename=(silo_files[0]).replace('_level00_0000.00000000.silo','')


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
        os.chdir(cwd)

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
                
        # Bundle silo files of different levels of same time instant into a snapshot.
        evolution = np.array(catalog).T
        print(f"Shape of evolution array: {evolution.shape}")
        return evolution