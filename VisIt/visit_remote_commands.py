import argparse
import os

# class ArgumentInputs:
#     def __init__(self):
#         parser = argparse.ArgumentParser()
#         parser.add_argument('path')
#         parser.add_argument('basename')
#         parser.add_argument('quantity')
#         args = parser.parse_args()
        
#         self.path = args.path
#         self.basename = args.basename
#         self.quantity = args.quantity

# cmdargs = ArgumentInputs()

hostname = "mimir.ap.dias.ie"
# file_path = "/mnt/massive-stars/data/thomas_simulations/wr140-sims/periast_start/new-mass-loss/wr140-mhd-n128/"
# file_basename = "wr140_mhd_cool_d3l7n128"

# Input parameters
# file_path = "/mnt/massive-stars/data/thomas_simulations/wr140-sims/covertex_start/red_z_res/wr140-hydro-l7n128/"
file_path = "/mnt/massive-stars/data/thomas_simulations/wr140-sims/late_start/wr140-hydro-n256/"
file_basename = "wr140_hydro_cool_d3l7n256"
quantity = "Temperature"

# file_path = cmdargs.path
# file_basename = cmdargs.basename
# quantity = cmdargs.quantity

# List to save full names of all databases
database_list = []

# Opening database for each level
for i in range(7):
    database = f"{hostname}:{os.path.join(file_path, file_basename)}_level0{i}_0000.*.silo database"
    database_list.append(database)
    OpenDatabase(database, 0)
    AddPlot("Pseudocolor", quantity, 1, 1)

# Creating a correlation between databases
CreateDatabaseCorrelation("Correlation01", database_list, 0)


# Creating 3 slice plot
AddOperator("ThreeSlice", 1)
ThreeSliceAtts = ThreeSliceAttributes()
ThreeSliceAtts.x = 0
ThreeSliceAtts.y = 0
ThreeSliceAtts.z = 0
ThreeSliceAtts.interactive = 1
SetOperatorOptions(ThreeSliceAtts, -1, 1)


