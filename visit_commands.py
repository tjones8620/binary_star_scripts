import os


file_directory = "/mnt/massive-stars/data/thomas_simulations/wr140-sims/covertex_start/red_z_res/wr140-mhd-l7n128"

level0 = os.path.join(file_directory, "wr140_mhd_cool_d3l7n128_level00*.silo")


OpenDatabase(f"mimir:{level0}", 0)
