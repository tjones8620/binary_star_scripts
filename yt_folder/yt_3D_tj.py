
from silo_to_yt import make_snapshots, get_ds
import yt
import numpy as np
import os

yt.set_log_level("ERROR")

cwd = os.path.abspath(os.path.dirname(__file__))

# Get the list of silo files
data_path = "/mnt/massive-stars/data/thomas_simulations/wr140-sims/covertex_start/orig_res/wr140-hydro-cool-n064/"


# Get the list of silo files
evolution = make_snapshots(data_path)


# Load the desired snapshots
ds_list = []
for i in range(100,125):
    ds_list.append(get_ds(evolution[i]))

# Get list of fields
fields = ds_list[0].field_list
print("\nList of fields: ", fields)

# Create time series object
ts = yt.DatasetSeries(ds_list)
del ds_list

#In[]
# Creating the image directory
img_dir = os.path.join(cwd, "yt_images/VolumeRenders")
if not os.path.exists(img_dir):
    os.makedirs(img_dir)
else:
    print("Images directory already exists")

# Creating volume renderings at each time step in ts
for i in range(len(ts)):
    sc = yt.create_scene(ts[i])
    
    # identifying the source
    source = sc[0]
    source.set_field('density')
    source.set_log(True)

    # building transfer function
    bounds = (1e-18, 1e-12)
    tf = yt.ColorTransferFunction(x_bounds=np.log10(bounds))
    tf.add_layers(5, w=0.005, colormap='viridis')
    
    source.tfh.tf = tf
    source.tfh.bounds = bounds
    
    cam = sc.camera
    # cam.zoom(4)
    # cam.switch_orientation(normal_vector=[1, 0, 0], north_vector=[0, 1, 0])

    sc.save(os.path.join(img_dir, f"wr140_hyd_n064_{i+1}"), sigma_clip=6.0)
    
    del sc, tf, cam, source