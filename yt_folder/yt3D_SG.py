# Author: Sam Green, Created: 07-02-22

#-------------------------------
# Libraries needed:
#-------------------------------
from pypion.ReadData import ReadData
from pypion.argparse_command import InputValues
#-------------------------------
#-------------------------------
import yt
import numpy as np
import math
#-------------------------------
# Astro stuff:
#-------------------------------
from astropy import units as u
#-------------------------------
#-------------------------------
line = InputValues()
time_dicts = line.time_dicts
dimen = line.dimen
#-------------------------------
#-------------------------------

#var1 = ["Density", -22, -27, "viridis", 'y', 127]

for files in time_dicts:

    arr = time_dicts[files]

    data_den = ReadData(arr).get_3Darray("Density")['data']
    data_temp = ReadData(arr).get_3Darray("Temperature")['data']
    #lim_max = (ReadData(arr).get_3Darray("Density")['max_extents'] * u.cm).to(u.pc)
    #lim_min = (ReadData(arr).get_3Darray("Density")['min_extents'] * u.cm).to(u.pc)
    #sim_age = ReadData(arr).get_3Darray("Density")['sim_time']

    '''
    for k in range(3):
            denk = data_den[k]
            tempk = data_temp[k]
            print(k)

            for j in range(256):
                denj = denk[j]
                tempj = tempk[j]
                    
                for i in range(256):
                    deni = denj[i]
                    tempi = tempj[i]

                    for m in range(256):
                        data_den[k][j][i][m] = (deni[m] / 2.3e-24) ** 2 * 2.0e-23 * math.exp(-1.0e7 / tempi[m])
                    #print(data_den[k][j][i])

    #print(data_den)
    '''

    grid_data = [
        dict(
            left_edge=[0.0, 0.0, 0.0],
            right_edge=[1.0, 1.0, 1.0],
            level=0,
            dimensions=[128, 128, 128],
            ),
        dict(
            left_edge=[0.5, 0.5, 0.5],
            right_edge=[1.0, 0.5, 0.5],
            level=1,
            dimensions=[128, 128, 128],
            ),
        dict(
            left_edge=[0.25, 0.25, 0.25],
            right_edge=[1.0, 0.25, 0.25],
            level=2,
            dimensions=[128 ,128, 128],
            ),
    ]

    i = 0
    for g in grid_data:
        g[("gas", "density")] = (data_den[i], "g/cm**3")
        i += 1

    #print(grid_data)

    ds = yt.load_amr_grids(grid_data, [128, 128, 128], length_unit="7e14 * cm", refine_by=1024)


    #data = dict(density = (data_den[0], "g/cm**3"))
    #bbox = np.array([[lim_min[0][0].value, lim_max[0][0].value], [lim_min[0][1].value, lim_max[0][1].value], [lim_min[0][2].value, lim_max[0][2].value]])

    #ds = yt.load_uniform_grid(data, data_den[0].shape, length_unit="pc", bbox=bbox, nprocs=256, sim_time=sim_age)

    #slc = yt.SlicePlot(ds, "y", ("gas", "density"))
    #slc.set_cmap(("gas", "density"), "viridis")
    #slc.annotate_grids(cmap=None)
    #slc.save("/mnt/data/yt_test.png")

    sc = yt.create_scene(ds, field=("gas", "density"))

    #source = sc[0]
    #source.tfh.set_bounds((1.5e-23, 6e-23))
    #source.tfh.set_log(True)
    #source.tfh.grey_opacity = True
    #source.tfh.plot("/mnt/data/transfer_function.png", profile_field=("gas", "density"))

    
    source = sc[0]

    source.set_field(("gas", "density"))
    source.set_log(True)
    #source.grey_opacity = True

    #bounds = (1e-23, 3e-23)
    #bounds = (8e-29, 2e-26)
    bounds = (1.6e-23, 6e-23)
    
    tf = yt.ColorTransferFunction(np.log10(bounds))
    tf.sample_colormap(np.log10(1.6e-23), w=0.01, colormap="viridis")
    #tf.add_gaussian(np.log10(3e-23), width=0.005, height=[0.753, 1.0, 0.933, 1.0])

    #tf.add_layers(8, colormap="viridis")

    def linramp(vals, minval, maxval):
            return (vals - vals.min()) / (vals.max() - vals.min())


    tf.map_to_colormap(
        np.log10(1.6e-23), np.log10(6e-23), colormap="viridis", scale_func=linramp
    )

    source.tfh.tf = tf
    source.tfh.bounds = bounds
    

    #sc.annotate_domain(ds, color=[1, 1, 1, 0.001])
    sc.annotate_axes(alpha=0.5)

    cam = sc.camera

    # Rotates the camera to look at (x, z) = 1, y = 0.
    #cam.position = np.array([0.5, 1.0, 0.5])
    #cam.north_vector = [1.0, 0.0, 0.0]
    #cam.rotate(np.pi / 4.0)
    #cam.yaw(-np.pi / 4.0)

    frame = 0
    source.tfh.plot("/mnt/local/thomas/yt_test_tf.png", profile_field=("gas", "density"))
    # sc.save("/mnt/local/thomas/yt_test_%04i.png" % frame, sigma_clip=6.0)
    # frame += 1


    # Rotate by 180 degrees over 5 frames
    for _ in cam.iter_rotate(np.pi, 5):
       sc.save("/mnt/local/thomas/yt_3Dtest_%04i.png" % frame, sigma_clip=6.0)
       frame += 1

