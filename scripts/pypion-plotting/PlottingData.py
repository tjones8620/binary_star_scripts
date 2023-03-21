from pypion_plotter import Plot_Functions
from argparse_pypionplotter import ArgparseInputs

inputs = ArgparseInputs()

def main():
    # plot = Plot_Functions(inputs.path, inputs.img_dir, inputs.fluidquantity, inputs.tolerance, inputs.surface, period=7.992, start_time=-2.671e7)
    # plot = Plot_Functions(inputs.path, inputs.img_dir, inputs.fluidquantity, inputs.tolerance, inputs.surface, period=7.992, start_time=-7.25e6)

    # Sets up the plotter object. Note that the period and start_time are optional arguments. 
    # If not specified, the period is set to 1 and the start_time is set to 0. So, if you want to
    # have the correct orbital phase, you need to specify the period and start_time.
    
    # The tolerance (poorly name) is the vmin and vmax for the colorbar in the 3D surface plot.
    plot = Plot_Functions(path=inputs.path, 
                        basename=inputs.basename, 
                        image_dir=inputs.img_dir, 
                        fluid_quantity=inputs.fluidquantity, 
                        tolerance=inputs.tolerance, 
                        plane=inputs.surface, 
                        period=7.992, 
                        start_time=-1.239e7)
    
    # plot.evolution = plot.evolution[50:]

    if plot.N_dims==3:
        plot.ThreeDSurfacePlotter(movie=inputs.make_movie, 
                                  colormap=inputs.cmap, 
                                  plot_inset=True, 
                                  inset_zoom=4)

    # plot.three_time_slice(log=inputs.log, colormap=inputs.cmap, d_phase=0.007, zoom=1, plot_inset=True)

if __name__=="__main__":
    main()

