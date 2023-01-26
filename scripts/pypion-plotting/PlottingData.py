from pypion_plotter import Plot_Functions
from argparse_pypionplotter import ArgparseInputs

inputs = ArgparseInputs()

def main():
    # plot = Plot_Functions(inputs.path, inputs.img_dir, inputs.fluidquantity, inputs.tolerance, inputs.surface, period=7.992, start_time=-2.671e7)
    # plot = Plot_Functions(inputs.path, inputs.img_dir, inputs.fluidquantity, inputs.tolerance, inputs.surface, period=7.992, start_time=-7.25e6)



    plot = Plot_Functions(inputs.path, inputs.img_dir, inputs.fluidquantity, inputs.tolerance, inputs.surface, period=7.992, start_time=-1.239e7)
    plot.evolution = plot.evolution[0:1]

    if plot.N_dims==3:
        plot.ThreeDSurfacePlotter(movie=inputs.make_movie, colormap=inputs.cmap)

    # plot.three_time_slice(log=inputs.log, colormap=inputs.cmap, d_phase=0.01, zoom=1, plot_inset=False)

if __name__=="__main__":
    main()

