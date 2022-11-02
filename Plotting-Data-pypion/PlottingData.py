from pypion_plotter import Plot_Functions
from argparse_pypionplotter import ArgparseInputs

inputs = ArgparseInputs()

def main():
    plot = Plot_Functions(inputs.path, inputs.img_dir, inputs.fluidquantity, inputs.tolerance, inputs.surface)
    plot.evolution = plot.evolution[::50]

    if plot.N_dims==3:
        plot.ThreeDSurfacePlotter(movie=inputs.make_movie, colormap=inputs.cmap)

if __name__=="__main__":
    main()

