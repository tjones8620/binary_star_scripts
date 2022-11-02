import argparse
import numpy as np

class ArgparseInputs:
    def __init__(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("Path", help="Path to data folder")
        parser.add_argument("img_dir", help="Path to image directory")
        parser.add_argument("-q", "--fluidquantity", help="Desired fluid quantity to plot", default="Density")
        parser.add_argument("-s", "--surface", help ="3D Plot surface - XY, XZ, YZ", default="XY")
        parser.add_argument('--tol', nargs='+', type=float, default=[np.log10(1.0e-18), np.log10(2.0e-13)])
        parser.add_argument('--cmap', default='viridis')
        parser.add_argument("-mm", "--make_movie", type=bool, choices=[True, False], default=False)
        args = parser.parse_args()

        self.path = args.Path
        self.fluidquantity = args.fluidquantity
        self.surface = args.surface
        self.tolerance = args.tol
        self.cmap = args.cmap
        self.make_movie = args.make_movie
        self.img_dir = args.img_dir