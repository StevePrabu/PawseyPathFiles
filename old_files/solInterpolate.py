from mwapy import aocal
import numpy as np
from argparse import ArgumentParser

def nan_helper(y):
	return np.isnan(y), lambda z: z.nonzero()[0]

def main(args):

    ao = aocal.fromfile(args.inputFile)
    
    for tile in range(128):
        tileSolution = ao[0,tile,:,:]

        ## check if tile is flagged
        if np.all(np.isnan(tileSolution)):
            print("tile " + str(tile) + " is flagged")
            continue

        else:
            nans, x = nan_helper(tileSolution)
            tileSolution[nans] = np.interp(x(nans), x(~nans), tileSolution[~nans])
            ao[0,tile,:,:] = tileSolution

    ao.tofile(str(args.obs) + ".bin")








if __name__ == "__main__":
    parser = ArgumentParser("SolInterpolater", description="interpolates the calibration solution for the flagged frequencies")
    parser.add_argument("--inputFile", required=True, help="the name of the flagged solution file")
    parser.add_argument("--obs", required=True, type=int, help="the observation id")
    args = parser.parse_args()
    main(args)
