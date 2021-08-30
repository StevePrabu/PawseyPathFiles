#!/usr/bin/env python
from tqdm import tqdm
from astropy.io import fits
import numpy as np
from argparse import ArgumentParser
import csv


def main(args):
    
    ## make timeStep array by reading csv file
    timeSteps = []
    with open(str(args.obs) + "-" + str(args.noradid) + ".csv") as csv_file:
         csv_reader = csv.reader(csv_file, delimiter=",")
         for row in csv_reader:
             timeSteps.append(int(row[0]))
    
    cube = []
    for f in tqdm(range(args.channels)):
        global_data = []
        for t in timeSteps:
            hduH = fits.open(str(args.obs)+ str(args.band) + "-" + str(args.midName) + "-" + str(t) + "h-" + str(f).zfill(4) + "-dirty.fits" )
            hduT = fits.open(str(args.obs)+ str(args.band) + "-" + str(args.midName) + "-" + str(t) + "t-" + str(f).zfill(4) + "-dirty.fits" )
            diff = hduH[0].data[0,0,:,:] - hduT[0].data[0,0,:,:] 
            global_data.append(diff)
        stack = np.mean(np.array(global_data), axis=0)
        cube.append(stack)
    np.save(str(args.noradid) + "-" + str(args.obs) + ".npy", cube)

        
             


if __name__ == "__main__":
    parser = ArgumentParser("make cube", description="stackes images and makes a cube for later analysis")
    parser.add_argument("--obs", required=True,  help="the observation id")
    parser.add_argument("--band", default="", help="the band name")
    parser.add_argument("--noradid",required=True, type=int, help="the norad id")
    parser.add_argument("--channels",default=768,type=int, help="the number of channels")
    parser.add_argument("--midName", default="2m", help="the mid name of the fits files")
    args = parser.parse_args()

    main(args)
    
    
