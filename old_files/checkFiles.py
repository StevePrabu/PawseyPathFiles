#!/usr/bin/env python
from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from tqdm import tqdm
from astropy.io import fits
import os

def main(args):

    print("checking if files exists")

    timeSteps = np.arange(args.t1, args.t2 + 1)
    waterfall = np.zeros((args.channels, int(len(timeSteps))))
    for f in tqdm(range(args.channels)):
        counter = 0
        for t in timeSteps:
            #print("working on timestep {} and channel {}".format(t,f))
            if os.path.isfile(str(args.obs) + "-" + str(args.midName)+ "-" + str(t) + "-" + str(f).zfill(4) + "-dirty.fits"):
                try:
                    hdu = fits.open(str(args.obs) + "-" + str(args.midName)+ "-" + str(t) + "-" + str(f).zfill(4) + "-dirty.fits")
                except:
                    print("error")
                    print(str(args.obs) + "-" + str(args.midName)+ "-" + str(t) + "-" + str(f).zfill(4) + "-dirty.fits")
                    
                data = hdu[0].data[0,0,:,:]
                if np.all(data == 0):
                    waterfall[f,counter] = -1
                else:
                    waterfall[f,counter] = 1
                counter += 1
                continue
            else:
                print("file not found {}".format(str(args.obs) + "-" + str(args.midName)+ "-" + str(t) + "-" + str(f).zfill(4) + "-dirty.fits"))
                counter +=1
    plt.imshow(waterfall, origin="lower",aspect="auto")
    plt.show()
    np.save(waterfall, "filecheck.npy")
    

if __name__ == "__main__":
    parser = ArgumentParser("fileCheck", description="checks if all files exists")
    parser.add_argument("--obs",required=True,type=int,help="the obs id")
    parser.add_argument("--t1",required=True,type=int, help="the start timeStep")
    parser.add_argument("--t2",required=True,type=int,help="the end timeStep")
    parser.add_argument("--channels",required=True,type=int,help="the number of channels")
    parser.add_argument("--midName",default="2m",help="the midname of the file")
    args = parser.parse_args()

    main(args)

