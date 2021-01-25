#!/usr/bin/env python
from __future__ import division, print_function
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser
import csv




def getTimeSteps(args):
    timeSteps = []
    with open(str(args.obs) + "-" + str(args.noradid) + ".csv") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=",")
        for row in csv_reader:
            timeSteps.append(int(row[0]))

    return timeSteps

def main(args):

    timeSteps = getTimeSteps(args)

    stacked_img = []
    #stacked_img.append(np.zeros((200,200)))
    counter =1 
    for t in timeSteps:
        hduH = fits.open(str(args.obs)+ "-2m-" + str(t) + "h-" + str(args.channel).zfill(4) + "-dirty.fits" )
        hduT = fits.open(str(args.obs)+ "-2m-" + str(t) + "t-" + str(args.channel).zfill(4) + "-dirty.fits" )
        diff = hduH[0].data[0,0,:,:] - hduT[0].data[0,0,:,:]
        stacked_img.append(diff)

        plt.subplot(121)
        plt.title("img {}/{}".format(counter, len(timeSteps)))
        plt.imshow(diff, origin="lower")
        plt.colorbar()

        plt.subplot(122)
        plt.imshow(np.mean(np.array(stacked_img), axis=0), origin="lower")
        plt.colorbar()
        plt.savefig("img" + str(counter).zfill(4) + ".png")
        plt.clf()
        #plt.show()
        
        counter += 1
        



if __name__ == "__main__":
    parser = ArgumentParser("Test", description="plot individual image for the stacked freq")
    parser.add_argument("--obs", type=int, help="the observatin id")
    parser.add_argument("--noradid",type=int, help="the norad id of the satellite")
    parser.add_argument("--channel", type=int, help="the freq channel")
    args = parser.parse_args()
    main(args)





