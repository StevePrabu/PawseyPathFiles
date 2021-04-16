#!/usr/bin/env python
from __future__ import division, print_function
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from astropy.wcs import WCS
from tqdm import tqdm
from subprocess import call
import os

def floodfill(xs, ys, floodfillValue, data, noise, imgSize):
    """
    performs the forest fire algorithm to source find
    """
    
    q = []
    q.append([xs, ys])

    while q:

        x, y = q.pop()
        searchMap[x,y] = 1
        snrMap[x,y] = data[x,y]/noise

        if checkValid(x+1, y, floodfillValue, data, noise, imgSize):
            q.append([x+1, y])
        
        if checkValid(x, y+1, floodfillValue, data, noise, imgSize):
            q.append([x, y+1])
        
        if checkValid(x-1, y, floodfillValue, data, noise, imgSize):
            q.append([x-1, y])

        if checkValid(x, y-1, floodfillValue, data, noise, imgSize):
            q.append([x, y-1])

        if checkValid(x+1, y+1, floodfillValue, data, noise, imgSize):
            q.append([x+1, y+1])

        if checkValid(x+1, y-1, floodfillValue, data, noise, imgSize):
            q.append([x+1, y-1])

        if checkValid(x-1, y-1, floodfillValue, data, noise, imgSize):
            q.append([x-1, y-1])

        if checkValid(x-1, y+1, floodfillValue, data, noise, imgSize):
            q.append([x-1, y+1])
    
    return



def checkValid(x,y,floodfillValue, data, noise, imgSize):
    """
    checks validity of seed pixel
    """

    if 1 < x < (imgSize-1) and 1 < y < (imgSize-1) and \
    data[x,y]/noise >= floodfillValue and searchMap[x,y] == 0:
        return True
    else:
        return False

def main(args):
    
    ## open file and extract all parameters from header
    hdu = fits.open(args.fileName)
    data = hdu[0].data[0,0,:,:]
    wcs = WCS(hdu[0].header, naxis=2)
    imgSize = hdu[0].header["NAXIS1"]
    unit = hdu[0].header['BUNIT']

    ## estimate noise
    tmp = np.copy(data)
    tmp[np.abs(tmp) >3*np.std(tmp)] = 0
    tmp[np.abs(tmp) >3*np.std(tmp)] = 0
    tmp[np.abs(tmp) >3*np.std(tmp)] = 0
    noise = np.std(tmp)

    if args.verbose:
        print("estimated global noise in image {} {}".format(noise, unit))

    ## create global SNR map
    global snrMap, searchMap
    snrMap = np.zeros((imgSize, imgSize))
    searchMap = np.zeros((imgSize, imgSize))
    
    ## start source finding
    seeds = np.asarray(np.where(data >= args.seedSigma*noise)).T
    
    counter = 0
    for seed in tqdm(seeds):
        counter += 1
        floodfill(seed[0], seed[1], args.floodfillSigma, data, noise, imgSize)

    if args.verbose:
        print("{} seed events found".format(counter))

    ## save snrmap to fits file
    if os.path.isfile("snr-{}".format(args.fileName)):
        bashSyn1 = "rm snr-{}".format(args.fileName)
        bashExecute1 = call(bashSyn1,shell=True)
    hdu_ouput = fits.PrimaryHDU(snrMap,header=hdu[0].header)
    hdu_ouput.writeto("snr-{}".format(args.fileName))
    
    ## plot
    plt.subplot(211, projection=wcs)
    plt.title("SNR of detected events")
    plt.imshow(snrMap, origin="lower")
    plt.grid(color='white', linestyle="dotted")
    plt.colorbar().set_label("SNR")
    plt.xlabel("RA (Degrees)")
    plt.ylabel("DEC (Degrees)")

    plt.subplot(212, projection=wcs)
    plt.title("input image")
    plt.imshow(data, origin="lower")
    plt.grid(color='white', linestyle="dotted")
    plt.colorbar().set_label("Jy/beam")
    plt.xlabel("RA (Degrees)")
    plt.ylabel("DEC (Degrees)")


    plt.show()









if __name__ == "__main__":
    parser = ArgumentParser("SourceFinder", description="source finding script")
    parser.add_argument("--fileName", required=True, help="the name of the input file")
    parser.add_argument("--seedSigma", type=float, default=5, help="the seed value")
    parser.add_argument("--floodfillSigma", type=float, default=3, help="the floodfill value")
    parser.add_argument("--verbose",type=bool,default=True, help="prints out stuff. Default True")
    args = parser.parse_args()

    if args.verbose:
        print("running source finding with seed snr {} and floodfill snr {}".format(args.seedSigma, args.floodfillSigma))

    main(args)


