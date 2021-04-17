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
import csv

def floodfill(xs, ys, floodfillValue, data, noise, imgSize):
    """
    performs the forest fire algorithm to source find
    """
    
    q = []
    q.append([xs, ys])

    while q:

        x, y = q.pop()
        searchMap[x,y] = 1
        eventMask[x,y] = 1
        snrMap[x,y] = data[x,y]/noise[x,y]

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

def xy2RaDec(x_array, y_array, wcs):
    """
    converts pixel coords to ra dec
    Parameters
    ----------
    x   : x pixel coord
    y   : y pixel coord
    wcs : world coord system obj
    Returns
    -------
    ra  : ra in degrees
    dec : dec in degrees
    """
    pixcrd = np.array([x_array, y_array], dtype=np.float64).T
    world = wcs.wcs_pix2world(pixcrd, 0)
    #print(world)
    ra_array, dec_array = world.T
    
    return ra_array, dec_array


def checkValid(x,y,floodfillValue, data, noise, imgSize):
    """
    checks validity of seed pixel
    """

    if 1 < x < (imgSize[0]-1) and 1 < y < (imgSize[1]-1) and \
    data[x,y]/noise[x,y] >= floodfillValue and searchMap[x,y] == 0:
        return True
    else:
        return False


def calculateNoise(data, args, unit, imgSize):
    """
    function calculates noise
    if input noise map provided, uses proved map
    else, assumes uniform noise and calculates the noise
    through recursive masking
    """

    ### check if noise map provided
    if args.noiseMap is None:

        if args.verbose:
            print("noise not map provided. Estimating global uniform noise")

        ## calculate noise through recursive masking
        tmp = np.copy(data)
        tmp[np.abs(tmp) >3*np.std(tmp)] = 0
        tmp[np.abs(tmp) >3*np.std(tmp)] = 0
        tmp[np.abs(tmp) >3*np.std(tmp)] = 0

        if args.verbose:
            print("estimated global noise in image {} {}".format(np.std(tmp), unit))

        return np.std(tmp)*np.ones(imgSize)

    else:

        if args.verbose:
            print("using provided noise map file {}".args.noiseMap)

        hdu = fits.open(args.noiseMap)
        noise = hdu[0].data[0,0,:,:] 

        ## check if image size and noise map size are equal
        if imgSize != noise.shape:
            print("input noise map and image are not of the same dimensions")
            print("terminating.")
            exit()

        return noise



def main(args):
    
    ## open file and extract all parameters from header
    hdu = fits.open(args.fileName)
    data = hdu[0].data[0,0,:,:]
    wcs = WCS(hdu[0].header, naxis=2)
    unit = hdu[0].header['BUNIT']
    imgSize = data.shape
    noise = calculateNoise(data, args, unit, imgSize)
    
    ## create global SNR map
    global snrMap, searchMap, eventMask
    snrMap = np.zeros(imgSize)
    searchMap = np.zeros(imgSize)
    eventMask = np.zeros(imgSize)
    
    ## start source finding
    seeds = np.asarray(np.where(data >= args.seedSigma*noise)).T
    
    counter = 0
    x_array, y_array = [], [] ## store x,y pixel locations of all seed events
    snr_array, peakFlux_array = [], [] ## store snr and peak flux for all detected seed events

    for seed in tqdm(seeds):

        if searchMap[seed[0],seed[1]] == 1:
            continue
        eventMask = np.zeros(imgSize)
        counter += 1
        floodfill(seed[0], seed[1], args.floodfillSigma, data, noise, imgSize)
        x_array.append(seed[0])
        y_array.append(seed[1])
        snr_array.append(np.nanmax(eventMask*data/noise))
        peakFlux_array.append(np.nanmax(eventMask*data))

    if args.verbose:
        print("{} seed events found".format(counter))

    ## save snrmap to fits file
    if os.path.isfile("snr-{}".format(args.fileName)):
        bashSyn1 = "rm snr-{}".format(args.fileName)
        bashExecute1 = call(bashSyn1,shell=True)
    hdu_ouput = fits.PrimaryHDU(snrMap,header=hdu[0].header)
    hdu_ouput.writeto("snr-{}".format(args.fileName))

    ## write dected events to file
    ra_array, dec_array = xy2RaDec(x_array, y_array, wcs)

    if args.verbose:
        print("writing data to file...",end="")

    with open("eventCatalog.csv", "w") as vsc:
        thewriter = csv.writer(vsc)
        thewriter.writerow(["ra", "dec", "x", "y", "peakFlux", "snr"])
        for ra, dec, x, y, peakFlux, snr in zip(ra_array, dec_array, x_array, y_array, peakFlux_array, snr_array):
            line = [ra, dec, x, y, peakFlux, snr]
            thewriter.writerow(line)

    if args.verbose:
        print("done")

    
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
    parser.add_argument("--noiseMap", help="[Optional] input noise map. if not provided, assumes uniform noise")
    parser.add_argument("--verbose",type=bool,default=True, help="prints out stuff. Default True")
    args = parser.parse_args()

    if args.verbose:
        print("running source finding with seed snr {} and floodfill snr {}".format(args.seedSigma, args.floodfillSigma))

    main(args)


