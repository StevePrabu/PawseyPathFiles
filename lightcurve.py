#!/usr/bin/env python
from __future__ import division, print_function
from argparse import ArgumentParser
from astropy.io import fits
from astropy.wcs import WCS
from datetime import datetime, timedelta
import csv
import numpy as np
import pandas as pd


def getTimeSteps(fileName):
    timeSteps = []
    utc_array = []
    cutoff_array = []
    with open(fileName) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=",")
        for row in csv_reader:
            
            timeSteps.append(int(row[0]))
            utc_array.append(row[4])
            cutoff_array.append(float(row[3]))           

    return timeSteps, utc_array, cutoff_array

def xy2RaDec(x_array, y_array, wcs):

    pixcrd = np.array([x_array, y_array], dtype=np.float64).T
    world = wcs.wcs_pix2world(pixcrd, 0)
    ra_array, dec_array = world.T
    
    return ra_array, dec_array

def getDist2Source(cutoff):
    wavelength = 3.059106714 ## lambda at 98MHz
    alt = (cutoff**2)*2/wavelength
    return alt

def floodfill(srow, scol, diff, noise, seedSigma):
    q = []
    q.append([srow, scol])
    row_array = []
    col_array = []

    while q:
        row, col = q.pop()
        row_array.append(row)
        col_array.append(col)
        tmpBinaryMap[row,col] = 1

        if diff[row+1,col] >= noise*seedSigma and tmpBinaryMap[row+1, col] == 0:
            q.append([row+1, col])
        if diff[row,col+1] >= noise*seedSigma and tmpBinaryMap[row, col+1] == 0:
            q.append([row, col+1])
        if diff[row-1,col] >= noise*seedSigma and tmpBinaryMap[row-1, col] == 0:
            q.append([row-1, col])
        if diff[row,col-1] >= noise*seedSigma and tmpBinaryMap[row, col-1] == 0:
            q.append([row, col-1])

        if diff[row+1,col+1] >= noise*seedSigma and tmpBinaryMap[row+1, col+1] == 0:
            q.append([row+1, col+1])
        if diff[row+1,col-1] >= noise*seedSigma and tmpBinaryMap[row+1, col-1] == 0:
            q.append([row+1, col-1])
        if diff[row-1,col-1] >= noise*seedSigma and tmpBinaryMap[row-1, col-1] == 0:
            q.append([row-1, col-1])
        if diff[row-1,col+1] >= noise*seedSigma and tmpBinaryMap[row-1, col+1] == 0:
            q.append([row-1, col+1])

    return row_array, col_array

        

def main(args):

    timeSteps, utc_array, cutoff_array = getTimeSteps(args.timeStepFile)

    ## load diff map
    df = pd.read_pickle(args.freqDiffMap)
    
    ## initialse output arrays
    out_timestep = []
    out_utc = []
    out_channel = []
    out_freq = []
    out_channel2 = []
    out_channelDiff = []
    out_valueType = [] ## one of the choice ['detection', 'upper limit', 'flagged']
    out_value = []
    out_err = []
    out_noise = []
    out_deltaX = []
    out_deltaY = []
    out_deltaRA = []
    out_deltaDEC = []
    global tmpBinaryMap
    
    for t, utc, cutoff in zip(timeSteps, utc_array, cutoff_array):
        print("working on timestep {}".format(t))
        createMask = False
        
        for f in range(args.noChannels):
            tmpBinaryMap = np.zeros((200,200))

            hdu = fits.open("{}-2m-{}-{}-image.fits".format(args.obs, t, str(f).zfill(4)))
            f2diff = int(df['diffChannelIndex'][df['mwaChannelIndex']==f])
            hdu_diff = fits.open("{}-2m-{}-{}-image.fits".format(args.obs, t, str(f2diff).zfill(4)))
            data = hdu[0].data[0,0,:,:]
            data_diff = hdu_diff[0].data[0,0,:,:]
            out_timestep.append(t)
            out_utc.append(utc)
            out_channel.append(f)
            out_freq.append(hdu[0].header['CRVAL3']/10**6)
            out_channel2.append(f2diff)
            out_channelDiff.append(np.abs(f2diff-f))

            ## create mask
            if not createMask:
                ## mask source
                imgSize = hdu[0].header["NAXIS1"]
                imgScale = np.abs(hdu[0].header["CDELT2"])
                pix2deg = abs(imgScale)

                
                ## mask source
                LOS_range = getDist2Source(cutoff)
                radius = np.degrees(20000/LOS_range) # cone radius (s = r * theta)
                number_of_pixels = radius/imgScale
                y = np.linspace(0, (imgSize - 1), imgSize)
                x = np.linspace(0, (imgSize - 1), imgSize)
                x, y = np.meshgrid(x, y)

                array_radius = np.sqrt((x - 99.5)**2 + (y-99.5)**2 )
                array_radius[array_radius > number_of_pixels] = -10
                array_radius[array_radius != -10] = 1
                array_radius[array_radius != 1] = 0
                mask = array_radius
                createMask = True

            if np.all(data == 0) or np.all(data_diff==0):   
                
                out_valueType.append("flagged")
                out_value.append(np.nan)
                out_err.append(np.nan)
                out_noise.append(np.nan)
                out_deltaX.append(np.nan)
                out_deltaY.append(np.nan)
                out_deltaRA.append(np.nan)
                out_deltaDEC.append(np.nan)
                continue


            diff = data - data_diff
            tmp = np.copy(diff)
            noise = np.std(tmp)
            diff = diff*mask
            out_noise.append(noise)  

            maxSNR = np.nanmax(diff)/noise

            threshold = 6
            if maxSNR >= threshold:
                bmin = hdu[0].header['BMIN']
                bmaj = hdu[0].header['BMAJ']
                wcs = WCS(hdu[0].header, naxis=2)
                
                row, col = np.where(diff == np.nanmax(diff))
                row, col = int(row), int(col)

                ## floodfill and find source
                row_array, col_array = floodfill(row, col, diff, noise, 3)

                ra_cent , dec_cent = xy2RaDec([99.5], [99.5], wcs)
                ra_peak , dec_peak = xy2RaDec([np.mean(col_array)], [np.mean(row_array)], wcs)
                deltaRA = np.abs(ra_cent - ra_peak)
                deltaDEC = np.abs(dec_cent - dec_peak)
                deltaX = np.abs(np.mean(col_array) - 99.5)
                deltaY = np.abs(np.mean(row_array) - 99.5)
                
                ## calculate integrated flux density
                beamVolume = 1.1331*bmin*bmaj
                intflux = np.sum(diff*tmpBinaryMap)*pix2deg**2/beamVolume
                npix = np.sum(tmpBinaryMap)
                nbeams = npix*(pix2deg**2)/beamVolume
                err = np.sqrt(nbeams)*noise

                print("found measurement at t {} f {} snr {}".format(t,f,maxSNR))

                out_valueType.append("detection")
                out_value.append(intflux)
                out_err.append(err)
                out_deltaX.append(deltaX)
                out_deltaY.append(deltaY)
                out_deltaRA.append(deltaRA)
                out_deltaDEC.append(deltaDEC)                

            else:
                out_valueType.append("upper limit")
                out_value.append(noise*threshold)
                out_err.append(np.nan)
                out_deltaX.append(np.nan)
                out_deltaY.append(np.nan)
                out_deltaRA.append(np.nan)
                out_deltaDEC.append(np.nan)

    ## write output df
    
    d = {'timestep': out_timestep,
        'utc': out_utc,
        'channel': out_channel,
        'freq': out_freq,
        'valueType': out_valueType,
        'value': out_value,
        'err': out_err,
        'noise': out_noise,
        'deltaX': out_deltaX,
        'deltaY': out_deltaY,
        'deltaRA': out_deltaRA,
        'deltaDEC': out_deltaDEC,
        'channel2': out_channel2,
        'channelDiff': out_channelDiff}
        
    df_out = pd.DataFrame(data=d)
    df_out.to_pickle("lc_{}_{}.plk".format(args.obs, args.noradid))
                



if __name__ == "__main__":
    parser = ArgumentParser("lightcurve", description="extracts light curve")
    parser.add_argument("--obs",required=True, type=int, help="the obs id")
    parser.add_argument("--timeStepFile", required=True, help="output file from ligth_curve_track.py")
    parser.add_argument("--noChannels",default=768,type=int, help="the number of fine channels")
    parser.add_argument("--noradid", required=True, type=int, help="the norad id")
    parser.add_argument("--tleCatalog", required=True, help="the input tle catalog")
    parser.add_argument("--freqDiffMap", default="freqDiffMap.plk", help="plk file that contains the channels to diff")
    args = parser.parse_args()

    main(args)

    