#!/usr/bin/env python
from __future__ import division, print_function
from astropy.io import fits
import numpy as np
from argparse import ArgumentParser
import csv
from scipy.ndimage import  rotate
from datetime import datetime, timedelta
import spacetracktool as st
from spacetracktool import operations as ops
from skyfield.api import EarthSatellite
from skyfield.api import Topos, load
import json
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from os import path
import pandas as pd
import sys


def obtainTLE(noradid, obs):
    debug = True
    
    if path.exists(str(noradid) + str(obs) + ".txt") == False:

        print("tle file not found. Exiting...")
        sys.exit()

    else:

        if debug:
            print("tle file found on disk.")

        with open(str(noradid) + str(obs) + ".txt") as json_file:
            result = json.load(json_file)

        line1 = result[0]["OBJECT_NAME"]
        line2 = result[0]["TLE_LINE1"]
        line3 = result[0]["TLE_LINE2"]


    return line1 , line2, line3

def getRotation(sat, utc, wcs):
    time_array = [-4 , 4]
    x_array, y_array = [], []
    for t in time_array:
        local_utc = utc + timedelta(seconds=t)
        x, y = getSatXY(sat, local_utc, wcs)
        x_array.append(x)
        y_array.append(y)
    
    slope = getSlope(x_array, y_array)
    
    return slope
    
def getSatXY(sat, local_utc, wcs):
    local_ts = ts.utc(local_utc.year, local_utc.month, local_utc.day, local_utc.hour, local_utc.minute, local_utc.second + local_utc.microsecond/1000000)
    sat.at(local_ts)
    difference = sat - mwa
    topocentric = difference.at(local_ts)
    ra, dec, distance = topocentric.radec()
    ra, dec = np.degrees(ra.radians), np.degrees(dec.radians)
    px, py = wcs.all_world2pix([[ra, dec]], 1)[0]
    return px, py

def getSat(line1, line2, line3):

    return EarthSatellite(line2, line3, line1, ts)


def main(args):

    ### make timestep array by reading csv file
    timeSteps = []
    utc_array = []
    with open(str(args.obs) + "-" + str(args.noradid) + ".csv") as csv_file:
         csv_reader = csv.reader(csv_file, delimiter=",")
         for row in csv_reader:
             timeSteps.append(int(row[0]))
             utc_array.append(datetime.strptime(row[4], '%Y-%m-%dT%H:%M:%S.%f'))
    
    cube = []

    ## get tle
    line1, line2, line3 = obtainTLE(args.noradid, args.obs)
    sat = getSat(line1, line2, line3)

    ## load diff map
    df = pd.read_pickle(args.freqDiffMap)

    for f in range(args.channels):
        global_data = []
        w = []
        f2diff = int(df['diffChannelIndex'][df['mwaChannelIndex']==f])
        for t in timeSteps:
            hdu1 = fits.open('img-t{}-f{}.fits'.format(t, f))
            hdu2 = fits.open('img-t{}-f{}.fits'.format(t, f2diff))
            #diff = hdu1[0].data - hdu2[0].data
            diff = hdu1[0].data 

            ## the below logic is to avoid getting nans due to inverting noise = 0
            if np.std(diff) == 0:
                continue
            else:
                global_data.append(diff)
                w.append(np.std(diff))

        if not w:
            cube.append(np.zeros(diff.shape))
        else:
            w = np.array(w)
            weights = 1/w
            try:
                stack = np.average(np.array(global_data), axis=0, weights=weights)
                cube.append(stack)
            except:
                cube.append(np.zeros(diff.shape))

        np.save("weightedRotated"+ str(args.noradid) + "-" + str(args.obs) + ".npy", cube)

    ## make images of all 6sigma events
    cube = np.array(cube)

    ### check if cube is full of nans, and if so exit(1)
    if np.all(np.isnan(cube)):
        print("cube full of zeros. terminating..")
        sys.exit(1)

    for f in range(cube.shape[0]):
        temp1 = np.copy(cube[f,:,:])
        temp2 = np.copy(cube[f,:,:])
        signal = np.nanmax(np.abs(temp1))
        temp2[np.abs(temp2) > 3*np.std(temp2)] = 0
        temp2[np.abs(temp2) > 3*np.std(temp2)] = 0

        noise = np.std(temp2)
        snr = signal/noise

        if snr >= 6:
            plt.clf()
            plt.imshow(cube[f,:,:], origin="lower")
            plt.colorbar()
            plt.title("channel {}".format(f))
            plt.savefig("weightedEvent" + str(f).zfill(4) + ".png")
    


if __name__ == "__main__":
    parser = ArgumentParser('make cube', description='stackes images and makes a cube for later analysis')
    parser.add_argument('--obs', required=True, help='the observation id')
    parser.add_argument('--noradid', required=True, type=int, help='the norad id')
    parser.add_argument('--channels', default=768, type=int, help='the number of channels')
    parser.add_argument('--user', required=True, help='spacetrack.org user name')
    parser.add_argument('--passwd', required=True, help='spacetrack.org password')
    parser.add_argument('--freqDiffMap', required=True, help='plk object that has info about freq diff')
    args = parser.parse_args()

    global query, ts, mwa
    query = st.SpaceTrackClient(args.user, args.passwd)
    ts = load.timescale(builtin=True)
    mwa = Topos("26.701276778 S", "116.670846137 E", elevation_m= 377.827)

    main(args)