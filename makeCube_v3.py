#!/usr/bin/env python
from tqdm import tqdm
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

def obtainTLE(noradid,refUTC):
    debug = True
    time1 = refUTC + timedelta(hours=-24*2)
    time2 = refUTC 
    day1, month1, year1 = str(time1.day).zfill(2), str(time1.month).zfill(2), str(time1.year)
    day2, month2, year2 = str(time2.day).zfill(2), str(time2.month).zfill(2), str(time2.year)
    date_range = ops.make_range_string(year1 + "-" + month1 + "-" + day1, year2 + "-" + month2 + "-" + day2)

    if path.exists(str(noradid) + ".txt") == False:

        if debug:
            print("requesting file from server")

        result = query.tle_query(epoch=date_range,norad_cat_id=noradid)

        with open(str(noradid) + ".txt", "w") as outfile:
            json.dump(result.json(), outfile)

        line1 = result.json()[0]["OBJECT_NAME"]
        line2 = result.json()[0]["TLE_LINE1"]
        line3 = result.json()[0]["TLE_LINE2"]

    else:

        if debug:
            print("tle file found on disk. Not downloading.")

        with open(str(noradid) + ".txt") as json_file:
            result = json.load(json_file)

        line1 = result[0]["OBJECT_NAME"]
        line2 = result[0]["TLE_LINE1"]
        line3 = result[0]["TLE_LINE2"]


    return line1 , line2, line3

def getSat(line1, line2, line3):

    return EarthSatellite(line2, line3, line1, ts)


def getSatXY(sat, local_utc, wcs):
    local_ts = ts.utc(local_utc.year, local_utc.month, local_utc.day, local_utc.hour, local_utc.minute, local_utc.second + local_utc.microsecond/1000000)
    sat.at(local_ts)
    difference = sat - mwa
    topocentric = difference.at(local_ts)
    ra, dec, distance = topocentric.radec()
    ra, dec = np.degrees(ra.radians), np.degrees(dec.radians)
    px, py = wcs.all_world2pix([[ra, dec]], 1)[0]
    return px, py

def getSlope(x_array, y_array):
    slope = np.rad2deg(np.arctan2(y_array[1]- y_array[0], x_array[1] - x_array[0]))
    return slope - 90

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
    
        

def main(args):
    
    ## make timeStep array by reading csv file
    timeSteps = []
    with open(str(args.obs) + "-" + str(args.noradid) + ".csv") as csv_file:
         csv_reader = csv.reader(csv_file, delimiter=",")
         for row in csv_reader:
             timeSteps.append(int(row[0]))
    
    cube = []
    
    ### get tle
    hdu = fits.open(str(args.obs)+ str(args.band) + "-" + str(args.midName) + "-" + str(timeSteps[0]) + "h-" + str(0).zfill(4) + "-dirty.fits" )
    ref_utc = datetime.strptime(hdu[0].header["DATE-OBS"], '%Y-%m-%dT%H:%M:%S.%f')
    line1, line2, line3 = obtainTLE(args.noradid,ref_utc)
    sat = getSat(line1, line2, line3)
    
    
    ## get rotation array
    print("calculating rotations...")
    rotation_array = []
    for t in timeSteps:
        hduH = fits.open(str(args.obs)+ str(args.band) + "-" + str(args.midName) + "-" + str(t) + "h-" + str(0).zfill(4) + "-dirty.fits" )
        utc = datetime.strptime(hduH[0].header["DATE-OBS"], '%Y-%m-%dT%H:%M:%S.%f')
        wcs = WCS(hduH[0].header, naxis=2)
        slope = getRotation(sat, utc, wcs)
        rotation_array.append(slope)
    print("done")
    
    for f in tqdm(range(args.channels)):
        global_data = []
        w = []
        for t, slope in zip(timeSteps,rotation_array):
            hduH = fits.open(str(args.obs)+ str(args.band) + "-" + str(args.midName) + "-" + str(t) + "h-" + str(f).zfill(4) + "-dirty.fits" )
            hduT = fits.open(str(args.obs)+ str(args.band) + "-" + str(args.midName) + "-" + str(t) + "t-" + str(f).zfill(4) + "-dirty.fits" )
            diff = hduH[0].data[0,0,:,:] - hduT[0].data[0,0,:,:] 
        
            ## the below logic is to avoid getting nans due to inverting noise = 0
            if np.std(diff) == 0:
                continue
            else:
                diff = rotate(diff, slope, order=5, reshape=False)
                global_data.append(diff)
                w.append(np.std(diff))

        if not w:
            cube.append(np.zeros(diff.shape))
        else:
            w = np.array(w)
            weights = 1/w
            stack = np.average(np.array(global_data), axis=0, weights=weights)
            cube.append(stack)
            
    np.save("weightedRotated"+ str(args.noradid) + "-" + str(args.obs) + ".npy", cube)

    ## make images of all 6sigma events
    cube = np.array(cube)

    for f in range(cube.shape[0]):
        temp1 = np.copy(cube[f,:,:])
        temp2 = np.copy(cube[f,:,:])
        signal = np.nanmax(temp1)
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
    parser = ArgumentParser("make cube", description="stackes images and makes a cube for later analysis")
    parser.add_argument("--obs", required=True,  help="the observation id")
    parser.add_argument("--band", default="", help="the band name")
    parser.add_argument("--noradid",required=True, type=int, help="the norad id")
    parser.add_argument("--channels",default=768,type=int, help="the number of channels")
    parser.add_argument("--midName", default="2m", help="the mid name of the fits files")
    parser.add_argument("--user",required=True, help="spacetrack.org user name")
    parser.add_argument("--passwd",required=True, help="spaceTrack.org password")
    args = parser.parse_args()
    
    global query, ts, mwa
    query = st.SpaceTrackClient(args.user, args.passwd)
    ts = load.timescale(builtin=True)
    mwa = Topos("26.701276778 S", "116.670846137 E", elevation_m= 377.827)

    main(args)
    
    
