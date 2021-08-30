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

def obtainTLE(noradid,refUTC, obs):
    debug = True
    time1 = refUTC + timedelta(hours=-24*2)
    time2 = refUTC 
    day1, month1, year1 = str(time1.day).zfill(2), str(time1.month).zfill(2), str(time1.year)
    day2, month2, year2 = str(time2.day).zfill(2), str(time2.month).zfill(2), str(time2.year)
    date_range = ops.make_range_string(year1 + "-" + month1 + "-" + day1, year2 + "-" + month2 + "-" + day2)

    if path.exists(str(noradid) + str(obs) + ".txt") == False:

        if debug:
            print("requesting file from server")

        result = query.tle_query(epoch=date_range,norad_cat_id=noradid)

        with open(str(noradid) + str(obs) + ".txt", "w") as outfile:
            json.dump(result.json(), outfile)

        line1 = result.json()[0]["OBJECT_NAME"]
        line2 = result.json()[0]["TLE_LINE1"]
        line3 = result.json()[0]["TLE_LINE2"]

    else:

        if debug:
            print("tle file found on disk. Not downloading.")

        with open(str(noradid) + str(obs) + ".txt") as json_file:
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
    hdu = fits.open(str(args.obs)+ str(args.band) + "-" + str(args.midName) + "-" + str(timeSteps[0]) + "h-" + str(args.channel1) + "-dirty.fits" )
    ref_utc = datetime.strptime(hdu[0].header["DATE-OBS"], '%Y-%m-%dT%H:%M:%S.%f')
    line1, line2, line3 = obtainTLE(args.noradid,ref_utc, args.obs)
    sat = getSat(line1, line2, line3)

    ## get rotation array
    rotation_array = []
    for t in timeSteps:
        hduH = fits.open(str(args.obs)+ str(args.band) + "-" + str(args.midName) + "-" + str(t) + "h-" + str(args.channel1) + "-dirty.fits" )
        utc = datetime.strptime(hduH[0].header["DATE-OBS"], '%Y-%m-%dT%H:%M:%S.%f')
        wcs = WCS(hduH[0].header, naxis=2)
        slope = getRotation(sat, utc, wcs)
        rotation_array.append(slope)
    print("done")

    cube1, cube2 = [], []
    snr1_array = []
    snr2_array = []
    
    w1l = []
    w2l = []
    counter = 1
    for t, slope in zip(timeSteps,rotation_array):
        hduH = fits.open(str(args.obs)+ str(args.band) + "-" + str(args.midName) + "-" + str(t) + "h-" + str(args.channel1) + "-dirty.fits" )
        hduT = fits.open(str(args.obs)+ str(args.band) + "-" + str(args.midName) + "-" + str(t) + "t-" + str(args.channel1) + "-dirty.fits" )
        diff = hduH[0].data[0,0,:,:] - hduT[0].data[0,0,:,:] 

        ## the below logic is to avoid getting nans due to inverting noise = 0
        if np.std(diff) == 0:
            continue
        else:
            diff = rotate(diff, slope, order=5, reshape=False)
            cube1.append(diff)
            w1l.append(np.std(diff))
            
        
        ## make plot
        fig, ax = plt.subplots(figsize=(10,15))
        plt.subplot(321)
        plt.imshow(diff, origin="lower")
        plt.colorbar()
        plt.title("timestep {}/{} channel {}".format(counter, len(timeSteps), args.channel1))

        
        w1 = np.array(w1l)
        weights1 = 1/w1
        stack1 = np.average(np.array(cube1), axis=0, weights=weights1)
        plt.subplot(323)
        plt.title("stacked img for channel {}".format(args.channel1))
        plt.imshow(stack1, origin="lower")
        plt.colorbar()

        tmp1 = np.copy(stack1)
        tmp1[np.abs(tmp1) > 3*np.std(tmp1)] = 0
        tmp1[np.abs(tmp1) > 3*np.std(tmp1)] = 0
        snr1_array.append(np.nanmax(stack1)/np.std(tmp1))
        plt.subplot(325)
        plt.ylabel("stacked snr")
        plt.plot(snr1_array)
        plt.xlabel("timeSteps")

        ## do the above for channel 2
        hduH = fits.open(str(args.obs)+ str(args.band) + "-" + str(args.midName) + "-" + str(t) + "h-" + str(args.channel2) + "-dirty.fits" )
        hduT = fits.open(str(args.obs)+ str(args.band) + "-" + str(args.midName) + "-" + str(t) + "t-" + str(args.channel2) + "-dirty.fits" )
        diff = hduH[0].data[0,0,:,:] - hduT[0].data[0,0,:,:] 

        ## the below logic is to avoid getting nans due to inverting noise = 0
        if np.std(diff) == 0:
            continue
        else:
            diff = rotate(diff, slope, order=5, reshape=False)
            cube2.append(diff)
            w2l.append(np.std(diff))
            
        
        ## make plot
        plt.subplot(322)
        fig.suptitle('obs {} norad {}'.format(args.obs, args.noradid))
        plt.imshow(diff, origin="lower")
        plt.colorbar()
        plt.title("timestep {}/{} channel {}".format(counter, len(timeSteps), args.channel2))

        
        w2 = np.array(w2l)
        weights2 = 1/w2
        stack2 = np.average(np.array(cube2), axis=0, weights=weights2)
        plt.subplot(324)
        plt.title("stacked img for channel {}".format(args.channel2))
        plt.imshow(stack2, origin="lower")
        plt.colorbar()

        tmp2 = np.copy(stack2)
        tmp2[np.abs(tmp2) > 3*np.std(tmp2)] = 0
        tmp2[np.abs(tmp2) > 3*np.std(tmp2)] = 0
        snr2_array.append(np.nanmax(stack2)/np.std(tmp2))
        plt.subplot(326)
        plt.ylabel("stacked snr")
        plt.plot(snr2_array)
        plt.xlabel("timeSteps")






        plt.savefig("fineChannelImage" + str(counter).zfill(4) + ".png")
        plt.clf()
        plt.close()


        


        counter += 1
        


    



if __name__ == "__main__":
    parser = ArgumentParser("fine channel investigation", description="makes fine channel image inspection pngs")
    parser.add_argument("--obs", required=True, type=int, help="the observation id")
    parser.add_argument("--noradid",required=True, type=int, help="the norad id")
    parser.add_argument("--channel1", type=int, help="the brightest channel")
    parser.add_argument("--channel2", type=int, help="the second brightest channel")
    parser.add_argument("--user", required=True, help="spacetrack.org user name")
    parser.add_argument("--passwd", required=True, help="spacetrac.org password")
    parser.add_argument("--midName", default="2m", help="the mid name of the fits files")
    parser.add_argument("--band", default="", help="the band name")
    args = parser.parse_args()

    global query, ts, mwa
    query = st.SpaceTrackClient(args.user, args.passwd)
    ts = load.timescale(builtin=True)
    mwa = Topos("26.701276778 S", "116.670846137 E", elevation_m= 377.827)

    main(args)
