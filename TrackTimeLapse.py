#!/usr/bin/env python
from __future__ import print_function, division
from astropy.io import fits
import numpy as np
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from skyfield.api import Topos, load, EarthSatellite
import spacetracktool as st
from spacetracktool import operations as ops
from os import path
import json
from tqdm import tqdm
from argparse import ArgumentParser
plt.style.use("dark_background")

def obtainTLE(noradid):

    time2 = startUTC + timedelta(hours=24)
    day1, month1, year1 = str(startUTC.day).zfill(2), str(startUTC.month).zfill(2), str(startUTC.year)
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





def main(args):
    ### create large wcs
    hdu_metafits = fits.open(str(args.obs) + ".metafits")
    pointing_ra = float(hdu_metafits[0].header["RA"])
    pointing_dec = float(hdu_metafits[0].header["DEC"])
    big_wcs = WCS(naxis=2)
    big_wcs.wcs.ctype = ['RA---SIN','DEC--SIN']
    big_wcs.wcs.crval = [pointing_ra, pointing_dec]
    big_wcs.wcs.crpix = [701.0, 701.0]
    big_wcs.wcs.cdelt = [-0.0833333333333333, 0.0833333333333333]
  
    ## make timeSteps by reading csv file
    
    timeSteps = np.arange(args.t1, args.t2 + 1)

    font = {'family' : 'normal','size'   : 1}
    hfont = {'fontname':'Helvetica', 'size':15, 'color': "white"}
    if debug:
        print("the selected timeSteps are " + str(timeSteps))

    line1, line2, line3 = obtainTLE(args.noradid)
    ts = load.timescale(builtin=True)
    sat = EarthSatellite(line2, line3, line1, ts)
    mwa = Topos("26.701276778 S", "116.670846137 E", elevation_m= 377.827)
    trackGlobalData = np.zeros((imgSize, imgSize))
    fullSkyGlobalData = np.zeros((1400, 1400))
    sat_x_array, sat_y_array = [], []
    for t in timeSteps:
        hdu = fits.open(args.prefix + "SigmaRFIBinaryMapSeed-t" + str(t).zfill(4) + ".fits")
        wcs = WCS(hdu[0].header, naxis=2)
        time = datetime.strptime(hdu[0].header['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f')
        ts_time = ts.utc(time.year, time.month, time.day, time.hour, time.minute, time.second + time.microsecond/1000000)
        ### determine satellite position
        sat.at(ts_time)
        difference = sat - mwa
        topocentric = difference.at(ts_time)
        ra, dec, distance = topocentric.radec()
        ra = np.degrees(ra.radians)
        dec = np.degrees(dec.radians)
        sat_x, sat_y = big_wcs.all_world2pix(np.array([[ra], [dec]]).T, 1).T
        sat_x_array.append(sat_x)
        sat_y_array.append(sat_y)



    for t in tqdm(timeSteps):
        if debug:
            print("plotting data for timeStep {}".format(t))

        hdu = fits.open(args.prefix + "SigmaRFIBinaryMapSeed-t" + str(t).zfill(4) + ".fits")
        data = hdu[0].data
        wcs = WCS(hdu[0].header, naxis=2)
        time = datetime.strptime(hdu[0].header['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f')
        ts_time = ts.utc(time.year, time.month, time.day, time.hour, time.minute, time.second + time.microsecond/1000000)
        ## determine satellite position
        #sat.at(ts_time)
        #difference = sat - mwa
        #topocentric = difference.at(ts_time)
        #ra, dec, distance = topocentric.radec()
        #ra = np.degrees(ra.radians)
        #dec = np.degrees(dec.radians)
        #sat_x, sat_y = big_wcs.all_world2pix(np.array([[ra], [dec]]).T, 1).T
        #sat_x_array.append(sat_x)
        #sat_y_array.append(sat_y)
        
     
        trackGlobalData += hdu[0].data
        wcs = WCS(hdu[0].header, naxis=2)
        fig = plt.figure(figsize = (10,5))
        st = fig.suptitle("UTC {}".format(time), fontsize="x-large")

        ax = plt.subplot(1,2,1, projection=wcs)

        plt.imshow(data, origin="lower",vmax=1,vmin=-1,cmap=plt.cm.seismic)
        plt.xlabel("RA", **hfont)
        plt.ylabel("DEC", **hfont)
        plt.grid(color="black",linestyle="dotted")

        ax = plt.subplot(1,2,2, projection=big_wcs)
        ## translate points from local fov to sky
        points = np.array(np.where(data > 0)).T
        #print(points)
        row_points, col_points = points.T
        #print(np.array([col_points, row_points]).T)
        points = np.array([col_points, row_points]).T
        
        world = wcs.wcs_pix2world(points, 0)
        px, py = big_wcs.wcs_world2pix(world,0 ).T
        for x, y in zip(px, py):
            fullSkyGlobalData[int(y), int(x)] = 1
        
        

   
        plt.imshow(fullSkyGlobalData, origin="lower", vmax=1, vmin=-1, cmap=plt.cm.seismic)
        scatter_row, scatter_col = np.where(fullSkyGlobalData == 1)
        plt.scatter(scatter_col, scatter_row, marker=".", color="red", s=1)
        plt.xlabel("RA")
        plt.ylabel("DEC")
        ### plot fov
        pixcrd = np.array([[0, 0], [0, 199], [199, 0], [199, 199]], dtype=np.float64)
        world = wcs.wcs_pix2world(pixcrd, 0)

        x, y = big_wcs.wcs_world2pix(world, 0).T
       
        #print(borders)
        plt.scatter(sat_x_array, sat_y_array, marker=".",color="green",alpha=0.8, s=0.5)
        plt.scatter(x, y, marker="x", c ="black")
        plt.grid(color="black",linestyle="dotted")
        plt.savefig("timeLapseObs" + str(args.obs)+ "norad" + str(args.noradid) + "t" + str(t).zfill(4) + ".png")

        


if __name__ == "__main__":
    parser = ArgumentParser("timelapse", description="time lapse for track pipeline")
    parser.add_argument("--obs",required=True, type=int, help="the obs id")
    parser.add_argument("--noradid",required=True, type=int, help="the noradid")
    parser.add_argument("--t1", required=True, type=int, help="the starting timeStep")
    parser.add_argument("--t2", required=True, type=int, help="the last timeStep")
    parser.add_argument("--user",required=True, help="the user name for space-track.org")
    parser.add_argument("--passwd", required=True, help="The password for space-track.org")
    parser.add_argument("--prefix", required=True, help="the file prefix")
    parser.add_argument("--debug", default=False, type=bool, help="run script in debug mode")
    args = parser.parse_args()

    global debug
    debug = args.debug

    global query
    query = st.SpaceTrackClient(args.user, args.passwd)

    ## get header info and make them global
    hdu = fits.open(str(args.prefix) + "SigmaRFIBinaryMap-t" + str(args.t1).zfill(4) + ".fits")

    global imgSize, pixel_scale,startUTC
    
    imgSize = hdu[0].header["NAXIS1"]
    pixel_scale = hdu[0].header["CDELT2"]
    startUTC = datetime.strptime(hdu[0].header['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f')

    if debug:
        print("running in debug mode")

    main(args)

