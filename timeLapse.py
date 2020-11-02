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

def getSatXY(line1, line2, line3, time, mwa, ts):
    
    sat = EarthSatellite(line2, line3, line1, ts)
    sat.at(time)
    difference = sat - mwa
    topocentric = difference.at(time)
    ra, dec, distance = topocentric.radec()
    ra, dec = np.degrees(ra.radians), np.degrees(dec.radians)
    px, py = wcs.all_world2pix([[ra, dec]], 1)[0]

    return px, py, distance


def getTLE(UTC):

    time1 = UTC + timedelta(hours=-24*args.tlepast)
    time2 = UTC + timedelta(hours=24*args.tlefuture)
    day1, month1, year1 = str(time1.day).zfill(2), str(time1.month).zfill(2), str(time1.year)
    day2, month2, year2 = str(time2.day).zfill(2), str(time2.month).zfill(2), str(time2.year)
    custom_name = year1 + "-" + month1 + "-" + day1 + "__" +  year2 + "-" + month2 + "-" + day2
    entries = 0

    if path.exists("../TLE_catalog" + custom_name + ".txt") == False:
        
        if debug:
            print("requesting file from server...",end="")

        date_range = ops.make_range_string(year1 + "-" + month1 + "-" + day1, year2 + "-" + month2 + "-" + day2)
        result = query.tle_publish_query(publish_epoch=date_range)

        if debug:
            print("done")
        
        ## write the data to file for latter use
        with open("../TLE_catalog" + custom_name + ".txt", "w") as outfile:
            json.dump(result.json(), outfile)
    
        entries = len(result.json()[:])
        catalog = result.json()

    else:

        if debug:
            print("tle file found on disk. Not downloading.")
            print("using file ../" + str(custom_name) + ".txt")
        
        with open("../TLE_catalog" + custom_name + ".txt") as json_file:
            result = json.load(json_file)
        
        entries = len(result[:])
        catalog = result

    
    return catalog, entries


def main(args):

    timeSteps = np.arange(args.t1, args.t2 + 1)
    font = {'family' : 'normal','size'   : 1}

    if debug:
        print("the selected timeSteps are " + str(timeSteps))

    ## get tle catalog
    catalog, noObjects = getTLE(startUTC)

    if debug:
        print("obtained tle for {0} objects".format(noObjects))

    ts = load.timescale()
    mwa = Topos("26.703319405555554 S", "116.91558083333334 E", elevation_m= 377.827)
    globalData = np.zeros((imgSize, imgSize))
    globalData2 = np.zeros((imgSize, imgSize))
    track_grid = np.zeros((imgSize, imgSize))

    for t in timeSteps:

        if debug:
            print("searching for sats in timeStep {0}".format(t))
        
        hdu_seed = fits.open(str(args.prefix) + "SigmaRFIBinaryMapSeed-t" + str(t).zfill(4) + ".fits")
        hdu_binary = fits.open(str(args.prefix) + "SigmaRFIBinaryMap-t" + str(t).zfill(4) + ".fits")
        hdu_peakFlux = fits.open(str(args.prefix) + "SigmaRFIBinaryMapPeakFlux-t" + str(t).zfill(4) + ".fits")

        time = datetime.strptime(hdu_seed[0].header['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f')
        ts_time = ts.utc(time.year, time.month, time.day, time.hour, time.minute, time.second)
        data_seed = hdu_seed[0].data
        data_binary = hdu_binary[0].data

        globalData += data_binary
        globalData += hdu_peakFlux[0].data

        ## start plot
        ax = plt.subplot(1,1,1, projection=wcs)
        plt.imshow(globalData, origin="lower", cmap=plt.cm.seismic, vmax=1, vmin=-1)
        plt.grid()
        plt.title("UTC " + str(time))
        ax.set_xlabel("RA (Degrees)")
        ax.set_ylabel("DEC (Degrees)")

        sats_searched = []

        for satNo in tqdm(range(noObjects)):

            line1 = "sat"
            line2 = catalog[satNo]["TLE_LINE1"]
            line3 = catalog[satNo]["TLE_LINE2"]
            norad = line2[2:7]            

            if norad in sats_searched:
                continue

            else:

                sats_searched.append(norad)

                px, py, distance = getSatXY(line1, line2, line3, ts_time, mwa, ts)

                if 0 < px < imgSize and 0 < py < imgSize:
                    track_grid[int(py), int(px)] = 1
                    #plt.scatter(px, py, marker=".", color="black", s=0.5)
                    plt.text(px, py, str(norad) + "d"+ str(distance), color="black", **font)

        points = np.asarray(np.where(track_grid == 1))
        ax.scatter(points[1,:], points[0,:], marker="o", c="blue", s=0.01)
        plt.savefig(str(args.obs)+"timeLapseHD" + str(t).zfill(4) + ".png", dpi=800)
        plt.savefig(str(args.obs)+"timeLapseSD" + str(t).zfill(4) + ".png")
        plt.clf()

        ax = plt.subplot(1,1,1, projection=wcs)
        plt.imshow(globalData2, origin="lower", cmap=plt.cm.seismic)
        plt.colorbar()
        plt.grid()
        plt.title("UTC " + str(time))
        ax.set_xlabel("RA (Degrees)")
        ax.set_ylabel("DEC (Degrees)")

        plt.savefig(str(args.obs)+"timeLapseSDFluxMap" + str(t).zfill(4) + ".png")
        plt.clf()





if __name__ == "__main__":
    parser = ArgumentParser("combined plot", description="used for identifying satellites")
    parser.add_argument("--obs",required=True, type=int, help="the observation id")
    parser.add_argument("--t1", required=True, type=int, help="the starting timeStep")
    parser.add_argument("--t2", required=True, type=int, help="the last timeStep")
    parser.add_argument("--user",required=True, help="the user name for space-track.org")
    parser.add_argument("--passwd", required=True, help="The password for space-track.org")
    parser.add_argument("--prefix", required=True, help="the file prefix")
    parser.add_argument("--tlepast", default=1, type=int, help="number of days from past to search tle from")
    parser.add_argument("--tlefuture", default=0, type=int, help="number of days from future to search tle from")
    parser.add_argument("--debug", default=False, type=bool, help="run script in debug mode")
    args = parser.parse_args()

    global debug
    debug = args.debug

    global query
    query = st.SpaceTrackClient(args.user, args.passwd)

    ## get header info and make them global
    hdu = fits.open(str(args.prefix) + "SigmaRFIBinaryMap-t" + str(args.t1).zfill(4) + ".fits")

    global wcs, imgSize, pixel_scale
    wcs = WCS(hdu[0].header, naxis=2)
    imgSize = hdu[0].header["NAXIS1"]
    pixel_scale = hdu[0].header["CDELT2"]
    startUTC = datetime.strptime(hdu[0].header['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f')

    if debug:
        print("running in debug mode")

    main(args)
