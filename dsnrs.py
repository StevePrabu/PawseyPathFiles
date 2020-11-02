#!/usr/bin/env python
import numpy as np
from datetime import datetime, timedelta
import spacetracktool as st
from spacetracktool import operations as ops
from skyfield.api import EarthSatellite
from skyfield.api import Topos, load
import tqdm
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from skyfield.api import EarthSatellite
from skyfield.api import Topos, load
import spacetracktool as st
from os import path
import json

def obtainTLE(noradid, refUTC, args):
    """
    for a given norad id, obtains the tle for the relevant date range

    Parameters
    ----------
    noradid : the noradid of the satellite 

    Returns
    ------
    returns the tle of the satellite (tle = line1 + line2 + line3)
    """

    query = st.SpaceTrackClient(args.user, args.passwd)
    time2 = refUTC + timedelta(hours=24)
    day1, month1, year1 = str(refUTC.day).zfill(2), str(refUTC.month).zfill(2), str(refUTC.year)
    day2, month2, year2 = str(time2.day).zfill(2), str(time2.month).zfill(2), str(time2.year)
    date_range = ops.make_range_string(year1 + "-" + month1 + "-" + day1, year2 + "-" + month2 + "-" + day2)

    if path.exists(str(noradid) + ".txt") == False:



        result = query.tle_query(epoch=date_range,norad_cat_id=noradid)

        with open(str(noradid) + ".txt", "w") as outfile:
            json.dump(result.json(), outfile)

        line1 = result.json()[0]["OBJECT_NAME"]
        line2 = result.json()[0]["TLE_LINE1"]
        line3 = result.json()[0]["TLE_LINE2"]

    else:


        with open(str(noradid) + ".txt") as json_file:
            result = json.load(json_file)

        line1 = result[0]["OBJECT_NAME"]
        line2 = result[0]["TLE_LINE1"]
        line3 = result[0]["TLE_LINE2"]


    return line1 , line2, line3

def getSatMWA(line1, line2, line3):
    """
    creates a satellite object and observer object
    
    Paramters
    ---------
    line1   : tle line 1
    line2   : tle line 2
    line3   : tle line 3

    Returns
    -------
    satellite   : the sat object
    mwa         : the observer object
    ts          : the time object
    """

    ts = load.timescale()
    satellite = EarthSatellite(line2, line3, line1, ts)
    mwa = Topos("26.701276778 S", "116.670846137 E", elevation_m= 377.827)

    return satellite, mwa, ts


def main(args):

    ## get the utc times
    hdu = fits.open(args.inputFile)
    timeSteps = list(set(hdu[1].data["timeStep"]))
    wcs = WCS(hdu[1].header, naxis=2)
    utc_array = []
    imgSize = 1400
    y = np.linspace(0, (imgSize - 1), imgSize)
    x = np.linspace(0, (imgSize - 1), imgSize)
    x, y = np.meshgrid(x, y)
    imgScale = np.abs(hdu[1].header["CDELT2"])
    ## create utc array
    print("obtaining utc times")
    for t in tqdm.tqdm(timeSteps):
        
        positions = np.where(hdu[1].data["timeStep"] == t)[0][0]
        #print(positions)
        utc_array.append(hdu[1].data["utc"][positions])

        # for entry in hdu[1].data:
        #     temp_t = entry[10]
        #     if temp_t == t:
        #         utc_array.append(entry[8])
        #         break
    print("done")
    #print(utc_array)
    ## get tle for satellite
    refUTC = datetime.strptime(utc_array[0], '%Y-%m-%d %H:%M:%S')
    line1, line2, line3 = obtainTLE(args.noradid, refUTC, args)
    sat, mwa, ts = getSatMWA(line1, line2, line3)

    waterfall = np.zeros((len(timeSteps), 768))
    for utc, t in tqdm.tqdm(zip(utc_array, timeSteps)):
        print("working on timeStep {0}".format(t))
        time = datetime.strptime(utc, '%Y-%m-%d %H:%M:%S')
        ts_time = ts.utc(time.year, time.month,  time.day, time.hour, time.minute, time.second)
        sat.at(ts_time)
        difference = sat - mwa
        topocentric = difference.at(ts_time)
        ra, dec, distance = topocentric.radec()
        ra = np.degrees(ra.radians)
        dec = np.degrees(dec.radians)
        px, py = wcs.all_world2pix([[ra, dec]], 1)[0]
        LOS_range = distance.m
        radius = np.degrees(50000/LOS_range)
        number_of_pixels = radius/imgScale

        for f in range(768):

            positions = np.where(hdu[1].data["timeStep"] == t)[0]
            temp_data = np.zeros((1400,1400))
            rfi_found = False
            for p in positions:
                freq = hdu[1].data["channelNo"][p]
                if freq == f:

                    flux = hdu[1].data["fluxDensity"][p]
                    xp = int(hdu[1].data["x"][p])
                    yp = int(hdu[1].data["y"][p])
                    temp_data[yp,xp] = flux
                    rfi_found = True

            
            ## mask temp_data arround sat
            if rfi_found:
                array_radius = np.sqrt((x - xp)**2 + (y-yp)**2 )
                array_radius[array_radius > number_of_pixels] = -10
                array_radius[array_radius != -10] = 1
                array_radius[array_radius != 1] = 0
                temp_data2 = temp_data*array_radius
                waterfall[int(t),f] = np.max(temp_data2)
        plt.imshow(temp_data, origin="lower")
        plt.show()
        f_axis = np.linspace(134.215,164.335, 768 )
        
        plt.plot(f_axis,waterfall[int(t), :])
        plt.xlabel("MHz")
        plt.show()

    plt.imshow(waterfall)
    np.save(args.noradid + "waterfall.npy", waterfall)
    plt.save(args.noradid + "waterfall.png")



    



if __name__ == "__main__":
    parser = ArgumentParser("dsnrs", description="does dsnrs using vo table data")
    parser.add_argument("--noradid", required=True,type=int, help="the noradid of the satellite")
    parser.add_argument("--inputFile", required=True, help="the vo table.fits")
    parser.add_argument("--user", required=True, help="space-track user")
    parser.add_argument("--passwd", required=True, help="space-track passwd")
    args = parser.parse_args()

    main(args)

