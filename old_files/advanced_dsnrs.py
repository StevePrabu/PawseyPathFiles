#!/usr/bin/env python
from __future__ import print_function, division
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from datetime import datetime, timedelta
from argparse import ArgumentParser
import spacetracktool as st
from spacetracktool import operations as ops
from skyfield.api import EarthSatellite
from skyfield.api import Topos, load
from tqdm import tqdm
from astropy.nddata import Cutout2D
import matplotlib.pyplot as plt
from astropy.coordinates import EarthLocation
from datetime import datetime, timedelta
import tqdm
from os import path
import json


def obtainTLE(noradid, refUTC):
    """
    for a given norad id, obtains the tle for the relevant date range

    Parameters
    ----------
    noradid : the noradid of the satellite 

    Returns
    ------
    returns the tle of the satellite (tle = line1 + line2 + line3)
    """

    time2 = refUTC + timedelta(hours=24)
    day1, month1, year1 = str(refUTC.day).zfill(2), str(refUTC.month).zfill(2), str(refUTC.year)
    day2, month2, year2 = str(time2.day).zfill(2), str(time2.month).zfill(2), str(time2.year)
    date_range = ops.make_range_string(year1 + "-" + month1 + "-" + day1, year2 + "-" + month2 + "-" + day2)

    if tle_download:
        if debug:
            print("requesting file from server")

        result = query.tle_query(epoch=date_range,norad_cat_id=noradid)

        with open(str(noradid) + ".txt", "w") as outfile:
            json.dump(result.json(), outfile)

        line1 = result.json()[0]["OBJECT_NAME"]
        line2 = result.json()[0]["TLE_LINE1"]
        line3 = result.json()[0]["TLE_LINE2"]

    elif path.exists(str(noradid) + ".txt") == False:

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


def getSatXY(line1, line2, line3, UT, mwa, ts, pixel_scale, imgSize, wcs):
    """
    determines of the satellite is within fov and if so returns t
    he x y pixel location of the satellite
    Paramters
    ---------
    line1   : line1 of tle
    line2   : line2 of tle
    line3   : line3 of tle
    UT      : UTC time
    mwa     : skyfield observer object
    ts      : skyfield time object
    Returns
    -------
    px                  : x location of satellite
    py                  : y location of satellite
    number_of_pixels    : search cone radius
    visible             : visible? (bool)
    """

    imgSize = 1400

    visible = False
    number_of_pixels = 0
    sat = EarthSatellite(line2, line3, line1, ts)

    # propogate satellite to utc
    time = ts.utc(UT.year, UT.month,  UT.day, UT.hour, UT.minute, UT.second)
    sat.at(time)
    difference = sat - mwa
    topocentric = difference.at(time)

    # determine angular and pixel space location of satellite
    ra, dec, distance = topocentric.radec()
    ra, dec = np.degrees(ra.radians), np.degrees(dec.radians)
    px, py = wcs.all_world2pix([[ra, dec]], 1)[0]

    ## check if satellite within image

    if 0 < px < imgSize and 0 < py < imgSize:

        if np.sqrt((px - (imgSize/2))**2 + (py -(imgSize/2))**2) < imgSize/3:
            LOS_range = distance.m

            radius = np.degrees(70000/LOS_range) # 25 km offset cone in orbit (can be made smaller) s=rxtheta
            number_of_pixels = radius/pixel_scale


    return px, py, number_of_pixels

def getDSNRS(px, py, radius, x, y, data):
    array_radius = np.sqrt((x - px)**2 + (y - py)**2)
    array_radius[array_radius > radius] = -10
    array_radius[array_radius != -10] = 1
    array_radius[array_radius != 1] = 0
    masked_data = array_radius * data
    masked_noise = data * np.abs(array_radius - 1)

    #numerator = np.sum(masked_data)
    numerator = np.max(masked_data)
    denominator = np.std(masked_noise)

    #dsnrs = numerator/ (denominator*np.sum(array_radius))
    dsnrs = numerator/ denominator

    return dsnrs


def main(args):

    timeSteps = np.arange(args.t1, args.t2 + 1)
    waterfall = np.zeros((len(timeSteps), 768))


    ts = load.timescale()
    mwa = Topos("26.703319405555554 S", "116.91558083333334 E", elevation_m= 377.827)

    trigger = False

    for i in timeSteps:

        if debug:
            print("working in timeStep {0}".format(i))

        for f in tqdm.tqdm(range(768)):

            hdu1 = fits.open(str(args.obs) + "-2m-"+str(i)+"-"+str(f).zfill(4)+"-dirty.fits")
            hdu2 = fits.open(str(args.obs) + "-2m-"+str(i+1)+"-"+str(f).zfill(4)+"-dirty.fits")
            utc = datetime.strptime(hdu2[0].header["DATE-OBS"], '%Y-%m-%dT%H:%M:%S.%f')
            data = hdu2[0].data[0,0,:,:] - hdu1[0].data[0,0,:,:]
            if not trigger:
                wcs = WCS(hdu2[0].header, naxis=2)
                imgSize = hdu2[0].header["NAXIS1"]
                pixel_scale = hdu2[0].header["CDELT2"]
                line1 , line2, line3 = obtainTLE(args.noradid, utc)
                x = np.linspace(0, (imgSize-1), imgSize)
                y = np.linspace(0, (imgSize-1), imgSize)
                x, y = np.meshgrid(x, y)
                trigger = True
            

            px, py, number_of_pixels = getSatXY(line1, line2, line3, utc, mwa, ts, pixel_scale, imgSize, wcs)
            dsnrs = getDSNRS(px, py, number_of_pixels, x, y, data)

            waterfall[i,f] = dsnrs

    np.save(str(args.obs) + "_" +str(args.noradid) + "waterfall.npy", waterfall)


            



if __name__ == "__main__":
    parser = ArgumentParser("advanced dsnrs", description="extracts all the info for orbit determination")
    parser.add_argument("--obs", required=True, type=int, help= "the observation id")
    parser.add_argument("--noradid", required=True, type=int, help="the noradid of the satellite")
    parser.add_argument("--user", required=True, help="the user name for space-track.org")
    parser.add_argument("--passwd", required=True, help="the password for the space-track.org")
    parser.add_argument("--debug", default=False, type=bool, help="run script in debug mode")
    parser.add_argument("--t1", required=True, type=int, help="the starting timeStep")
    parser.add_argument("--t2", required=True, type=int, help="the last timeStep")
    parser.add_argument("--force_tle_download", default=True, type=bool, help="download tle from server every time")
    args = parser.parse_args()

    global debug, query, tle_download
    tle_download = args.force_tle_download
    debug = args.debug
    query = st.SpaceTrackClient(args.user, args.passwd)
    
    if debug:
        print("running in debug mode")

    main(args)