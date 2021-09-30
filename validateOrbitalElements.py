#!/usr/bin/env python
from __future__ import print_function, division
from argparse import ArgumentParser
import numpy as np
from datetime import datetime, timedelta
from skyfield.api import Topos, load, EarthSatellite
from astropy.io import fits
from astropy.wcs import WCS
import yaml
from tqdm import tqdm
import spacetracktool as st
from spacetracktool import operations as ops
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def dayOfYear(utc):
    """
    converts utc to day of year
    """
    time_delta = utc - datetime.strptime(str(utc.year)+"-01-01 00:00:00", '%Y-%m-%d %H:%M:%S')
    return time_delta.total_seconds()/(24.0*60*60) +1
    

def tle_query(utc, args):
    """
    returns the historic tle for the object
    """
    query = st.SpaceTrackClient(args.user, args.passwd)
    start_date = utc + timedelta(hours=-24*14)
    end_date = utc 
    date_range = ops.make_range_string(str(start_date.year) + "-" + str(start_date.month).zfill(2) + "-"+ str(start_date.day).zfill(2),
        str(end_date.year) + "-" + str(end_date.month).zfill(2) + "-"+ str(end_date.day).zfill(2))
    result = query.tle_query(epoch=date_range, norad_cat_id=args.norad)
    return result


def loadSolution(fileName):
    f = open(fileName)
    line = f.readline()
    line = f.readline()
        
    i = line.split(" ")[1]
    ra = line.split(" ")[3]
    e = line.split(" ")[5]
    aop = line.split(" ")[7]
    ma = line.split(" ")[9]
    mm = line.split(" ")[11][:-2]

    line = f.readline()
    line = f.readline()

    i_err = line.split(" ")[1]
    ra_err = line.split(" ")[3]
    e_err = line.split(" ")[5]
    aop_err = line.split(" ")[7]
    ma_err = line.split(" ")[9]
    mm_err = line.split(" ")[11][:-2]

    return i, ra, e, aop, ma, mm, i_err, ra_err, e_err, aop_err, ma_err, mm_err


def getHistoricTLE(result):
    ## get historic data
    e_array, i_array, ra_array = [], [], []
    aop_array, ma_array, mm_array = [], [], []
    doy_array = []
    

    for i in tqdm(range(len(result.json()[:]))):
        line1 = result.json()[i]["TLE_LINE1"]
        line2 = result.json()[i]["TLE_LINE2"]
        utc = datetime.strptime(result.json()[i]["EPOCH"], '%Y-%m-%d %H:%M:%S')
        days = dayOfYear(utc)
        e = float(line2[26:33])/10000000
        i =  float(line2[8:16])
        ra = float(line2[17:25])
        aop = float(line2[34:42])
        ma = float(line2[43:51])
        mm = float(line2[52:63])
        doy = float(line1[20:32])
        e_array.append(e)
        i_array.append(i)
        ra_array.append(ra)
        aop_array.append(aop)
        ma_array.append(ma)
        mm_array.append(mm)
        doy_array.append(doy)

    return i_array, ra_array, e_array, aop_array, ma_array, mm_array, doy_array


        


def main(args):

    ## obtain the wcs object and the datetime of the observation
    hdu = fits.open(args.wcsFile)
    obs_date = datetime.strptime(hdu[0].header["DATE-OBS"],'%Y-%m-%dT%H:%M:%S.%f')
    obs_dayofyear = dayOfYear(obs_date)

    ## load orbital elements
    i, ra, e, aop, ma, mm, i_err, ra_err, e_err, aop_err, ma_err, mm_err = loadSolution(args.solutionFile)

    ## perform the tle query
    result = tle_query(obs_date,  args)

    i_array, ra_array, e_array, aop_array, ma_array, mm_array, doy_array = getHistoricTLE(result)

    fig, ax = plt.subplots(figsize=(15, 15))

    plt.subplot(321)
    plt.scatter(doy_array, i_array, label="historic TLE", color="black")
    plt.errorbar(float(obs_dayofyear), float(i), yerr=float(i_err), label="solution", barsabove=True, capsize=10)
    plt.legend()
    plt.xlabel("Day of year")
    plt.ylabel("inclination (deg)")

    plt.subplot(322)
    plt.scatter(doy_array, ra_array, label="historic TLE", color="black")
    plt.errorbar(float(obs_dayofyear), float(ra), yerr=float(ra_err), label="solution", barsabove=True, capsize=10)
    plt.legend()
    plt.xlabel("Day of year")
    plt.ylabel("raan (deg)")

    plt.subplot(323)
    plt.scatter(doy_array, e_array, label="historic TLE", color="black")
    plt.errorbar(float(obs_dayofyear), float(e), yerr=float(e_err), label="solution", barsabove=True, capsize=10)
    plt.legend()
    plt.xlabel("Day of year")
    plt.ylabel("eccentricity")

    plt.subplot(324)
    plt.scatter(doy_array, aop_array, label="historic TLE", color="black")
    plt.errorbar(float(obs_dayofyear), float(aop), yerr=float(aop_err), label="solution", barsabove=True, capsize=10)
    plt.legend()
    plt.xlabel("Day of year")
    plt.ylabel("aop (deg)")

    plt.subplot(325)
    plt.scatter(doy_array, ma_array, label="historic TLE", color="black")
    plt.errorbar(float(obs_dayofyear), float(ma), yerr=float(ma_err), label="solution", barsabove=True, capsize=10)
    plt.legend()
    plt.xlabel("Day of year")
    plt.ylabel("ma (deg)")

    plt.subplot(326)
    plt.scatter(doy_array, mm_array, label="historic TLE", color="black")
    plt.errorbar(float(obs_dayofyear), float(mm), yerr=float(mm_err), label="solution", barsabove=True, capsize=10)
    plt.legend()
    plt.xlabel("Day of year")
    plt.ylabel("mm (rev/day)")

    plt.savefig("validationOrbitalElements{}.png".format(args.norad))










if __name__ == "__main__":
    parser = ArgumentParser("validate orbital elements", description="validates the converted orbital elements agains the hisotoric tle")
    parser.add_argument("--norad", type=int, required=True, help="the norad id")
    parser.add_argument("--wcsFile", required=True, help="a dummy file that has the required wcs object")
    parser.add_argument("--user", required=True, help="the user name for space-track.org")
    parser.add_argument("--passwd", required=True, help="the password for space-track.org")
    parser.add_argument("--solutionFile", required=True, help="the orbital element solution file")
    args =  parser.parse_args()

    main(args)






