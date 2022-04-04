#!/usr/bin/env python
from __future__ import division, print_function
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from datetime import datetime, timedelta
from argparse import ArgumentParser
from skyfield.api import EarthSatellite
from skyfield.api import Topos, load
from argparse import ArgumentParser
from tqdm import tqdm
import json
import csv


def dayOfYear(utc):
    """
    converts utc to day of year
    """
    time_delta = utc - datetime.strptime(str(utc.year)+"-01-01 00:00:00", '%Y-%m-%d %H:%M:%S')
    return time_delta.total_seconds()/(24.0*60*60) +1

def getTLE(args, start_utc):
    obs_doy = dayOfYear(start_utc)

    ## read catalog
    with open(args.inputTLEcatalog) as json_file:
        catalog = json.load(json_file)

    entries = len(catalog[:])
    print("{} entries found in input catalog".format(entries))

    norad_index_array = []
    ref_doy_array = []
    line1_array = []
    line2_array = []
    print("sorting through catalog to find tles")
    for i in tqdm(range(entries)):
        line1 = catalog[i]["TLE_LINE1"]
        norad = line1[2:7]
        if int(norad) == args.noradid:
            norad_index_array.append(i)
            doy = float(line1[20:32])
            ref_doy_array.append(np.abs(obs_doy - doy))
            line2 = catalog[i]["TLE_LINE2"]
            line1_array.append(line1)
            line2_array.append(line2)

    min_pos = np.where(ref_doy_array == np.min(ref_doy_array))[0][0]
    print("delta time bw obs and tle {}".format(ref_doy_array[min_pos]))
    
    return line1_array[min_pos], line2_array[min_pos]

def getCutOff(distance):
    ## using nea-field equation from the xiang's metoer paper
    wavelength = 3.059106714 ## lambda at 98MHz
    d = np.sqrt(wavelength*distance/2)
    return d

def main(args):

    ## load metafits
    hdu = fits.open(args.metafits)
    duration = hdu[0].header["EXPOSURE"]
    quack = 3 ## to overcome the hardcoded behaviour of cotter

    try:
        start_utc = datetime.strptime(hdu[0].header["DATE-OBS"], '%Y-%m-%dT%H:%M:%S.%f')
    except:
        start_utc = datetime.strptime(hdu[0].header["DATE-OBS"], '%Y-%m-%dT%H:%M:%S')
    
    start_utc =  start_utc + timedelta(seconds=quack) ## update start time for quack time

    line1, line2 = getTLE(args, start_utc)
    print("found TLEs\nline1\n{}\nline2\n{}".format(line1, line2))
    ts = load.timescale(builtin=True)
    mwa = Topos("26.701276778 S", "116.670846137 E", elevation_m= 377.827) 
    sat = EarthSatellite(line1, line2, "sat", ts)
    pointing_ra = float(hdu[0].header["RA"])
    pointing_dec = float(hdu[0].header["DEC"])

    ## create delta array (in seconds)
    search_timeSteps, search_ra, search_dec = [], [], [] ## contains the timesteps with satellite signal
    baseline_cutoff = []
    time_array = []
    for timeStep in range(int(duration/args.integration)):
        ### skip the first 2 timeSteps and the last timeStep
        if timeStep in [0, 1, int(duration/args.integration)-1]:
            continue

        local_utc = start_utc + timedelta(seconds=timeStep*args.integration)

        local_ts = ts.utc(local_utc.year, local_utc.month, local_utc.day, local_utc.hour, local_utc.minute, local_utc.second + local_utc.microsecond/1000000.0)

        sat.at(local_ts)
        difference = sat - mwa
        topocentric = difference.at(local_ts)
        ra, dec, distance = topocentric.radec()
        ra_deg, dec_deg = np.degrees(ra.radians), np.degrees(dec.radians)
        dist = np.sqrt((ra_deg- pointing_ra)**2 + (dec_deg - pointing_dec)**2)
        cutOff = getCutOff(distance.m)
        print("ra {} dec {} timeStep {} dist {} los dist {} cutoff {}".format(ra, dec, timeStep, dist, distance.m, cutOff))

        ## update ra and dec in chgcentre format
        ra = str(ra).replace(" ", "")
        dec = str(dec).replace(" ", "").replace('"', 's').replace("'", "m").replace("deg", "d")
        search_timeSteps.append(timeStep)
        search_ra.append(ra)
        search_dec.append(dec)
        baseline_cutoff.append(cutOff)
        time_array.append(str(local_utc))

    with open("lightCurve" + str(args.obs) + "-" + str(args.noradid)+".csv", "w") as vsc:
        thewriter = csv.writer(vsc)
        for t, ra, dec, c, ut in zip(search_timeSteps,search_ra, search_dec, baseline_cutoff, time_array):
            ut1, ut2 = str(ut).split(" ")
            ut = ut1 + "T" + ut2
            line = [t, ra, dec, c, ut, "dummy"]
            thewriter.writerow(line)

    


if __name__ == "__main__":
    parser = ArgumentParser("light_curve_trak", description="finds time-steps where satellite is within horizon")
    parser.add_argument("--obs", required=True, type=int, help="the obs id")
    parser.add_argument("--metafits", required=True, help="the observation metafits file")
    parser.add_argument("--noradid", required=True, type=int,help="the noradid of the satellite")
    parser.add_argument("--inputTLEcatalog", required=True, help="the noradid of the satellite")
    parser.add_argument("--integration", default=2, type=float, help="the image integration time")
    args = parser.parse_args()
    
    main(args)